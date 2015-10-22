#!/usr/bin/perl -w
use strict;
use List::Util qw(max min);
use File::Basename;
use File::Copy;
use Data::Dumper;
use Getopt::Long;
use JSON;
use GD;


my $BAM2FQ_BIN   = "bamToFastq";
my $SAMTOOLS_BIN = "samtools";
my $BWA_BIN      = "bwa";
my $USE_BWA      = 1;


# Parse command line options
my %opt = get_options( );

# Get plugin settings
my $settings = read_json( $opt{'output-dir'}."/startplugin.json" );
my $samples =  read_json( $opt{'output-dir'}."/barcodes.json" );

my $basecaller_dir = $settings->{runinfo}->{basecaller_dir};


# Collect necessary sample information
my( %ref_files, %sample_data );
foreach my $id (keys %$samples) {

    # FIXME!!! Select samples/references somehow
    next unless $samples->{$id}->{sample} =~ /(amplicon|amplikon)/;

    my $ref = $samples->{$id}->{reference};
    my $ref_path = '/results/referenceLibrary/tmap-f3/' . $ref . '/' . $ref . '.fasta';
    my $bam = $samples->{$id}->{bam};
    $bam =~ s/_rawlib\./_rawlib\.basecaller\./;
    my $bam_path = $basecaller_dir."/".$bam;

    unless (-s $bam_path) { 
	print "WARNING: Could not find bam file for $id, omitting...\n";
	next;
    }
    unless (-s $ref_path) { 
	print "WARNING: Could not find reference file for $ref, omitting...\n";
	next;
    }

    $sample_data{$id} = {'bam' => $bam_path,
			 'ref' => $ref_path,
			 'name'=> $samples->{$id}->{sample} };
    
    $ref_files{ $ref_path } = 1;

    # If BWA, copy reference fasta to output dir and build index, unless index has already
    # been built manuallt. (Necessary due to lack of writing permissions in ref-path)
    if ($USE_BWA and !-s "$ref_path.bwt") {
	my $ref_copy = $opt{'output-dir'} . "/" . basename( $ref_path );

	# Only build each reference once
	unless ( -s "$ref_copy.bwt" ) {
	    copy( $ref_path, $ref_copy );
	    system( $BWA_BIN." index ".$ref_copy );
	}

	# Change reference path for this sample
	$sample_data{$id} -> {ref} = $ref_copy;
    }
}


# Get amplicon sequence lengths
my %len = get_lengths( keys %ref_files );


# Process raw sequence data 
my( %mapped_reads_per_amplicon, %coverage_by_position );
foreach my $id (sort keys %sample_data) {
    print "Processing $id...";
    
    my $in_bam = $sample_data{$id}->{bam};
    my $ref    = $sample_data{$id}->{ref};

    my $filebase = $opt{'output-dir'} . "/" . basename( $in_bam );
    $filebase =~ s/\.bam$//;

    my $in_fq       = $filebase.".fq";
    my $out_sam     = $filebase.".map.sam";
    my $out_bam     = $filebase.".map.bam";
    my $out_sortbam = $filebase.".map.sort.bam";
    my $out_cas     = $filebase.".map.cas";

    # Convert unmapped BAM to FastQ
    system( $BAM2FQ_BIN. " -i $in_bam -fq $in_fq" );

    # Map reads to reference sequences, and create BAM output
    if ($USE_BWA) {
	system( $BWA_BIN." mem -t 4 $ref $in_fq > $out_sam 2> bwa.log");
	system( $SAMTOOLS_BIN." view -Sb $out_sam  > $out_bam 2> /dev/null");
    }
    else {
	system( "clc_mapper_beta -d $ref -q $in_fq -g 1 -e 1 -G 2 -E 2 -s 0.6 --cpus 4 -o $out_cas > /dev/null");
	system( "clc_cas_to_sam -a $out_cas -u -o $out_bam > /dev/null" );
    }

    # Sort BAM file
    system( $SAMTOOLS_BIN. " sort -f $out_bam $out_sortbam" );

    # Get coverage data from sorted bam file
    $mapped_reads_per_amplicon{$id} = count_per_amplicon( $out_sortbam );
    $coverage_by_position{$id}      = get_coverage_by_position( $out_sortbam );

    # Clean up intermediate files
    unlink( $in_fq );
    unlink( $out_cas ) unless $USE_BWA;
    unlink( $out_sam ) if $USE_BWA;
    unlink( $out_bam );
    unlink( $in_fq );
    print " done.\n";
}    


print "Creating output.\n";
# Output average coverage table
print "seq";
print "\t".$_ foreach ( sort keys %sample_data );
print "\n";
my %avg_cov_per_amplicon;
foreach my $seq_id ( sort keys %len ) {
    if( length $seq_id > 30 ) {
	print substr( $seq_id, 0, 30 ) . "...";
    }
    else {
	print $seq_id;
    }
    
    foreach my $sample_id ( sort keys %sample_data ) {
	my $sum = 0;
	foreach my $pos (keys %{$coverage_by_position{$sample_id}->{$seq_id}}) {
	    $sum += ( $coverage_by_position{$sample_id}->{$seq_id}->{$pos} or 0 );
	}
	my $avg_cov = sprintf "%.2f", $sum/$len{$seq_id};
	print "\t$avg_cov";
	$avg_cov_per_amplicon{$sample_id}->{$seq_id} = $avg_cov
    }
    print "\n";
}
print "\n";


# Output read count table
print "seq";
print "\t".$_ foreach ( sort keys %sample_data );
print "\n";
foreach my $seqid ( sort keys %len ) {
    $seqid =~ /^(.*?):/;
    print $1;
    foreach my $sample_id ( sort keys %sample_data ) {
	printf "\t%d", ( $mapped_reads_per_amplicon{$sample_id}->{$seqid} or 0 );
    }
    print "\n";
}


# Make read count table in "block" html. 
open( HTML_BLOCK, ">ampliconCoverage_block.html" );

print HTML_BLOCK '<html><head></head><body>
<style>
    *{font-family:Arial; font-size:8pt}
    .alnright { text-align: right; }
    .fail { background-color: #FF6666; }
    .low { background-color: #FFFF66; }
    .header{ background-color: #BBBBCC; font-weight:bold; }
    .odd   { background-color: #FFFFFF; }
    .even  { background-color: #DDDDDD; }
    td { padding:4px; }
    table {
      border-collapse: collapse;
    }
    table, td, tr { border: 1px solid #707070 }
</style>
';



# Create read count HTML table
print HTML_BLOCK "<b>Number of mapped reads per amplicon:</b><br>";
print HTML_BLOCK "<table><tr class='header'><td>amplicon</td>";
print HTML_BLOCK "<td>".$_."</td>" foreach ( sort keys %sample_data );
print HTML_BLOCK "</tr>";

my $row = 0;
foreach my $seqid ( sort keys %len ) {
    $row++;

    my $show_id = $seqid;
    if( length $seqid > 30 ) {
        $show_id = substr( $seqid, 0, 30 ) . "...";
    }

    $seqid =~ /^(.*?):/;
    print HTML_BLOCK "<tr class='".($row%2==0?"even":"odd")."'><td><span title='$seqid'>$show_id</span></td>";
    foreach my $sample_id ( sort keys %sample_data ) {
	printf HTML_BLOCK "<td class='alnright'>%d</td>", ( $mapped_reads_per_amplicon{$sample_id}->{$seqid} or 0 );
    }
    print HTML_BLOCK "</tr>";
}
print HTML_BLOCK "</table>";


# Create average coverage HTML table
print HTML_BLOCK "<p><b>Average coverage (x) per amplicon:</b><br>";
print HTML_BLOCK "<table><tr class='header'><td>amplicon</td>";
print HTML_BLOCK "<td>".$_."</td>" foreach ( sort keys %sample_data );
print HTML_BLOCK "</tr>";

$row = 0;
foreach my $seqid ( sort keys %len ) {
    $row++;

    my $show_id = $seqid;
    if( length $seqid > 30 ) {
        $show_id = substr( $seqid, 0, 30 ) . "...";
    }


    print HTML_BLOCK "<tr class='".($row%2==0?"even":"odd")."'><td><span title='$seqid'>$show_id</span></td>";
    foreach my $sample_id ( sort keys %sample_data ) {
	my $cov = ( $avg_cov_per_amplicon{$sample_id}->{$seqid} or 0 );
	printf HTML_BLOCK "<td class='alnright %s'>%.2f</td>", ( $cov==0 ? "fail" : ($cov<10 ? "low" : "") ), $cov;
    }
    print HTML_BLOCK "</tr>";
}
print HTML_BLOCK "</table>";
close HTML_BLOCK;




# Make mapped coverage plot images
my %plot_names = create_coverage_plots_scale( \%coverage_by_position, \%len, \%sample_data );

# Output coverage profiles in HTML file
open( HTML, ">ampliconCoverage_results.html" );

foreach my $samp_id (sort keys %plot_names) {
    print HTML "<b>$samp_id (".$sample_data{$samp_id}->{name}.")</b><p>";
    foreach my $seq_id ( sort keys %{$plot_names{$samp_id}} ) {
	print HTML "<img src='".$plot_names{$samp_id}->{$seq_id}."'><br>";
    }
    print HTML "<hr>";
}
close HTML;


#####################################################


sub get_lengths {
    my @files = @_;

    my %lens;
    foreach my $fn (@files) {

	open( FA, $fn );
	while (<FA>) {
	    if (/^>(.+)$/) {
		chomp;
		my $id = $1;
		chomp (my $seq = <FA>);
		my $len = length($seq);
		$lens{$id} = $len;
	    }
	}
	close FA;
    }
    return %lens;
}


sub count_per_amplicon {
    my $bam = shift;

    my %cnt;
    open( BAM, $SAMTOOLS_BIN." view $bam|" );
    while (<BAM>) {
	next if /^@/;
	chomp;
	my @bam = split /\t/;
	$cnt{$bam[2]}++;
    }
    close BAM;

    return \%cnt;
}



sub get_coverage_by_position {
    my $bam = shift;

    my %cov;
    open( PILE, $SAMTOOLS_BIN." mpileup -d 100000 $bam 2> /dev/null|" );
    while (<PILE>) {
	my @pile = split /\t/;
	$cov{$pile[0]}->{$pile[1]} = $pile[3];
    }
    close PILE;
    
    return \%cov;
}


sub get_options {
    my %opt;
    GetOptions( \%opt, 'install-dir=s', 'output-dir=s', 'output-url=s', 'report-dir=s' );
    die "No output-dir specified..." unless $opt{'output-dir'};
    return %opt;
}



sub read_json {
    my $json_fn = shift;

    local $/;
    open( my $fh, '<', $json_fn );
    my $json_text = <$fh>;
    my $decoded   = decode_json( $json_text );
    return $decoded;
}



sub create_coverage_plots_scale {
    my( $cov, $len, $samp ) = @_;

    mkdir "images";

    my %plot_files;
    my $MAX_PLOT_WIDTH = 800;
    my $PLOT_HEIGHT = 75;
    my $OFS = 20;
    my $IMG_HEIGHT = $PLOT_HEIGHT + $OFS * 2;

    foreach my $samp_id ( sort keys %$samp ) {
	foreach my $seq_id ( sort keys %$len ) {
	    my $seq_len = $len->{$seq_id};
	    my $plot_width = min( $seq_len, $MAX_PLOT_WIDTH ); 
	    my $x_scale = $plot_width / $seq_len;
	    my ($short_seq_id) = ($seq_id =~ /^(.*?):/);
	    my $max_cov = max( values %{$cov->{$samp_id}->{$seq_id}}, 20 );
	    my $IMG_WIDTH = $plot_width + $OFS * 2;

	    my $im = new GD::Image( $IMG_WIDTH, $IMG_HEIGHT, 1 );
	    my $white        = $im->colorAllocate(255,255,255);
	    my $gray         = $im->colorAllocate(100,100,110);
	    my $plot_bg_col  = $im->colorAllocate(235,235,255);
	    my $black        = $im->colorAllocate(0,0,0);       
	    my $curve_col    = $im->colorAllocate(50,80,200);
	    $im->filledRectangle( 0, 0, $IMG_WIDTH, $IMG_HEIGHT, $white );
	    $im->filledRectangle( $OFS, $OFS, $IMG_WIDTH - $OFS, $IMG_HEIGHT - $OFS, $plot_bg_col );
	    $im->setAntiAliased( $curve_col );

	    my $prev_y;
	    my( $pos_plot, $prev_pos_plot ) = ( 0,0 );
	    my @covs;
	    foreach my $pos_seq (1..$seq_len) {
		$pos_plot += $x_scale;
		my $cov = ( $cov->{$samp_id}->{$seq_id}->{$pos_seq} or 0 );
		push @covs, $cov;

		# If new x pixel reached, draw line
		if (int $pos_plot > int $prev_pos_plot) {
		    my $y = ($PLOT_HEIGHT+$OFS) - (avg(@covs)/$max_cov)*$PLOT_HEIGHT;
		    my $x = ($pos_plot-1) + $OFS;

		    unless ($prev_y) { $prev_y = $y; } # First round

		    $im->line( $x, $prev_y, $x+1, $y, gdAntiAliased );
		    $prev_y = $y;
		    @covs = ();
		}

		$prev_pos_plot = $pos_plot;
	    }


	    # Draw axes
	    $im->line( $OFS, $OFS, $OFS, $IMG_HEIGHT - $OFS, $black );
	    $im->line( $OFS, $IMG_HEIGHT - $OFS, $IMG_WIDTH - $OFS, $IMG_HEIGHT - $OFS, $black );

	    # Draw X ticks and labels
	    my $x_tick_step = 1;
	    $x_tick_step *= 10 until ( $seq_len / $x_tick_step < 10 );
	    for my $x_tick ( 0..int($seq_len / $x_tick_step) ) {
		my $x = $x_tick * $x_tick_step;
		$im->line( $OFS + $x*$x_scale, $IMG_HEIGHT - $OFS, $OFS + $x*$x_scale, $IMG_HEIGHT - $OFS + 3, $black );
		$im->string( gdTinyFont, $OFS + $x*$x_scale - (length($x)*3) + 1, $IMG_HEIGHT - $OFS + 5, $x, $black );
	    }

	    # Show max value of Y axis
	    $im->string( gdTinyFont, $OFS - length($max_cov) * 6, $OFS - 2, $max_cov, $black );

	    # Print amplicon sequence name
	    $im->string( gdTinyFont, $OFS, $OFS - 10, $seq_id, $gray );

	    # Output png
	    my $out_png = "images/".$samp_id.".".$short_seq_id.".png";
	    open( PNG, ">".$out_png );
	    print PNG $im->png;
	    close PNG;
	    
	    $plot_names{$samp_id}->{$seq_id} = $out_png;
	}
    }
    return %plot_names;
}


sub avg {
    my $sum;
    $sum+=$_ foreach @_;
    return $sum / @_;
}
