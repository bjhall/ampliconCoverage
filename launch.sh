#!/bin/bash
VERSION="0.0.4"

echo 'running launch.sh'
${DIRNAME}/amplicon_coverage_plugin.pl \
  --install-dir   ${DIRNAME} \
  --output-dir    ${TSP_FILEPATH_PLUGIN_DIR} \
  --output-url    ${TSP_URLPATH_PLUGIN_DIR} \
  --report-dir    ${ANALYSIS_DIR}
