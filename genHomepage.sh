#!/usr/bin/env bash

#paths
FORMAT='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/signatureFormat.txt'
WEB_SCRIPT='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/generateTDFrontpage.py'
BASE_SITE='/home/mike/workspace/PellegriniResearch/frontend/makeCheckboxes/iSigPreCheckboxes.html'
OUTPUT_SITE='/home/mike/workspace/PellegriniResearch/frontend/iSigHomepage.html'
echo "Generating the new homepage at $OUTPUT_SITE using the base at $BASE_SITE and the format file $FORMAT"
python $WEB_SCRIPT -f $FORMAT -b $BASE_SITE -o $OUTPUT_SITE
echo -e "Done generating website"
