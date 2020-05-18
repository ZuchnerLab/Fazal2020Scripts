#!/bin/bash
#after running eh-denovo on all the files, combine them together
EHDn-v0.6.2_HelperScripts/combine_counts.py --manifest manifest.txt --combinedCounts EHDnResults.combinedCounts.json

# use the compare script to convert to bed format
EHDn-v0.6.2_HelperScripts/compare_anchored_irrs.py --manifest manifest.txt --inputCounts EHDnResults.combinedCounts.json --outputRegions EHDnResults.combinedCounts.bed --minCount 5

# filter out the non-canonical chromosomes
grep -wv '^hs37d5' EHDnResults.combinedCounts.bed | grep -v '^GL' > EHDnResults.combinedCounts.filtered.bed

# reformat the filtered bed file into a matrix of anchored IRR counts per sample per site
Rscript ./reformatEHDenovoOutput.R
