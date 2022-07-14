#!/bin/bash -ue
mkdir peptide_level psm_level
Rscript /Users/chuang8/proteomics/DeepRescore/bin/got_pga_input.R features.txt comet ./d2-rawPSMs.txt
Rscript /Users/chuang8/proteomics/DeepRescore/bin/calculate_fdr.R ./ d2 /Users/chuang8/proteomics/DeepRescore/bin/protein.pro-ref.fasta
