library(tidyverse)
library(PGA)

args <- commandArgs(T)
psm <- args[1]
out = args[2]
database <- args[3]
decoyP <- args[4]

calculateFDR(psmfile=psm,
	db=database,
	fdr=0.01,
	decoyPrefix=decoyP,
	better_score_lower=FALSE,
	remap=FALSE,
	peptide_level=TRUE,
	score_t = 0,
	protein_inference=FALSE,
	out_dir=paste(out,"_peptide_level", sep=""),
	xmx=20)

calculateFDR(psmfile=psm,
        db=database,
        fdr=0.01,
        decoyPrefix=decoyP,
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=FALSE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(out, "_psm_level", sep=""),
        xmx=20)
