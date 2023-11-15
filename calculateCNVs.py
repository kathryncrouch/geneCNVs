#!/usr/bin/env python3
# testing GitHub
from TPMtool.TPMtool import GffParser
from TPMtool.TPMtool import HtSeqCounts
from TPMtool.TPMtool import TPM
from CNVs.CNVs import CNV
import argparse

def __main__():

    parser = argparse.ArgumentParser(description='Estimate chromosome and gene copy number from read counts', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--genome', required=True, help='Reference genome annotation in GFF3 format.' )
    parser.add_argument('--input', required=True, help='Read counts per gene or transcript in htseq-count format' )
    parser.add_argument('--orthologFile', help='File containing gene to ortholog mappings' )
    parser.add_argument('--outputPrefix', required=True, help='Prefix for all output files' )
    parser.add_argument('--ploidy', default=2, type=int, help='Canonical ploidy for this organism')
    args = parser.parse_args()

    # hack to use TPM module as this is for RNAseq and expects a stranded arg
    args.stranded = False

    GeneModels = GffParser(args.genome)
    GeneCounts = HtSeqCounts()
    GeneCounts.senseCounts = GeneCounts.getCounts(args.input)
    Tpm = TPM(GeneModels, GeneCounts)
    Tpm.doTPMCalculation()
    Tpm.writeTPM(args)

    cnvs = CNV(GeneModels, Tpm, args.orthologFile, args.ploidy)
    cnvs.ploidyPrinter(args.outputPrefix)
    cnvs.cnvPrinter(args.outputPrefix)

    exit()

if __name__=="__main__": __main__()
