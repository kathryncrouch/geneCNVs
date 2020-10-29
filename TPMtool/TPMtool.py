#!/usr/bin/env python3

import argparse
import sys
import os, shutil, subprocess, sys, tempfile

#####Class to parse GFF for gene lengths and exon/transcript/parent relationships#####
class GffParser(object):

    def __init__(self, gff):
        self.featureTypeCol = 2
        self.startCol = 3
        self.endCol = 4
        self.attributeCol = 8
        self.gff = gff
        self.transcriptParents = self._getTranscriptParents()
        self.geneLengths, self.geneChromosomes = self._getGeneLengths()
                   

    #####Private methods called by constructor or other class methods#####
    def _getTranscriptParents(self):
        transcriptParents = {}
        try:
            file = open(self.gff)
            for line in file:
                line = line.rstrip()

                if line.startswith('##'):
                    continue

                fields = line.split('\t')

                if fields[self.featureTypeCol] == 'transcript' or fields[self.featureTypeCol].endswith('RNA'):
                    id = None;
                    parent = None;
                    attributeString = fields[self.attributeCol]

                    for attr in attributeString.split(';'):
                        if attr.startswith('ID='):
                            id = attr[3:]
                        if attr.startswith('Parent='):
                            parent = attr[7:]

                    if not id:
                        raise GffFormatError('No ID found in transcript attribute string  "' + attributeString + '"\n')
                    if not parent:
                        raise GffFormatError('No parent found in transcript attribute string "' + attributeString + '"\n')

                    transcriptParents[id] = parent

            file.close()
        except IOError:
            raise SystemExit('File ' + self.gff + ' could not be opened for reading.\n')
        except GffFormatError as e:
            raise SystemExit(e)

        return transcriptParents


    def _getGeneLengths(self):
        geneLengths = {}
        geneChromosomes = {}

        try:
            file = open(self.gff)
            for line in file:
                line = line.rstrip()

                if line.startswith('##'):
                    continue
                
                fields = line.split('\t')

                if fields[self.featureTypeCol] == 'exon':
                    try:
                        length = int(fields[self.endCol]) - int(fields[self.startCol]) + 1
                        length = length/1000
                    except ValueError:
                        raise SystemExit('Start value ' + fields[self.startCol] + ' or end value '  + fields[self.endCol] + ' could not be converted to a numeric type.\n')

                    attributeString = fields[self.attributeCol]
                    chromosome = fields[0]

                    for attr in attributeString.split(';'):
                        if attr.startswith('Parent='):
                            parent = attr[7:]

                            # if we get a comma-separated list of transcripts as parents, use the first one
                            if ',' in parent:
                                parent = parent[:parent.find(',')]

                            if not parent in self.transcriptParents:
                                raise GffFormatError('No parent gene found for transcript "' + parent + '"\n')

                            gene = self.transcriptParents[parent]

                            if not gene in geneLengths:
                                geneLengths[gene] = 0

                            geneLengths[gene] = geneLengths[gene] + length
                            if gene in geneChromosomes:
                                if geneChromosomes[gene] != chromosome:
                                    raise SystemExit("Gene {} appears to have exons on different chromosomes. Please check the GFF file.\n".format(gene))
                            else:
                                geneChromosomes[gene] = chromosome
        
            file.close()
        except IOError:
            raise SystemExit('File ' + self.gff + ' could not be opened for reading.\n')
        except GffFormatError as e:
            raise SystemExit(e)

        return (geneLengths, geneChromosomes)



#####Class to parse htseq-count output for gene counts#####
class HtSeqCounts(object):

    def getCounts (self, input):
        specialCounters = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']
        geneCounts = {}
        try:
            file = open(input)
            for line in file:
                line = line.rstrip()
                if line.startswith('#ID'):
                    continue
                
                fields = line.split('\t')
                id = fields[0]

                # if we get a comma-separated list of transcripts as parents, use the first one
                if ',' in id:
                    id = id[:id.find(',')]

                if id in specialCounters:
                    continue
                
                try:
                    count = int(fields[1])
                except ValueError:
                    raise SystemExit('Count value ' + fields[1] + ' could not be converted to a numeric type.\n')

                geneCounts[id] = count
            
            file.close()
        except IOError:
            raise SystemExit('File ' + self.countFile + ' could not be opened for reading.\n')

        return geneCounts



#####class to calculate TPM#####
class TPM(object):

    def __init__(self, geneModels, geneCounts, stranded=False):
        self.geneModels = geneModels
        self.geneCounts = geneCounts
        self.stranded = stranded


    def doTPMCalculation(self):
        self.rpkSum, senseRpkHash = self._calcRPK(self.geneCounts.senseCounts)
        if self.stranded:
            self.rpkSum, antisenseRpkHash = self._calcRPK(self.geneCounts.antisenseCounts, self.rpkSum)
            self.antisenseTPM = self._calcTPM(antisenseRpkHash)
        else: 
            antisenseTPM = None
        #note: rpk must be caluculated for both sense and antisense before tpm is calculated.
        self.senseTPM = self._calcTPM(senseRpkHash)


    def writeTPM(self, args):
        if self.senseTPM:
            self._writeTPMFile(self.senseTPM, args.outputPrefix + ".tpm")
        else:
            raise SystemExit('TPM values have not been calculated yet. Please use TPM.doTPMCalculation() to calulate TPM values before writing to file.\n')

        if self.stranded and self.antisenseTPM:
            self._writeTPMFile(self.antisenseTPM, args.antisense_outputPrefix + ".tpm")
        elif args.stranded and not self.antisesnseTPM:
            raise SystemExit('TPM values have not been calculated yet. Please use TPM.doTPMCalculation() to calulate TPM values before writing to file.\n')


    
    #####Private methods called by constructor or other class methods#####
    def _writeTPMFile(self, tpmHash, file):
        try:
            outFile = open(file, 'w')
            outFile.write('gene_id\tTPM\n')
            for gene, tpm in tpmHash.items():
                outFile.write("%s\t%.12f\n" % (gene, tpm))
            outFile.close()
        except IOError:
            raise  SystemError('Output file ' + outFile + ' cannot be opened for writing.\n')
                

    def _calcRPK(self, countHash, rpkSum=0):
        rpkHash = {}

        for gene, count in countHash.items():
            length = None

            if gene in self.geneModels.geneLengths:
                length = self.geneModels.geneLengths[gene]
            elif gene in self.geneModels.transcriptParents:
                gene = self.geneModels.transcriptParents[gene]
                length = self.geneModels.geneLengths[gene]
            else:
                raise SystemExit('ID "' + gene + '" in counts file not recognized as gene or transcript ID')
        
            try:
                rpk = count / length
            except ZeroDivisionError as e:
                raise SystemExit(e + '\nGene length for gene ' + gene + ' is ' + str(length) + '. Gene length should not be 0.\n')

            rpkSum += rpk
            rpkHash[gene] = rpk

        return rpkSum, rpkHash


    def  _calcTPM(self, rpkHash):
        tpmHash = {}
        for gene, rpk in rpkHash.items():
            tpm = rpk/(self.rpkSum/1000000)
            tpmHash[gene] = tpm

        return tpmHash

            


#####Custom exception classes#####
class GffFormatError(Exception):

    def __init__(self, data):
        self.data = data

    def __str__(self):
        return self.data


#####main#####
def __main__():

    # Parse Command Line
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', required=True, help='Reference genome annotation in GFF3 format.' )
    parser.add_argument('--input', required=True, help='Read counts per gene or transcript in htseq-count format (sense strand if stranded).' )
    parser.add_argument('--output', required=True, help='File to write TPM counts (sense strand if stranded)')
    parser.add_argument('--stranded', action='store_true', help='Use this flag if your dataset is stranded and you want TPM data for both strands')
    parser.add_argument('--antisense_input', help='Antisense read counts per gene or transcript in htseq-count format. Must be declared if --stranded is used.' )
    parser.add_argument('--antisense_output', help='File to write antisense TPM counts. Must be declared if --stranded is used.')
    args = parser.parse_args()

    if args.stranded and (not args.antisense_input or not args.antisense_output):
        parser.print_help()
        raise SystemExit('\nERROR: If the dataset is double stranded, an input and output files for antisense data must be given at the command line.\n')

    if args.antisense_input and args.input == args.antisense_input:
        raise SystemExit('ERROR: Output file names for sense and antisense strands are identical. Please specify different names for these two output files.\n')


    GeneModels = GffParser(args.genome)
    GeneCounts = HtSeqCounts()
    GeneCounts.senseCounts = GeneCounts.getCounts(args.input)
    if args.stranded:
        GeneCounts.antisenseCounts = GeneCounts.getCounts(args.antisense_input)

    Tpm = TPM(GeneModels, GeneCounts, args.stranded)
    Tpm.doTPMCalculation()
    Tpm.writeTPM(args)
    exit()

if __name__=="__main__": __main__()
