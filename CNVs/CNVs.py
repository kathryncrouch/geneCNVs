#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
import statistics
import re

class CNV(object):

    def __init__(self, GeneModels, TPM, orthologFile, basePloidy):
        self.TPM = TPM
        self.GeneModels = GeneModels
        self.orthologFile = orthologFile
        self.basePloidy = basePloidy
        self.tpmChrMedians, self.tpmGenomeMedian = self._calculateTpmMedians()
        self.chrPloidy = self._calculatePloidy()
        self.geneCopyNumbers = self._calculateGeneCopyNumbers()
        self.annotations = self._parseOrthologs()
        self._calculateCNVForOrthologGroups()




    def _calculateTpmMedians(self):
        tpmChrData = defaultdict(list)
        tpmChrMedians = {}
        medians = []
        for gene, chr in self.GeneModels.geneChromosomes.items():
            if gene in self.TPM.senseTPM:
                tpm = self.TPM.senseTPM[gene]
                tpmChrData[chr].append(tpm)
        for chr, values in tpmChrData.items():
            median = statistics.median(values)
            tpmChrMedians[chr] = median
            medians.append(median)
        tpmGenomeMedian = statistics.median(medians)
        return (tpmChrMedians, tpmGenomeMedian)



    def _calculatePloidy (self):
        chrPloidy = {}
        for chr, median in sorted(self.tpmChrMedians.items()):
            try:
                ploidy = (median/(self.tpmGenomeMedian/self.basePloidy))
            except ZeroDivisionError:
                print ("Median TPM value for chromosome {1} is {0:1.2f}. Cannot divide by {0:1.2f}. Setting ploidy value for chromosome {1} to 0.  It is recommended to check your htseq-count output and tpm values.\n".format(tpmGenomeMedian, chr))
                ploidy = 0
            chrPloidy[chr] = ploidy
        return chrPloidy


    
    def _calculateGeneCopyNumbers (self):
        geneCopyNumbers = defaultdict(list)
        for gene, data in sorted(self.TPM.senseTPM.items()):
            chr = self.GeneModels.geneChromosomes[gene]
            try:
                copyNumber = data/(self.tpmChrMedians[chr]/self.chrPloidy[chr])
            except ZeroDivisionError:
                print ("Median TPM value for {0} is {1:1.2f}.  Cannot divide by {1:1.2f}.  Setting gene copy number for {2} to 0.  It is recommended to check your htseq-count output and tpm values.\n".format(chr, self.tpmChrMedians[chr], gene))
                copyNumber = 0
            geneCopyNumbers[gene] = [data, chr, copyNumber]
        return geneCopyNumbers


    def _parseOrthologs(self):
        annotations = {}
        if self.orthologFile is None:
            print ("Ortholog mapping file not supplied. Gene copy numbers have already been calculated, but pooling over gene families wil not be carried out\n")
            return annotations
        try:
            with open (self.orthologFile) as f:
                next(f)
                for line in f:
                    data = line.rstrip().split('\t')
                    if data[2].startswith('OG6'):
                        annotations[data[0]] = data[2]
                f.close()
        except FileNotFoundError:
            print ("File {} containing gene to ortholog mappings could not be found. Gene copy numbers have already been calculated, but pooling over gene families will not be carried out\n".format(self.orthologFile))
        return annotations


    def _calculateCNVForOrthologGroups (self):
        clusterData = defaultdict(lambda: defaultdict(list))
        for gene, data in self.geneCopyNumbers.items():
            tpm, chr, copyNumber = data
            if gene in self.annotations:
                orthoMclId = self.annotations[gene]
            else:
                orthoMclId = gene #assume single copy
  
            clusterData[orthoMclId]['genes'].append(gene)
            clusterData[orthoMclId]['tpms'].append(data[0])
            clusterData[orthoMclId]['copyNumbers'].append(data[2])
            try:
                haploidNo = data[2]/self.chrPloidy[chr]
            except ZeroDivisionError:
                print("Cannot calculate haploid numbers for genes on chromosome {} because ploidy is estimated to be 0\nSetting haploid number for gene {} to NaN".format(chr, gene))
                haploidNo = float('NaN')
            clusterData[orthoMclId]['haploidNos'].append(haploidNo)

        for cluster, data in clusterData.items():
            geneCount = len(data['genes'])
            try:
                haploidNo = sum(data['haploidNos'])
            except TypeError:
                print ("Haploid numbers could not be calculated for some genes in cluster {}. Haploid number for this cluster has been set to NaN".format(cluster))
                haploidNo = float('NaN')
            geneDose = sum(data['copyNumbers'])
            
            for gene in data['genes']:
                self.geneCopyNumbers[gene] = self.geneCopyNumbers[gene] + [cluster, geneCount, haploidNo, geneDose, ', '.join(data['genes'])]             


    def ploidyPrinter(self, outputPrefix):
        outputFile = "{}_chrCNVs.tsv".format(outputPrefix)
        try:
            with open (outputFile, 'wt') as out:
                out.write("Chromosome\tRaw Copy Number\tRounded Copy Number\n")
                for chr, ploidy in sorted(self.chrPloidy.items()):
                    out.write("{0}\t{1:4.3f}\t{2}\n".format(chr, ploidy, int(ploidy+0.5)))
                out.close()
        except FileNotFoundError:
            raise SystemExit("File {} could not be opened for writing\n".format(outputFile))


    def cnvPrinter(self, outputPrefix):
        outputFile = "{}_geneCNVs.tsv".format(outputPrefix)
        try:
            with open (outputFile, 'wt') as out:
                out.write("Gene Id\tChromosome\tEstimated Ploidy\tTPM\tGene Dose\tHaploid Number\tOrtholog Id\tCopies in Reference\tCluster Gene Dose\tCluster Haploid Number\tGene List\n")
                for gene, data in self.geneCopyNumbers.items():
                    ploidy = self.chrPloidy[data[1]]
                    try:
                        haploidNo = data[2]/ploidy
                    except ZeroDivisionError:
                        haploidNo = float('NaN')
                    out.write("{0}\t{1}\t{2:4.3f}\t{3:4.3f}\t{4:4.3f}\t{5:4.3f}\t{6}\t{7}\t{8:4.3f}\t{9:4.3f}\t{10}\n".format(gene, data[1], ploidy, float(data[0]), data[2], haploidNo, data[3], data[4], data[6], data[5], data[7]))
                out.close()
        except FileNotFoundError:
            raise SystemExit("File {} could not be opened for writing\n".format(outputFile))
