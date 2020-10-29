# calculateCNVs.py
Estimate chromosome ploidy and gene copy numbers from coverage data. Provision of ortholog group assignments for each gene enables pooling of data over gene families.

In organisms with very plastic genomes where copy number in tandem arrays is unstable and tandem arrays may also be collapsed in the reference, methods for finding segmental duplications based on clear regions of higher than expected coverage are often not informative due to noise. This method focusses instead on estimating the copy number for each gene based on read counts per gene.


## Input
The input for this script is read counts in htseq-count format. To generate this file, it is recommended to align your reads in such a way that multi-mapping reads are aligned only once and ensure that these are included in the counts. Inclusion of an ortholog mapping file ensures that the contribution of these reads is distributed over all members of the gene family from which they originated.

It is also recommended to filter your data to exclude genes on small unassigned contigs, particularly if these are repetitive, and to include only protein-coding genes.

In the example given, ortholog assignments for genes are obtained from orthoMCL v6. Other methods for grouping genes into families can be used as long as the file format is consistent with that given in the example data, and the identifier for each group is unique.

See exampleData for examples of input file formats.

## Usage

```
usage: calculateCNVs.py [-h] --genome GENOME --input INPUT
                        [--orthologFile ORTHOLOGFILE] --outputPrefix
                        OUTPUTPREFIX [--ploidy PLOIDY]

Estimate chromosome and gene copy number from read counts

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       Reference genome annotation in GFF3 format. (default:
                        None)
  --input INPUT         Read counts per gene or transcript in htseq-count
                        format (default: None)
  --orthologFile ORTHOLOGFILE
                        File containing gene to ortholog mappings (default:
                        None)
  --outputPrefix OUTPUTPREFIX
                        Prefix for all output files (default: None)
  --ploidy PLOIDY       Canonical ploidy for this organism (default: 2)
  ```
  
  Example command line using example data:
  ``./calculateCNVs.py --genome exampleData/TriTrypDB-48_TbruceiTREU927.gff --input exampleData/counts_example.txt --orthologFile exampleData/GeneOrthologMap.txt --outputPrefix example``
  
By default, the script assumes you are working with a diploid organism. For organisms where this isn't true, the ploidy can be changed using the ```--ploidy``` flag. 

## Output
Three tab-delimited files are produced, using the prefix given in the command line.

A file showing the calculated TPM values for each gene.

A file showing the estimated ploidy for each chromosome based on coverage.  Note that large segmental duplications may result in estimated ploidy values that are between expected values (e.g., 2.5 in the case of a segmental duplication on a diploid chromosome)

A file showing the estimated copy number value for each gene. The columns are:
* Gene Id: the identifier for this gene
* Chromosome: the chromosome or contig on which this gene resides
* Estimated ploidy: the estimated ploidy for this chromosome
* TPM: the TPM value calculated for this gene
* Gene Dose: the total number of estimated copies of this gene taking into account the estimated ploidy
* Haploid Number: the estimated number of copies of this gene per chromsome copy
* Ortholog Id: the ortholog identifier for this gene
* Copies in reference: the number of genes sharing this ortholog identifier appearing in the reference annotation
* Cluster Gene Dose: the gene dose estimated across all genes with this ortholog identifier (note for single copy genes, this will be the same as the per gene estimate)
* Cluster Haploid Dose: the haploid number estimated across all genes with this ortholog identifier (note for single copy genes, this will be the same as the per gene estimate)
* Gene List: List of genes in this ortholog group

## Known Issues
* Copy numbers are often underestimated on smaller chromosomes. I need to find a way to normalise for chromosome length.
* Care should be taken interpreting genes or gene families where the estimated number of copies is predicted to be lower than the number of copies in the reference, especially where these genes are in multigene families. Poor mappability in these regions may result in low coverage.
* We have been using non-parametric tests (Mann Witney U) accompanied by empirical p-value estimation by permutation to compare estimated copy numbers between two groups of replicates. A better approach might be to estimate copy numbers once for each condition based on mean or median TPM values across replicates.
