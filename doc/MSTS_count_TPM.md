# MSTS_count_TPM.py

This script counts TPM for each transcript from a Bam file and an annotation file in a valid gff3 format.

Counts are performed on each available transcript of a gene. No read distribution model are implemented, we count twice a read if 2 transcript isoforms are annotated (one per transcript on common regions). Counts are on exon per default, but you can choose cds with featType option. The mean fragment length is estimated for each transcript and then for the run. Only bases on feature (exon, cds) are taken into account for fragment length estimation. For example,in case where a fragment spans an intron, these bases will substracted of the whole length.
These script support paired-data and single(experimental - not validated)

## Usage and options

### Usage:

`MSTS_count_TPM.py input.bam genes.gff3 > counts.tpm`



### Options:

| Option | Description |
| ------ | ----------- |


## Outputs


