# MSTS_count_TPM.py

This script counts TPM for each transcript from a Bam file and an annotation file in a valid gff3 format.

Counts are performed on each available transcript for a gene. No read distribution model are implemented, we count twice a read if 2 transcript isoforms are annotated (one per transcript on common regions). Counts are on exon per default, but you can choose cds with featType option. The mean fragment length is estimated for each transcript and then for the run. Only bases on feature (exon, cds) are taken into account for fragment length estimation. For example,in case where a fragment spans an intron, these bases will substracted of the whole length.
This script supports paired-end and single (experimental - not validated) data.

## Usage and options

### Usage:

`MSTS_count_TPM.py input.bam genes.gff3 > counts.tpm`

or

`MSTS_count_TPM.py input.bam genes.gff3 -t cds -m 10 -s reverse -v 2 > cds_counts.tpm`

### Arguments:

| Argument | Description |
| --------- | ----------- |
| `BamFile` | mapping file in bam format, sorted and indexed |
| `AnnotFile` | annotation file in gff3 format |


### Options:

| Option | Description |
| ------ | ----------- |
| `-t, --featType` | Feature type choice for counts [exon,cds], default=exon |
| `-m, --minNbFrags` | Minimum number of fragment per trancript to compute mean fragment length, [default=3] |
| `-s, --stranded` | Define which fragments take into account for computing mean fragment length and counts. If your reads are not oriented, set to "no" as default, if not use "yes" or "reverse". [default=no] |
| `-v, --verbosity` | increase output verbosity 1=error, 2=info, 3=debug |
| `--version` | tool version |
| `-h, --help` | help message |

## Inputs:

The Bam file must be preferentially obtained from a splice aware mapper (ex: STAR, HISAT, TopHat, ...)

Below a valid gff3 format for this script:

```
##gff-version 3
seq1    bioinfobioger   gene    19966   20759   .       -       .       ID=G1;Name=G1
seq1    bioinfobioger   mRNA    19966   20759   .       -       .       ID=T1;Name=T1;Parent=G1
seq1    bioinfobioger   exon    19966   20130   .       -       .       ID=G1E1;Parent=T1;gene_id=G1
seq1    bioinfobioger   exon    20186   20759   .       -       .       ID=G1E2;Parent=T1;gene_id=G1
seq1    bioinfobioger   CDS     20186   20668   .       -       .       ID=G1C1;Parent=T1;gene_id=G1
seq2    bioinfobioger   gene    93012   93200   .       +       .       ID=G2;Name=G2
seq2    bioinfobioger   mRNA    93012   93200   .       +       .       ID=T2;Name=T2;Parent=G2
seq2    bioinfobioger   exon    93012   93200   .       +       .       ID=G2E1;Parent=T2;gene_id=G2
seq2    bioinfobioger   CDS     93012   93200   .       +       .       ID=G2C1;Name=G2C1;Parent=T2;gene_id=G2
seq2    bioinfobioger   gene    137377  138896  .    +       .       ID=G3;Name=G3;
seq2    bioinfobioger   mRNA    137377  138896  .    +       .       ID=T3;Name=T3;Parent=G3;
seq2    bioinfobioger   exon    137377  137570  .       +       .       ID=G3E1;Parent=T3;gene_id=G3
seq2    bioinfobioger   exon    137683  137882  .       +       .       ID=G3E2;Parent=T3;gene_id=G3
seq2    bioinfobioger   exon    138004  138896  .    +       .       ID=G3E3;Parent=T3;gene_id=G3
seq2    bioinfobioger   CDS     137377  137570  .       +       .       ID=G3C1;Name=G3C1;Parent=T3;gene_id=G3
seq2    bioinfobioger   CDS     137683  137882  .       +       .       ID=G3C1;Name=G3C1;Parent=T3;gene_id=G3
seq2    bioinfobioger   CDS     138004  138896  .       +       .       ID=G3C1;Name=G3C1;Parent=T3;gene_id=G3
seq2    bioinfobioger   gene    139543  140931  .       +       .       ID=G4;Name=G4;
seq2    bioinfobioger   mRNA    139543  140931  .       +       .       ID=T4;Name=T4;Parent=G4;
seq2    bioinfobioger   mRNA    139543  143300  .      +       .       ID=T4.2;Name=T4.2;Parent=G4;
seq2    bioinfobioger   exon    139543  139555  .       +       .       ID=G4E1;Parent=T4,T4.2;gene_id=G4
seq2    bioinfobioger   exon    139642  139709  .       +       .       ID=G4E2;Parent=T4;gene_id=G4
seq2    bioinfobioger   exon    139642  139700  .       +       .       ID=G4E2.2;Parent=T4.2;gene_id=G4
seq2    bioinfobioger   exon    139790  140054  .       +       .       ID=G4E3;Parent=T4;gene_id=G4
seq2    bioinfobioger   exon    140150  140931  .       +       .       ID=G4E4;Parent=T4;gene_id=G4
seq2    bioinfobioger   exon    140900  143300  .       +       .       ID=G4E4;Parent=T4.2;gene_id=G4
seq2    bioinfobioger   CDS     139543  139555  .       +       .       ID=G4C1;Name=G4C1;Parent=T4;gene_id=G4
seq2    bioinfobioger   CDS     139642  139709  .       +       .       ID=G4C1;Name=G4C1;Parent=T4;gene_id=G4
seq2    bioinfobioger   CDS     139790  140054  .       +       .       ID=G4C1;Name=G4C1;Parent=T4;gene_id=G4
seq2    bioinfobioger   CDS     140150  140931  .       +       .       ID=G4C1;Name=G4C1;Parent=T4;gene_id=G4
seq2    bioinfobioger   CDS     139543  139555  .       +       .       ID=G4C2;Name=G4C2;Parent=T4.2;gene_id=G4
seq2    bioinfobioger   CDS     139642  139700  .       +       .       ID=G4C2;Name=G4C2;Parent=T4.2;gene_id=G4
seq2    bioinfobioger   CDS     140900  143250  .       +       .       ID=G4C2;Name=G4C2;Parent=T4.2;gene_id=G4
```

All features have an ID. All CDS, exon, mRNA have a Parent. Gene feature is also necessary. It is a very simple gff parser, you must respect the order of declaration: Gene, mRNA, exon, CDS. More complex parsing configuration will be implemented on demand, but it always better to work with well structured gff files. 



## Outputs

Tab delimited output:

```
ID      length  EffectiveLength counts  TPM
T4.2    2423    2222.5  4448    164407.204235
T4      1128    927.5   1054    93352.1936696
T2      189     201.5   4       1630.73268566
T3      1287    1086.5  9730    735666.440371
T1      483     282.5   17      4943.4290396
```
