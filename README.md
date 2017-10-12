# MSTS: MAINE-Seq Tool Suite

MSTS has been developed to analyse NGS data in the frame of MAINE-Seq experiments and to propose some utilities to draw graph and perform simple statistics on your data.

MSTS is developed by the BioinfoBIOGER plateform (N.Lapalu, A.Simon) at [INRA-BIOGER](http://www.versailles-grignon.inra.fr/bioger). Please do not hesitate to contact us (nlapalu at inra dot fr) if you have any comments or questions.

## Dependencies and external tools



## Install

## Protocole to analyze MAINE-Seq data

### Introduction

some comments on mapping step ...

some comments on quality criteria and metrics / graph to control quality

### Single reads

not fully supported at this time

### Paired reads

#### convert your mapping file

Be sure your mapping file is sorted and indexing, if not:

`samtools sort mapping.bam > mapping.sorted.bam`

`samtools index mapping.sorted.bam `

then:

`MSTS_converter.py mapping.sorted.bam -m fragment-middle -w 20 -p mapping --wig --size`

#### draw a phasogram 


## MSTS tools

* [MSTS_converter.py](doc/MSTS_converter.md)
* [MSTS_phasogram.py](doc/MSTS_phasogram.md)

## References

