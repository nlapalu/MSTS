# MSTS: MAINE-Seq Tool Suite

MSTS has been developed to analyse NGS data in the frame of MAINE-Seq experiments and to propose some utilities to draw graph and perform simple statistics on your data.

MSTS is developed by the BioinfoBIOGER plateform (N.Lapalu, A.Simon) at [INRA-BIOGER](http://www.versailles-grignon.inra.fr/bioger). Please do not hesitate to contact us (nlapalu at inra dot fr) if you have any comments or questions.

## Dependencies and external tools

### External tools

* samtools
* wigToBigWig (Kent's tools)

## Install  (DEV mode !!!)

### Prerequesites:

TODO

### Download

```
wget https://github.com/nlapalu/MSTS/archive/develop.zip
unzip develop.zip
cd MSTS-develop/
```

### Test the code

### Install as root

### Install as user

```
python setup.py install --prefix=/home/nlapalu/test/MSTS
export PYTHONPATH=/home/nlapalu/test/MSTS/lib/python2.7/site-packages/
export PATH=$PATH:/home/nlapalu/test/MSTS/bin
```

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

`MSTS_converter.py mapping.sorted.bam -m fragment-middle -w 20 -p mapping -g assembly.genome --wig --size`

#### convert wig file to bigWig file

`wigToBigWig mapping.wig assembly.genome mapping.bw`

#### draw a phasogram 

`MSTS_phasogram.py mapping.bw -w 1200 -o mapping.phasogram.png -t "phasogram - mapping" -v 2 --flush --regression > mapping.phaso`

#### analyze relation between Transcript Expression level and nucleosome occupancy

If you have RNA-Seq data and MAINE-Seq data, you could draw phasograms with intervals of expression levels

Generate TPM count file from RNA-Seq mapped data (sorted, indexed bam) and annotation file

`MSTS_count_TPM.py RNA_mapping.sorted.bam genes.gff3 -v 2 > counts.tpm`

Export transcript list of IDs for defined Expression Level intervals (x>50,50>=x>5,5>=x>1,x<=1)

`tail -n+2 counts.tpm | awk -F"\t" '{if($5 > 50){print $1}}' > 50.tpm` 

`tail -n+2 counts.tpm | awk -F"\t" '{if($5 <= 50 && $5 > 5){print $1}}' > 50-5.tpm` 

`tail -n+2 counts.tpm | awk -F"\t" '{if($5 <= 5 && $5 > 1){print $1}}' > 5-1.tpm` 

`tail -n+2 counts.tpm | awk -F"\t" '{if($5 < 1){print $1}}' > 1.tpm` 

Generate phasogram for each list

`MSTS_feature_phasogram.py mapping.bw genes.gff3 -v 2 -o myfeaturestartphasogram1.png -t "phasogram on transcript, start as pivot, TPM > 50"  -ft mRNA -l 50.tpm --context --GaussianSmoothing`

`MSTS_feature_phasogram.py mapping.bw genes.gff3 -v 2 -o myfeaturestartphasogram1.png -t "phasogram on transcript, start as pivot, 50>TPM>5"  -ft mRNA -l 50-5.tpm --context --GaussianSmoothing`

`MSTS_feature_phasogram.py mapping.bw genes.gff3 -v 2 -o myfeaturestartphasogram1.png -t "phasogram on transcript, start as pivot, 5>TPM>1"  -ft mRNA -l 5-1.tpm --context --GaussianSmoothing`

`MSTS_feature_phasogram.py mapping.bw genes.gff3 -v 2 -o myfeaturestartphasogram1.png -t "phasogram on transcript, start as pivot, TPM < 1"  -ft mRNA -l 1.tpm --context --GaussianSmoothing`

## MSTS tools

* [MSTS_converter.py](doc/MSTS_converter.md)
* [MSTS_phasogram.py](doc/MSTS_phasogram.md)
* [MSTS_feature_phasogram.py](doc/MSTS_feature_phasogram.md)
* [MSTS_count_TPM.py](doc/MSTS_count_TPM.md)
* [MSTS_dinuc_frequency.py](doc/MSTS_dinuc_frequency.md)
* [MSTS_detect_nucleosomes.py](doc/MSTS_dinuc_frequency.md)

## MSTS overview

<img src="doc/images/MSTS_overview.png" width="850">

## References

