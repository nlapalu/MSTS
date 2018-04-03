




# MSTS: MAINE-Seq Tool Suite

MSTS has been developed to analyse NGS data in the frame of MAINE-Seq experiments and to propose some utilities to draw graph and perform simple statistics on your data.

MSTS is developed by the BioinfoBIOGER plateform (N.Lapalu, A.Simon) at [INRA-BIOGER](http://www.versailles-grignon.inra.fr/bioger). Please do not hesitate to contact us (nlapalu at inra dot fr) if you have any comments or questions.

# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Tool documentations](#tool-documentations)
* [Protocole to analyze MAINE-Seq data](#protocole-to-analyze-maine-seq-data)
	* [Introduction](#introduction)
	* [Mapping sequencing reads](#mapping-sequencing-reads)
	* [Draw phasogram and get nucleosome spacing](#draw-phasogram-and-get-nucleosome-spacing)
	* [Draw phasogram and analyze feature specificity](#draw-phasogram-and-analyze-feature-specificity)
	* [Analyze relation between Transcript Expression level and nucleosome occupancy](#analyze-relation-between-transcript-expression-level-and-nucleosome-occupancy)
	* [Detect and classify nucleosomes](#detect-and-classify-nucleosomes)
	* [Analyze di-nucleotide composition](#analyze-di-nucleotide-composition)
		* [From mapped reads](#from-mapped-reads) 
		* [From detected nucleosome positions](from-detected-nucleosome-positions)
* [References](#references)




## Dependencies and external tools

### External tools

* samtools
* wigToBigWig (Kent's tools)
* bedToBigBed (Kent's tools)

## Install  (DEV mode !!!)

### Prerequesites:

* numpy
* scipy
* scikit-learn
* fastcluster (http://danifold.net/fastcluster.html)
* matplotlib
* pysam
* pyBigWig

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

# Protocole to analyze MAINE-Seq data

## Introduction

some comments on mapping step ...

some comments on quality criteria and metrics / graph to control quality

## Mapping sequencing reads

Be sure your mapping file is sorted and indexing, if not:

`samtools sort mapping.bam > mapping.sorted.bam`

`samtools index mapping.sorted.bam `

then:

### For single reads

`MSTS_converter.py mapping.sorted.bam -m single-expanded -p mapping -g assembly.genome --wig --size`

#### For paired reads

`MSTS_converter.py mapping.sorted.bam -m fragment-middle -w 20 -p mapping -g assembly.genome --wig --size --bed`

### Convert wig file to bigWig file

`wigToBigWig mapping.wig assembly.genome mapping.bw`

## Draw phasogram and get nucleosome spacing 

`MSTS_phasogram.py mapping.bw -w 1200 -o mapping.phasogram.png -t "phasogram - mapping" -v 2 --flush --regression > mapping.phaso`

## Draw phasogram and analyze feature specificity

...


## Analyze relation between Transcript Expression level and nucleosome occupancy

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

## Detect and classify nucleosomes


`MSTS_detect_nucleosomes.py`  


## Analyze di-nucleotide composition

The di-nucleotide pattern AT/GC with a frequence of 10 bp for nucleosome fixation site has been largely described (ref) and proposed to be used as data quality control (*S.Hu et al. 2017*). You can perform a such analysis on your mapping data converted in bigBed file with MSTS_converter.py

### From mapped reads 

```
# convert your mapping bed file to bigBed and run MSTS_dinuc_frequency
bedToBigBed mapping.bed assembly.genome mapping.bb
MSTS_dinuc_frequency.py mapping.bb 
```

<img src="doc/images/detect_dinuc_mapping_40_ATGC_Normalized.png" width="425"> <img src="doc/images/detect_dinuc_mapping_40_ATGC_Correlogram.png" width="425">

### From detected nucleosome positions

You can also control your nucleosome detection and classification with this tool. We expect a better signal for well positioned nucleosome compare to bad/loosely positioned. To do that run MSTS_detect_nucleosomes.py with --bed option and convert them to bigBed. We present below the dinucleotide frequency analyzed for 4 types of clusters (bad, fuzzy, well, very-well):

#### bad positioned nucleosomes:

`python MSTS_dinuc_frequency.py genome.fasta detect.k1-bad.cluster.sorted.bb --pAutocorMix --pFreqNormMix -p detect_dinuc_bad -ami -65 -amx 65`

<img src="doc/images/detect_dinuc_bad_ATGC_Normalized.png" width="425"> <img src="doc/images/detect_dinuc_bad_ATGC_Correlogram.png" width="425">

#### fuzzy positioned nucleosomes:*

`python MSTS_dinuc_frequency.py genome.fasta detect.k2-fuzzy.cluster.sorted.bb --pAutocorMix --pFreqNormMix -p detect_dinuc_fuzzy -ami -65 -amx 65`

<img src="doc/images/detect_dinuc_fuzzy_ATGC_Normalized.png" width="425"> <img src="doc/images/detect_dinuc_fuzzy_ATGC_Correlogram.png" width="425">

#### well positioned nucleosomes:*

`python MSTS_dinuc_frequency.py genome.fasta detect.k3-well.cluster.sorted.bb --pAutocorMix --pFreqNormMix -p detect_dinuc_well -ami -65 -amx 65`

<img src="doc/images/detect_dinuc_well_ATGC_Normalized.png" width="425"> <img src="doc/images/detect_dinuc_well_ATGC_Correlogram.png" width="425">

#### very-well positioned nucleosomes:*

`python MSTS_dinuc_frequency.py genome.fasta detect.k4-very-well.cluster.sorted.bb --pAutocorMix --pFreqNormMix -p detect_dinuc_very-well -ami -65 -amx 65`

<img src="doc/images/detect_dinuc_very-well_ATGC_Normalized.png" width="425"> <img src="doc/images/detect_dinuc_very-well_ATGC_Correlogram.png" width="425">




# Tool documentations

* [MSTS_converter.py](doc/MSTS_converter.md)
* [MSTS_phasogram.py](doc/MSTS_phasogram.md)
* [MSTS_feature_phasogram.py](doc/MSTS_feature_phasogram.md)
* [MSTS_count_TPM.py](doc/MSTS_count_TPM.md)
* [MSTS_dinuc_frequency.py](doc/MSTS_dinuc_frequency.md)
* [MSTS_detect_nucleosomes.py](doc/MSTS_detect_nucleosomes.md)

## MSTS overview

<img src="doc/images/MSTS_overview.png" width="850">

# References

* S. Hu, X. Chen, J. Liao, Y. Chen, C. Zhao, and Y. Zhang, “CAM: A quality control pipeline for MNase-seq data,” PLoS One, vol. 12, no. 8, p. e0182771, Aug. 2017.

