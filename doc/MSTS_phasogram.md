# MSTS_phasogram.py

This script allows you to draw a phasogram from a bigWig file.

## Usage and options

### Usage:

`MSTS_phasogram.py input.bw -w 1000 -o myphasogram.png -t "phasogram wild-type data" -v 2`

or 

`MSTS_phasogram.py input.bw  -w 1000 -o myphasogram.png -t "phasogram wild-type data" -v 2  --flush --regression`


We recommend not to exceed 1kb as analyzed window. Beyond this value, the signal becomes less reliable.

### Options:

| Option | Description |
| ------ | ----------- |
| `-w, --window` | window size to compute phases |
| `--flush` | print phases on stdout to save in file, > phases.out |   
| `--regression` | detect peaks and perform a regression. Regression curve drawn on the graph |
| `--norm` | normalize the signal with the mean coverage |
| `-o, --out` | name of output graph |
| `-t, --title` | title text |
| `-x, --xax` | x axis text |
| `-y, --yax` | y axis text |
| `-b, --bigBed` | bigBedFile, use to limit phasogram to specific regions |
| `-v, --verbosity` | increase output verbosity 1=error, 2=info, 3=debug |
| `--version` | tool suite version |
| `-h, --help` | help message |

## Outputs

#### simple phasogram
![image](images/myphasogram.png)

#### phasogram (--regression)
![image](images/myphasogram2.png)
