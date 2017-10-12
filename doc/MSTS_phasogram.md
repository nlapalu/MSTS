# MSTS_phasogram.py

This script allows you to draw a phasogram from a bigWig file.

## Usage and options

### Usage:

`MSTS_phasogram.py input.bw -w 1000 -o myphasogram.png -t "phasogram wild-type data" -v 2`


### Options:

| Option | Description |
| ------ | ----------- |
| `-w, --window` | window size to compute phases |
| `--flush` |print phases on stdout to save in file, > phases.out |
| `-o, --out` | name of output graph |
| `-t, --title` | title text |
| `-x, --xax` | x axis text |
| `-y, --yax` | y axis text |
| `-v, --verbosity` | increase output verbosity 1=error, 2=info, 3=debug |
| `--version` | tool suite version |
| `-h, --help` | help message |

## Outputs

![image](images/myphasogram.png)
