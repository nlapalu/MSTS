# MSTS_converter.py

blabla

## Usage and options

### Usage:

`MSTS_converter.py input.sorted.bam`

### Options

| Option | Description |
| ------ | ----------- |
| `-p, --prefix` | prefix for output files, default=[out.] |
| `-m, --mode` | counting mode type: [single, single-expanded, fragment, fragment-middle], default=single |
| `--bed` | output bed file |
| `--wig` | output wig file |
| `--cov` | output coverage file |
| `--size` | output size file |
| `-g, --genome` | genome file used to sort references in defined rank, expected format: \<refname\>\<TAB\>\<size\> |
| `-w, --window` | for fragment middle mode only, +- window size on each side |
| `-k, --keepPosBigBedFile` | bigBed file used to filter positions, keep only positions in file, warning only usable with seq < 100Mb  |
| `--minFSize` | minimum fragment size to keep |
| `--maxFSize` | maximum fragment size to keep |
| `--minDepCov` | minimum depth of coverage to keep, warning it is a posteriori filter, compute on reduce fragment if mode=fragment-middle |
| `--maxDepCov` | maximum depth of coverage to keep, warning it is a posteriori filter, compute on reduce fragment if mode=fragment-middle |
| `-v, --verbosity` | increase output verbosity 1=error, 2=info, 3=debug |
| `--version` | tool suite version |
| `-h, --help` | help message |
