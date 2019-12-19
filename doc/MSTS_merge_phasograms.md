Here a simple script to generate a merged phasogram from several runs of MSTS_feature_phasogram.py 

write in list.txt all the *.phaso file to merge with the associated legend (tab delimited).

__*list.txt*__:

```
1.phaso	1TPM>X
5-1.phaso	5TPM>X>1TPM
50-5.phaso	50TPM>X>5TPM
50.phaso	X>50TPM
```

(up to 6 conditions) and run the command:

`MSTS_merge_phasograms.py list.txt`

Your phasogram will be exported in phaso.merge.png

