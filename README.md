# CORGi

CORGi requires Matplotlib, Bokeh, and local installations of BLAST, SAMTools, and optionally IGV.

python corgi.py \
* -r reference.fa
* -b input.bam
* -c chromosome
* -p start & end coordinates
* -o outputDir/
* --skip-igv skip IGV screenshot
* --skip-bam skip read extraction

example:
python corgi.py -r hg38.fa -b my.bam -c chr1 -p 1000 2000 -o myOut/
