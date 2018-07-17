# snmC-seq2

snmC-seq2 is an improved technique for single nucleus methylome sequencing. The code used for the manuscript is in this repository.

## Computing bin/gene level methylation ratio

The scirpt allc2mat.new.py is for calculating bin/gene level methylation ratio. The output of allc2mat.new.py includes the following two files.

### sample.mC.txt
each rows is a region from the input bed file

each columns are
mCH basecalls
CH basecalls
mCH/CH ratio

mCG basecalls

CG basecalls

mCG/CG ratio

### sample.tot.txt
each columns are

sample name

global mCCC basecalls

global CCC basecalls

global mCH basecalls

global CH basecalls

global mCG basecalls

global CG basecalls

global mCCC/CCC ratio

global mCH/CH ratio

global mCG/CG ratio
