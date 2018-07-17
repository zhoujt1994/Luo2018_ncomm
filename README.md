# snmC-seq2

snmC-seq2 is an improved technique for single nucleus methylome sequencing. The code used for the manuscript is in this repository.

## Computing bin/gene level methylation ratio

The scirpt allc2mat.new.py is for calculating bin/gene level methylation ratio. It can be used as
```
python allc2mat.new.py allc_sample.tsv.gz #cpus regions.bed
```
It requires the bed file to be sorted, and allc file need to be put in the same directory with allc index file. The output of allc2mat.new.py includes the following two files.

### sample_region.mC.txt
each rows is a region from the input bed file

each columns are

1. mCH basecalls

2. CH basecalls

3. mCH/CH ratio

4. mCG basecalls

5. CG basecalls

6. mCG/CG ratio

### sample_region.tot.txt
each columns are

1. sample name

2. global mCCC basecalls

3. global CCC basecalls

4. global mCH basecalls

5. global CH basecalls

6. global mCG basecalls

7. global CG basecalls

8. global mCCC/CCC ratio

9. global mCH/CH ratio

10. global mCG/CG ratio

## Clustering and visualization

Code for clustering and visualization of gene level mCH ratio is in analysis.py, which takes the cell x bin/gene matrix of methylation level as input.

