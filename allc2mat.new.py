#python allc2mat.new.py ${allc_file_name} ${ncpus} ${bed_file_name}
#bed file need to be sorted, allc file need to be put in the same directory with allc index file
#this will generate sample.mC.txt with mCH and mCG basecalls and ratio, and sample.tot.txt with global mCCC, mCH, mCG level

import os
import sys
import gzip
import numpy as np
from multiprocessing import Pool

def add(x,mode):
	if mode==1:
		if x[5]>2 or x[3][1]=='N':
			y=[0,0,0,0]
		elif x[3][1]!='G':
			y=[x[4],x[5],0,0]
		else:
			y=[0,0,x[4],x[5]]
		return y
	if mode==0:
		y=[0 for i in range(6)]
		if x[5]>2:
			return y
		if x[3]=='CCC':
			y[0]+=x[4]
			y[1]+=x[5]
		if x[3][1]=='G':
			y[4]+=x[4]
			y[5]+=x[5]
		elif x[3][1]!='N':
			y[2]+=x[4]
			y[3]+=x[5]
		return y

def mapping(chr):
	global bin,out_folder,idx,sample,bed
	fin=gzip.open(sys.argv[1])
	fout=open(out_folder+bed+'_'+chr+'.mC.txt','w')
	fin.seek(idx[chr])
	x=bin[bin[:,0]==('chr'+chr),1:].astype(int)
	print(chr)
	ans0=[0 for i in range(6)]
	tmp=[chr]+[-1 for i in range(6)]
	y=[tmp[:]]
	for i in range(len(x)):
		ans=[0,0,0,0]
		if len(y)>0 and y[-1][1]>x[i,1]:
			j=len(y)
			while j-1>=0 and y[j-1][1]>x[i,1]:
				j-=1
			while j-1>=0 and y[j-1][1]>x[i,0]:
				j-=1
				result=add(y[j],1)
				ans=[ans[k]+result[k] for k in range(4)]
			y=y[j:]
		elif len(y)>0 and y[-1][1]>x[i,0]:
			j=len(y)
			while j-1>=0 and y[j-1][1]>x[i,0]:
				j-=1
				result=add(y[j],1)
				ans=[ans[k]+result[k] for k in range(4)]
			y=y[j:]
		elif len(tmp)==7 and tmp[0]==chr:
			y=[]
			while tmp[1]<=x[i,0]:
				tmp=fin.readline().strip().split('\t')
				if len(tmp)<7 or tmp[0]!=chr:
					break
				tmp[1],tmp[4],tmp[5]=int(tmp[1]),int(tmp[4]),int(tmp[5])
				result=add(tmp,0)
				ans0=[ans0[k]+result[k] for k in range(6)]
			if len(tmp)==7 and tmp[0]==chr:
				if tmp[1]<=x[i,1]:
					result=add(tmp,1)
					ans=[ans[k]+result[k] for k in range(4)]
				y.append(tmp)
		if len(tmp)==7 and tmp[0]==chr:
			while tmp[1]<=x[i,1]:
				tmp=fin.readline().strip().split('\t')
				if len(tmp)<7 or tmp[0]!=chr:
					break
				tmp[1],tmp[4],tmp[5]=int(tmp[1]),int(tmp[4]),int(tmp[5])
				if tmp[1]<=x[i,1]:
					result=add(tmp,1)
					ans=[ans[k]+result[k] for k in range(4)]
				result=add(tmp,0)
				ans0=[ans0[k]+result[k] for k in range(6)]
				y.append(tmp)
		if ans[1]>0 and ans[3]>0:
			result=[float(ans[0])/ans[1],float(ans[2])/ans[3]]
		elif ans[1]>0 and ans[3]==0:
			result=[float(ans[0])/ans[1],0.0]
		elif ans[1]==0 and ans[3]>0:
			result=[0.0,float(ans[2])/ans[3]]
		else:
			result=[0.0,0.0]
		fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(ans[0],ans[1],result[0],ans[2],ans[3],result[1]))
	while len(tmp)==7 and tmp[0]==chr:
		tmp=fin.readline().strip().split('\t')
		if len(tmp)<7 or tmp[0]!=chr:
			break
		tmp[1],tmp[4],tmp[5]=int(tmp[1]),int(tmp[4]),int(tmp[5])
		result=add(tmp,0)
		ans0=[ans0[k]+result[k] for k in range(6)]
	fin.close()
	fout.close()
	return(ans0)

sample=sys.argv[1].split('/')[-1][5:-7]
bed=sys.argv[3].split('/')[-1][:-4]
ncpus=int(sys.argv[2])
idx={}
out_folder='sc_matrix/'+sample+'_'
fin=open(sys.argv[1]+'.idx')
for line in fin:
	tmp=line.strip().split()
	if len(tmp)!=2:
		break
	idx[tmp[0]]=int(tmp[1])
fin.close()
chr=[str(i+1) for i in range(19)]+['X']
bin=np.loadtxt(sys.argv[3],dtype=np.str)
p=Pool(ncpus)
result=p.map(mapping,chr)
p.close()
result=np.sum(np.array(result),axis=0)
fout=open(out_folder+bed+'.tot.txt','w')
fout.write(sample+'\t')
fout.write('\t'.join([str(result[i]) for i in range(6)])+'\t')
fout.write('{0}\t{1}\t{2}\n'.format(float(result[0])/result[1],float(result[2])/result[3],float(result[4])/result[5]))
fout.close()
os.system('cat $(ls '+out_folder+bed+'_*.mC.txt | sort) > '+out_folder+bed+'.mC.txt')
os.system('rm '+out_folder+bed+'_*.mC.txt')
