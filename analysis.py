#change the data matrix file name, #dimensions for pca, #clusters before using

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from multiprocessing import Pool
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, fcluster

#load data matrix
rate_bin_ch=np.load(binlevel_mCH_ratio_matrix_file)
rate_gene_ch=np.load(genelevel_mCH_ratio_matrix_file)
read_bin_ch=np.load(binlevel_read_count_matrix_file)
read_gene_ch=np.load(genelevel_read_count_matrix_file)
meta=np.load(meta_data)

#keep bins with >100bp coverage in >97.5% cells
filter=np.sum(read_bin_ch>100,axis=0)>0.975*len(read_bin_ch)
print(sum(filter))
rateb_ch=rate_bin_ch[:,filter]
rateb_ch=np.divide(rateb_ch.T,meta[:,8].astype(float)).T
bin=bin_all[filter]

#keep bins with >100bp coverage in >90% cells
filter=np.sum(read_gene_ch>100,axis=0)>0.9*len(read_gene_ch)
print(sum(filter))
rateg_ch=rate_gene_ch[:,filter]
rateg_ch=np.divide(rateg_ch.T,meta[:,8].astype(float)).T
gene=gene_all[filter]

#pca of data normalized matrix
pca=PCA(n_components=len(rateb_ch)-1)
rateb_reduce_ch=pca.fit_transform(rateb_ch)

#exclude outliers as the single cells who merged into the hierarchy at the last nc branch points
dim=number_of_dimensions
nc=number_of_clusters
cluster=fcluster(Z,t=nc,criterion='maxclust')-1
count=np.array([sum(cluster==i) for i in range(nc)])
outlier=np.array([(count[cluster[i]]>5) for i in range(len(cluster))])

#tsne and hierarchical clustering
tsne=TSNE(n_components=2,perplexity=50,random_state=0)
y=tsne.fit_transform(rateb_reduce_ch[outlier,:50])
Z=linkage(rateb_reduce_ch[outlier,:dim],'ward')
nc=number_of_clusters
cluster=fcluster(Z,t=nc,criterion='maxclust')-1
count=np.array([sum(cluster==i) for i in range(nc)])
tot=0
plt.figure()
for i in range(nc):
	if count[i]>5:
		cell=(cluster==i)
		plt.scatter(y[cell,0],y[cell,1],s=5,c=color[tot],alpha=0.5,edgecolors='none',label='cluster'+str(tot+1))
		tot+=1

plt.legend()
plt.legend(bbox_to_anchor=(-0.3,1), loc="upper left")
plt.tight_layout()
plt.savefig(file_dir+'.pca'+str(dim)+'.hierar'+str(nc)+'.pdf',bbox_inches="tight")
plt.close()
print(nc,tot,count)

#plot marker gene mCH ratio to identify specific layers
marker=['Satb2','Cux1','Cux2','Rorb']
for k in marker:
	mch=rateg_ch[:,gene[:,-1]==k][:,0]
	plt.figure()
	plt.scatter(y_ch[:,0],y_ch[:,1],s=5,c=mch,alpha=0.5,edgecolors='none',cmap=cm.jet)
	plt.clim([np.percentile(mch,1),np.percentile(mch,99)])
	plt.colorbar()
	plt.title(k)
	plt.savefig(file_dir+'.mrk.'+k+'.pdf')
	plt.close()

