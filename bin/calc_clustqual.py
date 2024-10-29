#TimesVector
#
#Created by Inuk Jung on 2016-08-3.
#Copyright (c) 2016 Inuk Jung. All rights reserved.
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import math
from scipy.spatial import distance
import numpy as np

class Cluster(object):
	def __init__(self, id, p, tp):
		self.kclid=id
		self.pheno={}
		self.base_cent=[]
		self.n_pheno=p
		self.n_tp=tp

class Phenotype(object):
	def __init__(self):
		self.vect=[]
		self.cent=[]

def magnitude(v):
	return math.sqrt(sum(v[i]*v[i] for i in range(len(v))))

def normalize(v):
	vmag=magnitude(v)
	return [v[i]/vmag for i in range(len(v))]

def vect_mu(v_list, l):
	v=np.array(v_list)
	vsum=np.ndarray.sum(v,axis=0)
	return normalize(vsum)

def getgenexpr(fname):
	f=open(fname, "r")
	dict={}

	for lidx, line in enumerate(f):
		if lidx==0: continue
		tok=line.strip().split()
		gid=tok[0]
		# expr=map(float,tok[1:])
		# expr=[x+1.0 for x in expr]
		expr=[float(x)+1.0 for x in tok[1:]]
		profile="dummy"

		dict[gid]=[expr, profile]
	f.close()

	return dict
	
def read_clusters(fname, g, p, tp): # cluster file, gene expression file, phenotypes, time points
	f=open(fname, "r")
	kc={}
	base_cent=[]	

	for line in f:
		if "x" in line: continue
		tok=line.strip().split()
		kclid=tok[1]
		gid=tok[0]

		if kclid not in kc:
			kc[kclid]=Cluster(kclid, p, tp)

		for i in range(0,p):
			if i not in kc[kclid].pheno:
				kc[kclid].pheno[i]=Phenotype()

			# build centroid per phenotype 
			expr=map(float, g[gid][0])
			v=normalize(expr[i*tp:(i*tp)+tp])
			kc[kclid].pheno[i].vect.append(v)
	f.close()

	for k in kc:
		for i in range(0,p):
			kc[k].pheno[i].cent=vect_mu(kc[k].pheno[i].vect, p*tp)

		cos_dist=1
		for i in range(0,p):
			for j in range(0,p):
				if i==j: continue
				cur_dist=distance.cosine(kc[k].pheno[i].cent, kc[k].pheno[j].cent)
				if cur_dist < cos_dist:
					cos_dist=cur_dist	# update minimum cosine distance
					base_cent=kc[k].pheno[i].cent

		kc[k].base_cent=base_cent


	return kc

def print_dist(c, e_f, c_f):

	p=len(c.pheno)
	tp=(len(c.pheno[0].vect[0])*p)/p
	cent=c.base_cent
	e_f.write("gene\t")
	count=1
	for i in c.pheno:
		for v in c.pheno[i].vect:
			e_f.write("%d\t"%(count))
			count+=1
	e_f.write("\n")

	e_f.write("%s\t"%(c.kclid))
	for i in c.pheno:
		for v in c.pheno[i].vect:
			dist=distance.cosine(v, cent)
			e_f.write("%f\t"%(dist))
	e_f.write("\n")

	# condition labels
	for i in c.pheno:
		for v in c.pheno[i].vect:
			c_f.write("%d\t"%(i))		
	c_f.write("\n")


### MAIN ###
expr_f=sys.argv[1]
fname=sys.argv[2]
outdir=sys.argv[3]
pheno=int(sys.argv[4])
tp=int(sys.argv[5])

gene_expr_dict=getgenexpr(expr_f)
clusters=read_clusters(fname, gene_expr_dict, pheno, tp) 
bdir=outdir+"/metadata/"

for c in clusters:
	out_expr_f=open(bdir+c+"_expr.dat", "w")
	out_cond_f=open(bdir+c+"_cond.dat", "w")
	print_dist(clusters[c], out_expr_f, out_cond_f)

out_expr_f.close()
out_cond_f.close()



