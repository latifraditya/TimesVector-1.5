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
import random
import numpy as np
import scipy.stats
from scipy.spatial import distance
#import refine


class Cluster(object):
	def __init__(self, p, tp):
		self.type="NEP"		# initialize to NEP
		self.sub_type="NEP"	# initialize to NEP
		self.pheno={}
		self.k_cent=[]
		self.base_cent=[]
		self.n_pheno=p
		self.n_tp=tp
		self.p_val=-1.0
		self.cdist=1.0
		self.gene_count=0

		self.k_skmeans_cent=[]	# the centroid of the phenotype appended expression values (long centroid)
		self.k_gid={}			# the whole phenotype & time expression values of a gene
		self.m=1.0
		self.ml=1.0
		self.mh=1.0
		self.max=1.0

		self.p_m=1.0
		self.p_ml=1.0

class Phenotype(object):
	def __init__(self):
		self.vect=[]
		self.cent=[]
		self.gid=[]
		self.mag=[]
		self.corrdat=[]
		self.mu_dist=1.0
		self.sd_dist=1.0
		self.max_dist=1.0
		self.conf_dist=1.0

def magnitude(v):
	return math.sqrt(sum(v[i]*v[i] for i in range(len(v))))

def normalize(v):
	vmag=magnitude(v)
	return [v[i]/vmag for i in range(len(v))]

def vect_mu(v_list, l):
	v=np.array(v_list)
	vsum=np.ndarray.sum(v,axis=0)
	return normalize(vsum)

def getgenexpr(fname, p, tp, dt):
	f=open(fname, "r")
	p_list=[]
	tp_list=[]
	dict={}

	if dt=="m": pseudo=0.00001
	elif dt=="n": pseudo=1

	for lidx, line in enumerate(f):
		# logging sample ids and timepoints
		if lidx==0:
			tok=line.strip().split("\t")
			for i in tok[1:p*tp]:
				pid, tid=i.split("_")
				p_list.append(pid.strip())
				tp_list.append(tid.strip())
			continue
		tok=line.strip().split("\t")
		gid=tok[0]
		expr=map(float,tok[1:])
		expr=[x+pseudo for x in expr]
		dict[gid]=expr

	f.close()

	p_set=sorted(set(p_list), key=lambda x:p_list.index(x))
	tp_set=sorted(set(tp_list), key=lambda x:tp_list.index(x))

	return p_set, tp_set, dict
	
def mean_confidence_interval(data, confidence):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t._ppf((1+confidence)/2., n-1)
	return m, m-h, m+h

	
def read_clusters(fname, g, pid, tpid, dep): # cluster file, gene expression file, phenotypes, time points
	f=open(fname, "r")
	kc={}
	p=len(pid)
	tp=len(tpid)

	for line in f:
		if "x" in line: continue
		gid, kclid=line.strip().split()
		if kclid not in kc:
			kc[kclid]=Cluster(p, tp)
			if kclid in dep:
				kc[kclid].type="DEP"
				kc[kclid].sub_type="DEP"

		# archive k vector data
		if gid not in kc[kclid].k_gid:
			kc[kclid].k_gid[gid]=[]
		kc[kclid].k_gid[gid]=g[gid]

		# archive phenotype vector data
		for idx, i in enumerate(pid):
			if i not in kc[kclid].pheno:
				kc[kclid].pheno[i]=Phenotype()
			expr=[float(i) for i in g[gid]]
			v=expr[idx*tp:(idx*tp)+tp]

			kc[kclid].pheno[i].vect.append(v)
			kc[kclid].pheno[i].gid.append(gid)
	f.close()
		
	for k in sorted(kc, key=int):
		
		# build k centroid across all time points and phenotype
		k_cents=[]
		for gid in kc[k].k_gid:
			k_cents.append(kc[k].k_gid[gid])
		kc[k].k_skmeans_cent=vect_mu(k_cents, len(k_cents[0]))

		# build centroid per phenotype 
		p_cents=[]
		for i in pid:
			kc[k].gene_count=len(kc[k].pheno[i].gid)
			kc[k].pheno[i].cent=vect_mu(kc[k].pheno[i].vect, tp)
			p_cents.append(kc[k].pheno[i].cent)
		kc[k].k_cent=vect_mu(p_cents, tp)

	return kc 

def ODEP_test(kc, pid):

	for k in sorted(kc,key=int):

		if kc[k].sub_type != "NEP": continue

		i_set=[]
		for x in pid:
			cdist=[]
			for idx, i in enumerate(pid):
				if i == x: continue
				i_cent=kc[k].pheno[i].cent
				for j in pid:
					if j == x: continue

					for v in kc[k].pheno[j].vect:
						cdist.append(distance.cosine(i_cent, v))
					
			if len(cdist):
				i_set.append(cdist)

		# performing anova on i_set
		if len(i_set):
			f_val, p_val=scipy.stats.f_oneway(*i_set)
			if p_val <= 0.1**15:
				kc[k].p_val=p_val
				if kc[k].sub_type=="DEP" or kc[k].sub_type=="ODDEP":
					kc[k].sub_type="ODDEP"	
				else:
					kc[k].type="DEP"
					kc[k].sub_type="ODEP"	
			elif kc[k].type =="DEP":
				kc[k].p_val=p_val
			else:
				kc[k].type="NEP"
				kc[k].sub_type="NEP"

	return kc

# Identify clusters with constant expression pattern.
# The treshold for defining a constant expression pattern is given by the user.
# By default it is set to 1.5 fold.
def CEP_test(kc, pid, p_n, t_n):
	fc_cut=1.5
	cep_counter=0
	cep_gene_counter=0
	for k in sorted(kc,key=int):
		if kc[k].type != "NEP": continue
		max_c=max(kc[k].k_skmeans_cent)
		min_c=min(kc[k].k_skmeans_cent)
		fc=max_c/min_c
		if fc < fc_cut: 
			kc[k].type="NEP"
			kc[k].sub_type="CEP"
			cep_counter+=1
			cep_gene_counter+=kc[k].gene_count
	return kc


class Gene(object):
	def __init__(self, gid):
		self.gid=gid
		self.pheno={}

def SEP_test(kc, pid, p_n, t_n):
	gene_pool={}
	cos_dists=[]
	for k in sorted(kc,key=int):
		if kc[k].type != "NEP": continue
		cdist=[]
		for p in pid:
			for gidx, gid in enumerate(kc[k].pheno[p].gid):
				if gid not in gene_pool:
					gene_pool[gid]=Gene(gid)
				v=kc[k].pheno[p].vect[gidx]
				gene_pool[gid].pheno[p]=v
				dist=distance.cosine(kc[k].k_cent, v)
				cdist.append(dist)
				cos_dists.append(dist)

		kc[k].m_dist, kc[k].ml, kc[k].mh=mean_confidence_interval(cdist, 0.99)

	m, ml, mh = mean_confidence_interval(cos_dists, 0.99)
	ml=ml/(p_n)
	
	g_count=0
	k_count=0
	for k in sorted(kc,key=int):
		if kc[k].sub_type != "NEP": continue

		if kc[k].ml < ml:
			g_count+=kc[k].gene_count
			k_count+=1
			kc[k].type="SEP"
			kc[k].sub_type="SEP"
			
			z=(kc[k].ml - np.mean(cos_dists))/(np.std(cos_dists)/np.sqrt(len(cos_dists)))
			p_val=scipy.stats.norm.sf(abs(z))

			kc[k].p_val=p_val

	return kc

def Rescue_test(kc, pid, p_n, t_n):

	rg=set()

	for k in sorted(kc,key=int):
		if kc[k].type != "NEP": continue

		for g in kc[k].k_gid:
			v=kc[k].k_gid[g]

			min_dist=sys.float_info.max
			min_k=""
			for k_j in sorted(kc,key=int):
				if kc[k_j].type == "NEP": continue
				dist=distance.cosine(kc[k_j].k_skmeans_cent, v)

				if dist < kc[k_j].ml and dist <= min_dist:
					min_dist=dist
					min_k=k_j

			if min_k != "":
				rg.add(g)

				kc[min_k].k_gid[g]=kc[k].k_gid[g]
				for pidx, p in enumerate(kc[min_k].pheno):
					vp=v[pidx*t_n:(pidx*t_n)+t_n]
					kc[min_k].pheno[p].vect.append(vp)
					kc[min_k].pheno[p].gid.append(g)
					kc[min_k].pheno[p].mag.append(magnitude(vp))

	return kc

def DEP_test(kc, pid, p_n, t_n, dep):
	
	for k in sorted(kc,key=int):
		if kc[k].type != "NEP": continue

		#if k in dep[0]:
		if k in dep:
			kc[k].type="DEP"
			kc[k].sub_type="DEP"
			kc[k].p_val=dep[k]

	return kc

def classify_clusters(kc, pid, tpid, dep):

	p=len(pid)
	tp=len(tpid)

	k_cosdist_set=[]
	NEP_counter=0
	kc=update(kc)

	# searching for CEP clusters
	kc=CEP_test(kc, pid, p, tp)

	## searching for DEP clusters
	kc=DEP_test(kc, pid, p, tp, dep)

	# searching for ODEP clusters
	kc=ODEP_test(kc, pid)

	# searching for SEP clusters
	kc=SEP_test(kc, pid, p, tp)

	kc=update(kc)

	# rescuing genes 
	kc=Rescue_test(kc, pid, p, tp)

	kc=update(kc)

	return kc

def get_DEPs(fname):
	c={}
	f=open(fname, "r")
	for line in f:
		tok=line.strip().split()
		clid=tok[1]
		pval=float(tok[3])
		if pval < 0.05:
			c[clid]=pval
	return c


def update(kc):

	for k in sorted(kc, key=int):

		kcents=[]
		# updating whole cluster centroid
		for gid in kc[k].k_gid:
			kcents.append(kc[k].k_gid[gid])
		kc[k].k_skmeans_cent=vect_mu(kcents, kc[k].n_pheno*kc[k].n_tp)

		# updating confidence interval of whole cluster centroid
		kdists=[]
		for gid in kc[k].k_gid:
			dist=distance.cosine(kc[k].k_skmeans_cent, kc[k].k_gid[gid])
			kdists.append(dist)
		kc[k].m, kc[k].ml, kc[k].mh = mean_confidence_interval(kdists, 0.99)
		kc[k].max = max(kdists)

		pcents=[]
		for p in kc[k].pheno.values():
			if len(p.gid) < 2:
				kc[k].gene_count=len(p.gid)
				break
			p.cent=vect_mu(p.vect, kc[k].n_tp)			# updating p centroid
			pcents.append(p.cent)
			kc[k].gene_count=len(p.gid)

			# measuring cos dist. statistics
			dists=[]
			for v in p.vect:
				dists.append(distance.cosine(p.cent, v))
			p.mu_dist=np.mean(dists)
			p.sd_dist=np.std(dists)
			m,ml,mh=mean_confidence_interval(dists, 0.95)
			#p.conf_dist=mh
			p.conf_dist=max(dists)

		if kc[k].gene_count < 2: continue
		kc[k].k_cent=vect_mu(pcents, kc[k].n_tp)
	
	return kc

def print_cluster(kc):
	cls={}
	gids=set()
	for k in kc:
		k_gids=set()
		for p in kc[k].pheno.values():
			for g in p.gid:
				gids.add(g)
				k_gids.add(g)

		if kc[k].type not in cls:
			cls[kc[k].type]=0
		cls[kc[k].type]+=len(k_gids)

	for c in cls:
		print (c, cls[c], len(gids))

def get_minmax(kc):
	vs=[]
	for k in kc:
		for v in kc[k].k_skmeans_cent:
			vs.append(v)
	min_v=min(vs)
	max_v=max(vs)
	return [min_v, max_v]

def results_to_file(kc, path, pids, tps, cv, ge, dt, ge_f):

	# parameter configuration log
	param_f=open("%s/param.conf"%(path), "w")				# parameter configuration log
	param_f.write("class=")							# class
	for p in pids[0:-1]:
		param_f.write("%s,"%(p))
	param_f.write("%s\n"%(pids[-1]))
	
	param_f.write("timepoints=")						# timepoints
	for tp in tps[0:-1]:
		param_f.write("%s,"%(tp))
	param_f.write("%s\n"%(tps[-1]))

	param_f.write("expressionfile=%s\n"%(ge_f))				# expression file
	param_f.write("datatype=%s\n"%(dt))					# datatype
	param_f.close()
	
	
	# writing data to files
	depc_f=open("%s/DEP_clusters.dat"%(path), "w")				# DEP cluster output file
	depg_f=open("%s/DEP_genes.dat"%(path), "w")				# DEP gene output file
	depc_plot_f=open("%s/plots/DEP_clusters_plot.dat"%(path), "w")		# DEP cluster output file for plotting
	depg_plot_f=open("%s/plots/DEP_genes_plot.dat"%(path), "w")		# DEP gene output file
	depg_norm_plot_f=open("%s/plots/DEP_genes_norm_plot.dat"%(path), "w")	# DEP gene normalized output file
	
	sepc_f=open("%s/SEP_clusters.dat"%(path), "w")				# SEP cluster output file
	sepg_f=open("%s/SEP_genes.dat"%(path), "w")				# SEP gene file
	sepc_plot_f=open("%s/plots/SEP_clusters_plot.dat"%(path), "w")		# SEP cluster output file
	sepg_plot_f=open("%s/plots/SEP_genes_plot.dat"%(path), "w")		# SEP gene file
	sepg_norm_plot_f=open("%s/plots/SEP_genes_norm_plot.dat"%(path), "w")	# SEP gene normalized file
	
	cluster_files={"DEP":depc_f, "SEP":sepc_f}
	gene_files={"DEP":depg_f, "SEP":sepg_f}
	cluster_plot_files={"DEP":depc_plot_f, "SEP":sepc_plot_f}
	gene_plot_files={"DEP":depg_plot_f, "SEP":sepg_plot_f, "DEP_norm":depg_norm_plot_f, "SEP_norm":sepg_norm_plot_f}

	# writing headers
	for f in cluster_plot_files:
		cluster_plot_files[f].write("Clid\tSamples\tTimePoint\tValue\tType\tMin\tMax\n")
	for f in gene_plot_files:
		gene_plot_files[f].write("GeneID\tClid\tSamples\tTimePoint\tValue\tType\n")

	for f in cluster_files:
		cluster_files[f].write("Clid\tType\tGeneCount\tP-value\n")
	for f in gene_files:
		gene_files[f].write("GeneID\tClid\tType\n")

	for k in sorted(kc, key=int):

		if kc[k].type not in cluster_files: continue

		# writing cluster information
		cluster_files[kc[k].type].write("%s\t%s\t%d\t%.4f\n"%(k, kc[k].sub_type, kc[k].gene_count, kc[k].p_val))

		v=kc[k].k_skmeans_cent

		if kc[k].type == "SEP":
			max_v=0.6
			min_v=0.0
		else:
			max_v=max(v)
			min_v=min(v)

		for pidx, p in enumerate(pids):
			# writing cluster plotting data
			for tpidx, expr in enumerate(v[pidx*len(tps):(pidx*len(tps))+len(tps)]):
				tid=tps[tpidx]
				cluster_plot_files[kc[k].type].write("%s\t%s\t%s\t%f\t%s\t%f\t%f\n"%(k, p, tid, expr, kc[k].sub_type, min_v, max_v))
			
			# writing gene information
			if pidx==0:
				for g in kc[k].pheno[p].gid:
					gene_files[kc[k].type].write("%s\t%s\t%s\n"%(g, k, kc[k].sub_type))

			# writing gene plotting data
			for g in kc[k].pheno[p].gid:
				ge_norm=normalize(ge[g])
				for expr_idx in range(pidx*len(tps), (pidx*len(tps))+len(tps)):
					tidx=expr_idx%len(tps)
					tid=tps[tidx]

					gidx=kc[k].pheno[p].gid.index(g)

					gene_plot_files[kc[k].type].write("%s_%s\t%s\t%s\t%s\t%f\t%s\n"%(g, p, k, p, tid, ge[g][expr_idx], kc[k].sub_type))
					normtype="%s_norm"%(kc[k].type)
					gene_plot_files[normtype].write("%s_%s\t%s\t%s\t%s\t%f\t%s\n"%(g, p, k, p, tid, ge_norm[expr_idx], kc[k].sub_type))


	# closing files
	for f in list(cluster_files.values()) + list (gene_files.values())+ list (cluster_plot_files.values()) + list (gene_plot_files.values()):
		f.close()


def main():
	expr_f=sys.argv[1]
	kcluster=sys.argv[2]
	outdir=sys.argv[3]
	pheno=int(sys.argv[4])
	tp=int(sys.argv[5])
	datatype=sys.argv[6]
	mi_f="%s/mi_result.txt"%(outdir)

	# process gene expression file
	phenoids, timepoints, gene_expr_dict=getgenexpr(expr_f, pheno, tp, datatype)

	# archiving clusters
	DEP_cls=get_DEPs(mi_f) # archive DEPs solely based on MI results

	clusters=read_clusters(kcluster, gene_expr_dict, phenoids, timepoints, []) 

	# classify clusters further to DEP, ODDEP, ODP, SEP and NEP
	clusters=classify_clusters(clusters, phenoids, timepoints, DEP_cls) 
	clusters=update(clusters) # update clusters

	val_range=get_minmax(clusters)
	results_to_file(clusters, outdir, phenoids, timepoints, val_range, gene_expr_dict, datatype, expr_f)	# writing cluster results


if __name__=="__main__":
	main()

