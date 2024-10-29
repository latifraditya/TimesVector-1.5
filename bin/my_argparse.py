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

import argparse

class input_info:
	args = {}
	fin_name = ''
	fout_name = ''
	ntop = 0
	bin_num = 0
	kde = False
	outdir = ''
	cond = ''

def arg_parsing():
	info = input_info()
	parser = argparse.ArgumentParser(prog="deg_entropy.py")
	parser.add_argument('-b', '--bin', type=int, default=20, help="number of bin")
	parser.add_argument('-kde', '--kde', action='store_true', help="using kde")
	parser.add_argument('-o', '--outdir', help="output directory")
	parser.add_argument('-n', '--num', type=int, default=50, help="number of top genes")
	parser.add_argument('-i', '--input', help="txt file name of expression data", required=True)
	parser.add_argument('-c', '--condition', help="condition of each sample", required=True)

	args = vars(parser.parse_args())
	info.args = args
	#f_expr = open(args['expression'], 'r')
	#fout = open(args['output'], 'w')
	info.fin_name = args['input']
	info.cond = args['condition']
	info.outdir = args['outdir']
	info.ntop = args['num']
	info.bin_num = args['bin']
	info.kde = args['kde']
	
	return info

class fin_matrix:
	geneList = []
	input_matrix = []
	clsList = []
	numEachClass = {}
	clsCols = {}
	nClass = 0
	nSample = 0

def fin_parsing(fin_name, fcls_name):
	fmat = fin_matrix()
	fin = open(fin_name, 'r')
	fin.readline() # first line is header
	for line in fin:
		spl = line.strip().split()
		gene = spl[0]
		fmat.geneList += [gene] 
		fmat.input_matrix.append([float(x) for x in spl[1:]])
	fin.close()

	fcls = open(fcls_name, 'r')
	cls_spl = fcls.readline().strip().split()
	fmat.clsList = cls_spl
	fmat.nSample = len(cls_spl)
	for i, c in enumerate(cls_spl):
		if not c in fmat.numEachClass:
			fmat.numEachClass[c] = 0
		if not c in fmat.clsCols:
			fmat.clsCols[c] = []
		fmat.numEachClass[c] +=1
		fmat.clsCols[c].append(i)
	fmat.nClass = len(fmat.numEachClass)
	
	return fmat

