###########################################################################
#Hub-Explorer_v1.0.0/exec_goatools.py
#
#    Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#    This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1.Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
###########################################################################
import os
import shutil
import pandas as pd
import numpy as np
from goatools import obo_parser
import Bio.UniProt.GOA as GOA
import gzip
from goatools.go_enrichment import GOEnrichmentStudy

#This fuunction is designed to generate annotation files for using goatools.
#I referred to `https://www.kaggle.com/code/alexandervc/gene-ontology-python-tutorial` for implementation of this function.
def prep_for_goatools(input_directory, annotation_file, gaf_file, go_file):
	'''prep_for_directory'''
	if os.path.isdir(input_directory + '/out/' +'GO_result/'): shutil.rmtree(input_directory + '/out/' +'GO_result/')
	os.mkdir(input_directory + '/out/' +'GO_result/')
	'''file_loading(gaf_file>dict)'''
	with gzip.open(input_directory + '/' + gaf_file, 'rt') as gaf_fp:
		funcs = {}
		for entry in GOA.gafiterator(gaf_fp):
			uniprot_id = entry.pop('DB_Object_ID')
			funcs[uniprot_id] = entry
	pop = funcs.keys()
	'''file_loading(gaf_file>dict>GO terms)'''
	assoc = {}
	for x in funcs:
		if x not in assoc:
			assoc[x] = set()
		assoc[x].add(str(funcs[x]['GO_ID']))
	'''file_loading'''
	obo_file = os.path.join(input_directory,go_file)
	go = obo_parser.GODag(obo_file)
	'''Annotation_file'''
	annotation = pd.read_csv(input_directory + '/' + annotation_file).loc[:,["Gene symbol","Uniprot ID"]]
	return pop, assoc, go, annotation

#Conduct goatools.
def goatools(input_directory, annotation_file, gaf_file, go_file, GCN):
	print('loading data..')
	GCN = GCN.rename(columns={"gene_2":"Gene symbol"})
	pop, assoc, go, annotation = prep_for_goatools(input_directory, annotation_file, gaf_file, go_file)
	print('proceeding with goatools..')
	for gene in np.unique(GCN["gene_1"].to_list()):
		df_corr = GCN[GCN["gene_1"] == gene]
		list = df_corr.merge(annotation, how="left", on="Gene symbol")
		list = list.dropna()["Uniprot ID"].tolist() #Have to deal with '0'
		study = {}
		for x in list:
			if x not in study:
				study[x] = set()
		study = study.keys()
		g_fdr = GOEnrichmentStudy(pop, assoc,go,propagate_counts=True,alpha=0.05,methods=["fdr_bh"]) #I used `fdr_bh` method.
		g_fdr_res = g_fdr.run_study(study) #Apply GO analysis object to targeted genes involved in GCN.
		g_fdr.wr_xlsx(f'{input_directory}/out/GO_result/{gene}.xlsx', g_fdr_res)
