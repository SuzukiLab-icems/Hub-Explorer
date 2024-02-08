###########################################################################
#Hub-Explorer_v1.0.0/hub_extraction.py
#
#    Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#    This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
###########################################################################
import glob
import pandas as pd
import numpy as np
from multiprocessing import Pool, get_context, cpu_count

#This function is designed for significant GO terms (=hub components) from each gene.
def extract_hub_components(input_directory):
	genes = glob.glob(f'./{input_directory}/out/GO_result/*.xlsx')
	hub_components = pd.DataFrame(columns=["GO", "name", "p_fdr_bh", "study_items", "gene_info"])
	for gene in genes:
		gene = gene.replace('.xlsx', '').replace(f'./{input_directory}/out/GO_result/','')
		tmp = pd.read_excel(f'./{input_directory}/out/GO_result/{gene}.xlsx').loc[:,["GO", "name", "p_fdr_bh", "study_items"]]
		tmp = tmp[tmp["p_fdr_bh"] < 0.05]
		tmp["gene_info"] = gene
		hub_components = pd.concat([hub_components, tmp])
	hub_components.dropna().to_csv(input_directory + "/out/Hub_Components.csv")
	return hub_components

def jaccard_similarity(x, y):
	intersection = len(set.intersection(*[set(x), set(y)]))
	union = len(set.union(*[set(x), set(y)]))
	return intersection / float(union)

#This function is designed for calculating Jaccard Similarity Index score between genes according to hub components.
def go2go(args):
	hub_components, gene_1 = args
	jaccard_index = []
	x_GO = hub_components[hub_components["gene_info"] == gene_1]["GO"].tolist()
	for gene_2 in np.unique(hub_components["gene_info"].tolist()):
		y_GO = hub_components[hub_components["gene_info"] == gene_2]["GO"].tolist()
		try: index = jaccard_similarity(x_GO, y_GO)
		except "ZeroDivisionError": index = 'NaN'
		except "TypeError": index = 'NaN'
		jaccard_index.append([gene_1, gene_2, str(index)])
	return jaccard_index

#Multiprocessing go2go function.
#Finally, dataframe is generated.
def extract_similarity_of_genes(hub_components):
	gene_list = np.unique(hub_components["gene_info"].tolist())
	with get_context("fork").Pool(processes=cpu_count()-1) as pl:
		jaccard_index = pl.map(go2go, [(hub_components, gene) for gene in gene_list])
	gene_similarity = [item for sublist in jaccard_index for item in sublist]
	gene_similarity = pd.DataFrame(gene_similarity, columns=['gene_1','gene_2','jaccard_index'])
	return gene_similarity

#This fucntion converts gene_similarity dataframe into pivot_table (=gene_similarity_matrix).
def generate_gene_similarity_matrix(input_directory, hub_components):
	gene_similarity = extract_similarity_of_genes(hub_components)
	gene_similarity['combination'] = gene_similarity['gene_1'] + str('_vs_') + gene_similarity['gene_2'] #This script is implemented for exporting .csv file, but not required for downstream analysis.
	gene_similarity_matrix = gene_similarity.pivot_table(index="gene_1", columns="gene_2").transpose().reset_index().drop("level_0", axis=1).set_index("gene_2")
	gene_similarity_matrix = gene_similarity_matrix.transpose()
	gene_similarity_matrix.to_csv(input_directory + '/out/gene_similarity_matrix.csv')
	return gene_similarity_matrix
	

