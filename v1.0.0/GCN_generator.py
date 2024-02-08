###########################################################################
#Hub-Explorer_v1.0.0/GCN_generator.py
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
import pandas as pd
import numpy as np
import warnings
from scipy.stats import spearmanr
from multiprocessing import Pool, get_context, cpu_count
warnings.simplefilter('ignore')

#Caluculate Spearman's rank correlation coefficient between each targeted gene in the input expression matrix.
def corr2corr(args):
    expr_matrix, gene_1 = args
    opposite_list = expr_matrix.columns.to_list()
    opposite_list.remove(gene_1)
    corr = []
    for gene_2 in opposite_list:
        correlation, _ = spearmanr(expr_matrix[gene_1], expr_matrix[gene_2])
        corr.append([gene_1, gene_2, str(correlation)])
    return corr

#Conduct corr2corr function in multiprocessing manner, generating GCN.
def multi_corr2corr(expr_matrix, list_for_corr):
    try:
        expr_matrix = pd.read_csv(expr_matrix, index_col=0).drop(["groups"], axis=1) #for the matrix from scanpy.pl.paga function
    except KeyError:
        expr_matrix = pd.read_csv(expr_matrix, index_col=0) #for the matrix with Gene symbol as columns and Cell IDs as index.
    corr_list = pd.read_csv(list_for_corr, header=None)[0].tolist()
    with get_context("fork").Pool(processes=cpu_count()-1) as pool:
        corr_matrix = pool.map(corr2corr, [(expr_matrix, gene) for gene in corr_list])
    GCN = [item for sublist in corr_matrix for item in sublist]
    return GCN

#Main function for generating GCN.
#GCN is tissue (or cell type) specific GCN (it depends on your input expression matrix)
def generate_gene_coexpression_network(input_directory, matrix_file, list_file, targeted_gene):
	GCN = multi_corr2corr(
		expr_matrix = input_directory + '/' + matrix_file,
		list_for_corr = input_directory + '/' + list_file
		)
	GCN = pd.DataFrame(GCN,columns=["gene_1","gene_2","rho"])
	GCN = GCN[GCN["rho"] != "nan"]
	if targeted_gene == None:
		pass
	elif targeted_gene != None: #If you want to infer the pathway of interactors, you can input gene name used as immunoprecipitation (e.g., Rd3).
		GCN = GCN[GCN["gene_1"] != targeted_gene]
		GCN = GCN[GCN["gene_2"] != targeted_gene]
	GCN.to_csv(input_directory + '/out/result/prefiltered_GCN.csv')
	GCN["rho"] = GCN["rho"].astype("float")
	GCN = GCN[GCN["rho"] > 0.8] #Criteria for correlation is over 0.8.
	return GCN
