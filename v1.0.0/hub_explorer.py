###########################################################################
#Hub-Explorer_v1.0.0/hub_explorer.py
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
###########################################################################


import os
import glob
import shutil
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import hub_extraction
import hub_classification

def gene_similarity_map(input_directory, annotated_clustered,n_clusters,dendrogram,min,center,max,format_type,show):
    zeileis_colors = np.array(sc.pl.palettes.godsnot_102)
    col_colors = np.array(annotated_clustered["Cluster"]).astype("<U7")
    for color in np.arange(n_clusters):
      col_colors = np.where(col_colors == str(color), zeileis_colors[[color]], col_colors)
    sns.set(style="ticks",
        font_scale=2.5,
        font='arial',
        rc = {'figure.figsize':(15,8), 'axes.linewidth':1.5})
    g = sns.clustermap(
        annotated_clustered.drop("Cluster",axis=1),
        method="ward",
        cmap="GnBu",
        center=center,
        vmin=min, vmax=max,
        col_cluster=dendrogram,
        row_cluster=dendrogram,
        col_colors=col_colors,
        row_colors=col_colors,
        xticklabels=False,
        yticklabels=False)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    g.ax_cbar.set_title('J.I.', rotation=0)
    plt.tight_layout()
    plt.savefig(f'./{input_directory}/out/result/gene_similarity_matrix.{format_type}',format=f'{format_type}',dpi=150)
    if show:
        plt.show()
    plt.close()

def prep_for_hub_visualization(input_directory):
	data_dir = glob.glob(input_directory + '/out/' + 'Module_Summary/Summary_module*.csv')
	"""matrix_processing"""
	df = pd.DataFrame(columns=["Cluster","Gene","GO","name"])
	for data in data_dir:
		df_tmp = pd.read_csv(data, usecols=[1,2,3,4])
		df = pd.concat([df,df_tmp])
	df = df.sort_values(by="Cluster", ascending=True).set_index("Cluster")
	df["Count"] = 1
	table = df.pivot_table(index="Gene", columns="name")
	table = table.fillna(0).transpose().reset_index().drop("level_0", axis=1).set_index("name")
	label = df.reset_index().loc[:,["Cluster", "Gene"]]
	label = label.drop_duplicates()
	name_label = df.reset_index().loc[:,["Cluster","name"]]
	name_label = name_label.set_index("Cluster")["name"].drop_duplicates()
	anno_table = label.merge(table.transpose().reset_index(), how="left", on="Gene")
	tmp = pd.DataFrame(columns=anno_table.columns)
	for module in np.unique(anno_table["Cluster"].tolist()):
		df_module = anno_table[anno_table["Cluster"] == module]
		df_module["gene_scale"] = anno_table.sum(axis=1)
		df_module = df_module.sort_values(by="gene_scale", ascending=False).drop("gene_scale",axis=1)
		tmp = pd.concat([tmp, df_module])
	anno_table = tmp
	anno_table = anno_table.drop("Cluster",axis=1).set_index("Gene").transpose()
	anno_table = name_label.reset_index().merge(
				anno_table.reset_index().rename(columns={"index":"name"}),
				how="left",
				on="name").sort_values(by="Cluster", ascending=True)
	frame = pd.DataFrame(columns = anno_table.columns)
	for module in np.unique(anno_table["Cluster"].tolist()):
		df_module = anno_table[anno_table["Cluster"] == module]
		df_module["GO_scale"] = anno_table.sum(axis=1)
		df_module = df_module.sort_values(by="GO_scale", ascending=True).drop("GO_scale",axis=1)
		frame = pd.concat([frame,df_module])
	anno_table = frame.drop("Cluster", axis=1).set_index("name")
	return anno_table

def hub_visualization(input_directory, anno_table):
		sns.set(style="ticks",
			font_scale=2.5,
			font='arial',
			rc = {'figure.figsize':(40,18), 'axes.linewidth':1.5})
		sns.heatmap(
			anno_table.transpose(),
			cmap="Pastel1_r",
			vmin=0, vmax=1,
			cbar=False)
		plt.xticks(fontsize=34, rotation=25, ha ='right')
		plt.yticks(fontsize=34)
		plt.tight_layout()
		plt.savefig(f'./{input_directory}/out/result/summary_table_of_hub_components.eps',format='eps',dpi=300)
		plt.savefig(f'./{input_directory}/out/result/summary_table_of_hub_components.png',format='png',dpi=300)
		plt.close()

def data_process(input_directory,n_cluster):
	hub_components = hub_extraction.extract_hub_components(input_directory)
	gene_similarity_matrix = hub_extraction.generate_gene_similarity_matrix(input_directory, hub_components)
	annotated_cluster, clustered, df_annotation, go2cluster = hub_classification.finalization(input_directory,gene_similarity_matrix,hub_components,n_cluster)
	return annotated_cluster, go2cluster

def exec(input_directory,n_cluster):
	print('data processing..')
	annotated_cluster, go2cluster = data_process(input_directory,n_cluster)
	anno_table = prep_for_hub_visualization(input_directory)
	print('visualization..')
	gene_similarity_map(
		input_directory, annotated_clustered = annotated_cluster, n_clusters=n_cluster, dendrogram=True, min=0, center=0.375, max=1.0, format_type='eps', show=False)
	gene_similarity_map(
		input_directory, annotated_clustered = annotated_cluster, n_clusters=n_cluster, dendrogram=True, min=0, center=0.375, max=1.0, format_type='png', show=False)
	hub_visualization(input_directory, anno_table)
	
	
