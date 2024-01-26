###########################################################################
#Hub-Explorer_v1.0.0/hub_classification.py
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
import shutil
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.metrics import jaccard_score

def generate_kmeans_clustered_matrix(input_directory, corr_matrix, n_cluster):
	'''Define Model'''
	Kmodel = KMeans(n_clusters=n_cluster, random_state=0, init='random')
	Kmodel.fit(corr_matrix)
	'''Generate Clustered Matrix'''
	cluster = Kmodel.labels_.tolist()
	df_annotation = pd.DataFrame(data={'Cluster': cluster}, index=corr_matrix.index).reset_index()
	df_annotation.to_csv(input_directory + '/out/labels_of_cluster.csv')
	clustered = df_annotation.sort_values(by="Cluster", ascending=False).merge(corr_matrix.reset_index(), how="left", on="gene_1").set_index("gene_1")
	clustered = df_annotation.rename(columns = {"gene_1":"index"}) \
                             .merge(clustered.drop("Cluster",axis=1).transpose().reset_index(),how="left",on="index") \
                             .sort_values(by="Cluster", ascending=False).drop("Cluster",axis=1).set_index("index")
	annotated_cluster = df_annotation.rename(columns = {"gene_1":"index"}) \
                                     .merge(clustered.reset_index(), how="left", on="index") \
                                     .sort_values(by="Cluster", ascending=False).set_index("index").rename(columns={"index":"Gene"})
	annotated_cluster.to_csv(input_directory + '/out/clustered_matrix.csv')
	return annotated_cluster, clustered, df_annotation
    
def preparation_for_finalization(input_directory, corr_matrix, hub_components, n_cluster):
    """file_preparation"""
    if os.path.isdir(input_directory + '/out/' + 'Module_GO/'): shutil.rmtree(input_directory + '/out/' + 'Module_GO/')
    os.mkdir(input_directory + '/out/' + 'Module_GO/')
    if os.path.isdir(input_directory + '/out/' + 'Module_Summary/'): shutil.rmtree(input_directory + '/out/' + 'Module_Summary/')
    os.mkdir(input_directory + '/out/' + 'Module_Summary/')
    """data_preparation"""
    annotated_cluster, clustered, df_annotation = generate_kmeans_clustered_matrix(input_directory, corr_matrix, n_cluster)
    hub_components = hub_components.rename(columns = {"gene_info":"Gene"})
    go2cluster = hub_components.merge(df_annotation.rename(columns={"gene_1":"Gene"}), how="left", on="Gene") \
                        .loc[:,["Gene","Cluster","GO","name","p_fdr_bh"]] \
                        .sort_values(by="Cluster", ascending=True)
    return annotated_cluster, clustered, df_annotation, go2cluster

def overlapped_core_extraction(input_directory, corr_matrix, hub_components, n_cluster):
    annotated_cluster, clustered, df_annotation, go2cluster= preparation_for_finalization(input_directory, corr_matrix, hub_components, n_cluster)
    for module in np.unique(go2cluster["Cluster"].tolist()):
        go2cluster_module = go2cluster[go2cluster["Cluster"] == module]
        gene_1 = go2cluster_module["Gene"].tolist()[0]
        x = go2cluster_module[go2cluster_module["Gene"] == gene_1]["GO"].tolist()
        opposite_list = go2cluster_module["Gene"].tolist()
        opposite_list.remove(gene_1)
        shared = set(x)
        """shared_GO"""
        for gene_2 in opposite_list:
            y = go2cluster_module[go2cluster_module["Gene"] == gene_2]["GO"].tolist()
            shared = shared & set(y)
        if len(shared) != 0:
            df_shared = pd.DataFrame(data=shared, columns=["GO"]) \
                          .merge(hub_components.loc[:,["GO"]],how="left",on="GO") \
                          .groupby("GO").sum() \
                          .merge(hub_components.loc[:,["GO","name"]].drop_duplicates(),how="left",on="GO").loc[:,["GO","name"]]
            df_shared.to_csv(f'./{input_directory}/out/Module_GO/shared_GO_module{module}.csv')
        else:
            print(f'Module{module}: There is no shared GO')
    go2cluster.to_csv(input_directory + '/out/' + "overlapped_hubs.csv")
    return annotated_cluster, clustered, df_annotation, go2cluster

def finalization(input_directory, corr_matrix, hub_components, n_cluster):
    annotated_cluster, clustered, df_annotation, go2cluster = overlapped_core_extraction(input_directory, corr_matrix, hub_components, n_cluster)
    go2cluster = go2cluster.loc[:,["Cluster","Gene","GO","name","p_fdr_bh"]]
    for module in np.unique(go2cluster["Cluster"].tolist()):
        go2cluster[go2cluster["Cluster"] == module].to_csv(f'./{input_directory}/out/Module_Summary/Summary_module{module}.csv')
    annotated_cluster.reset_index().loc[:,["index","Cluster"]].groupby("Cluster").count().rename(columns={"index":"count"}).to_csv(input_directory + '/out/' + 'Module_Summary/module_count_summary.csv')
    return annotated_cluster, clustered, df_annotation, go2cluster
