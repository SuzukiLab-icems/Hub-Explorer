###########################################################################
#Hub-Explorer_v1.0.0/main.py
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
import sys
import glob
import parser
import GCN_generator
import exec_goatools
import hub_explorer

def main(input_directory, matrix_file, list_file, annotation_file, gaf_file, go_file, n_cluster, targeted_gene):
	if len(glob.glob(input_directory + '/out/Module_Summary/*')) == 0:
		print('generate gene coexpression network..')
		GCN = GCN_generator.generate_gene_coexpression_network(input_directory,matrix_file,list_file,targeted_gene)
		print('execute GO analysis..')
		exec_goatools.goatools(input_directory, annotation_file, gaf_file, go_file, GCN)
	elif len(glob.glob(input_directory + '/out/Module_Summary/*')) != 0:
		print("assumed you've finished analysis. only visualization step starts proceeding..")
		print("if you want to start reanalysis, please delete the out directory.")
	print('execute hub-explorer..')
	hub_explorer.exec(input_directory,n_cluster)
	
def execute():
    args = parser.parser_setting()
    main(**vars(args))

if __name__ == "__main__":
	execute()
	print('done.')


