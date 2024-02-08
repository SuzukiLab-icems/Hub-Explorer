###########################################################################
#Hub-Explorer_v1.0.0/parser.py
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
import argparse
import glob

def parser_setting():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_directory',help='Specify the input directry')
	parser.add_argument('-m', '--matrix_file',help='Specify the matrix file (CSV format)')
	parser.add_argument('-l', '--list_file',help='Specify the list file (CSV format)')
	parser.add_argument('-a', '--annotation_file',help='Specify the annotation table (CSV format)')
	parser.add_argument('-g', '--gaf_file',help='Specify the goa_{species}.gaf')
	parser.add_argument('-o', '--go_file',help='Specify the go-basic.obo')
	parser.add_argument('-k', '--n_cluster',type=int,help='Specify the number of clusters')
	parser.add_argument('-t', '--targeted_gene',type=str,help='Specify the targeted gene')
	args = parser.parse_args()
	return args
