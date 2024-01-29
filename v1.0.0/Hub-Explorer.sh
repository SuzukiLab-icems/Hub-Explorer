#!/usr/sh
###########################################################################
#Hub-Explorer_v1.0.0/Hub-Explorer.sh
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
#
##CMD for generating figures
#Figure.5E-5G: sh ./Hub-Explorer_v1.0.0/Hub-Explorer.sh -i IO -m mtx_from_Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv -l list_modified_from_Figure.5C_summary_table_of_spermatid_specific_RD3_interactors.csv -a annotation_file_from_scanpy_var.csv -o manual -k 4 -t Rd3
#Figure.S5D-S5E: sh ./Hub-Explorer_v1.0.0/Hub-Explorer.sh -i IO -m mtx_from_Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv -l list_from_Figure.S5C.csv -a annotation_file_from_scanpy_var.csv -o manual -k 7 -t None
###########################################################################

input_directory=""
matrix_file=""
list_file=""
annotation_file=""
go_annotation=""
n_cluster=""

function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          	Display this help message
    -i DIRECTORY	Specify the input directry
    -m FILE     	Specify the matrix file (CSV format)
    -l FILE     	Specify the list file (CSV format)
    -a FILE     	Specify the annotation table (CSV format)
    -o FILE     	Specify the GO annotation file (Manual: 'manual', Auto:'human' for Homo Sapiens, 'mouse' for Mus Musculus)
    -k VALUE    	Specify the number of clusters
    -t VALUE		Specify the NAME of targeted gene or None
EOM
    exit 2
}

function argparse {
    while getopts ":i:m:l:a:o:k:t:h" opt; do
        case "$opt" in
			i) input_directory="${OPTARG}" ;;
            m) matrix_file="${OPTARG}" ;;
            l) list_file="${OPTARG}" ;;
            a) annotation_file="${OPTARG}" ;;
            o) go_annotation="${OPTARG}" ;;
            k) n_cluster="${OPTARG}" ;;
			t) targeted_gene="${OPTARG}" ;;
            h) usage ;;
            \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
            :)  echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
        esac
    done
    shift "$((OPTIND - 1))"
}

function prep_go_annotation {
    if [ "${go_annotation}" = "manual" ]; then
        echo "Annotation file of GO analysis: MANUAL"
    else
        echo "Annotation file of GO analysis: AUTO"
        if ! command -v wget > /dev/null; then
            echo "Error: wget is not installed." >&2
            exit 1
        fi
        echo "Downloading go-basic.obo..."
        wget -P ./${input_directory} http://purl.obolibrary.org/obo/go/go-basic.obo || { echo "Failed to download go-basic.obo"; exit 1; }
        case "${go_annotation}" in
            "human")
                echo "Downloading goa_human.gaf.gz..."
                wget -P ./${input_directory} http://http.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz || { echo "Failed to download goa_human.gaf.gz"; exit 1; }
                ;;
            "mouse")
                echo "Downloading goa_mouse.gaf.gz..."
                wget -P ./${input_directory} http://http.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz || { echo "Failed to download goa_mouse.gaf.gz"; exit 1; }
                ;;
            *)
                echo "Invalid option for GO annotation: '${go_annotation}'. Please set your option to 'human' or 'mouse'."
                exit 1
                ;;
        esac
    fi
}

function config {
	dir=$(pwd)"/v1.0.0" #You can frexibly change this directory.
	gaf_file=$(basename ./${input_directory}/*.gaf.gz)
	obo_file=$(basename ./${input_directory}/*.obo)
	echo " "
	echo "Start Hub-Explorer.."
	echo " "
	echo "#####Configuration#####"
	echo "version				:Hub-Explorer_v1.0.0"
	echo "package directory 	:${dir}"
	echo "input directory 		:${input_directory}"
	echo "matrix file 			:${matrix_file}"
	echo "list file 			:${list_file}"
	echo "annotation file 		:${annotation_file}"
	echo "gaf file 				:${gaf_file}"
	echo "obo file 				:${obo_file}"
	echo "Number of Cluster		:${n_cluster}"
	echo "Targeted Gene			:${targeted_gene}"
	echo " "
}

function hub_explorer {
	config
	mkdir ${input_directory}/out ${input_directory}/out/result
	python ${dir}/main.py \
		-i ${input_directory} \
		-m ${matrix_file} \
		-l ${list_file} \
		-a ${annotation_file} \
		-g ${gaf_file} \
		-o ${obo_file} \
		-k ${n_cluster} \
		-t ${targeted_gene}
}

function main {
	argparse "$@"
	prep_go_annotation
	hub_explorer
}

main "$@"

exit=0

