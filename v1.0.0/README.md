# Hub-Explorer v1.0.0
## Overview
Hub-Explorer v1.0.0 is an advanced bioinformatics toolset tailored for the identification of the cell type-specific signaling pathway based on proteomic and transcriptomic data. It is particularly adept at constructing gene co-expression networks and conducting Gene Ontology (GO) analysis. The package integrates a suite of Python scripts and a shell script to offer a comprehensive approach to genomic data analysis and visualization.

## Installation (developer's environment)
### Prerequisites:
Python 3.8.12 \
Python libraries:\
-NumPy v1.21.2\
-pandas v1.5.3\
-scipy v1.8.0\
-biopython v1.80\
-goatools v1.2.4\
-sklearn v0.0.post1\
-scanpy v1.9.2
### Steps:
1. Clone the Hub-Explorer repository or download the ZIP file.
2. Prepare for the directory for your analysis, containing 'gene_expression_matrix.csv', 'list_of_gene.csv', 'gene annotation_table.csv', 'goa_{species of target}.gaf.gz', and 'go-basic.obo'.\
*You can automatically download 'goa_{species of target}.gaf.gz' and 'go-basic.obo'. However, please take care of file version used for analysis. Both files used for Cell Genomics (2024) were deposited to demo_data directory.
3. Install the required Python libraries.

## Usage
```bash
sh Hub-Explorer_v1.0.0/Hub-Explorer.sh [OPTIONS]
Options:
    -h          	Display this help message
    -i DIRECTORY	Specify the input directry
    -m FILE     	Specify the matrix file (CSV format)
    -l FILE     	Specify the list file (CSV format)
    -a FILE     	Specify the annotation table (CSV format)
    -o FILE     	Specify the GO annotation file (Manual: 'manual', Auto:'human' for Homo Sapiens, 'mouse' for Mus Musculus)
    -k VALUE    	Specify the number of clusters
    -t VALUE		Specify the NAME of targeted gene or None
```
### File Constructions:
1. gene expression matrix file (csv format)

| Cell Types  | Gene1       | Gene2      | ••• |
| ----------- | ----------- | ---------- | --- |
| 1           | expression  | expression | ••• |
| 2           | expression  | expression | ••• |
|•••          | •••         |•••         | ••• |

2. gene list file (txt format)

|Gene1|
|-----|
|Gene2|
|Gene3|
|•••••|

3. annotation file (csv format)

| Gene symbol | Uniprot ID  | ••• (* There is no problem if you add other columns) |
| ----------- | ----------- | ---------- |
| Xkr4        | Q5GH67      |••• |
|•••          |•••          |••• |

*You should use `Gene symbol` and `Uniprot ID`!
*You can easy to get these file via scanpy object.
*In `v1.0.0`, I used Uniprot ID for GO analysis with Goatools, but other ids are also available by changing Goatools command. (e.g., Ensemble gene_id etc..)

## Authors
Noguchi Yuki (Jun Suzuki lab)\
Email: nyuhki21@gmail.com, jsuzuki@icems.kyoto-u.ac.jp

## Citation
If you utilize Hub-Explorer in your research, please cite the following publication:\
Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki, J. 2024. "In vivo CRISPR screening directly targeting testicular cells." Cell Genomics.

## License
This software is released under the MIT License. Please refer to the LICENSE file for more details.
This expanded README provides a more detailed guide, explaining the purpose and usage of each script in the Hub-Explorer package. You can use this in your GitHub repository to offer users a thorough understanding of how to install, configure, and use the software effectively.

## References
1. Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
