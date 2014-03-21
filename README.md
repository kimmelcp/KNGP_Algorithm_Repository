KNGP_Algorithm_Repository
=========================


Overview: This repository contains the files for the KNGP algorithm. The KNGP code, written in R, can be found at the following GitHub URL: https://gist.github.com/kimmelcp/9025360. Alternatively, there is also an online app for the KNGP algorithm at the following URL: http://spark.rstudio.com/kimmelcp/KNGPApp/. Note that the online app cannot handle datasets larger than 2 GB. If more space is needed, a higher-end machine should be sought.

A brief description of the four inputs is given below.. These files should be provided in CSV format without any headers.

Link knowledge matrix: The link knowledge matrix is a n*n matrix where n is the number of nodes in the knowledge network, and an entry in it represents the link weight between the nodes specified by the row number and column number. The link weight matrix in this repository represents the similarity between all 17,691 proteins in UniProt. The actual link weights are derived by combining known protein interactions and predicted protein-protein interactions along with similarity values calculated from the Cellular Component Gene Ontology (GO). Note that this file is so large (>1 GB) that is was first compressed, so the user should un-zip the file.

Node knowledge vector: The node knowledge vector is a n dimensional vector where n is the number of nodes in the knowledge network, and an entry in it represents the node weight associated with the corresponding node. The node weights in this repository represent the number of Gene Ontology terms associated with each of the unique 17,691 proteins.

Gene identifiers: List of n protein identifiers that correspond to the nodes in the node knowledge vector. The proteins in this repository represent all 17,691 unique proteins in UniProt as of March 2011.

Root node set: A set of gene identifiers that is associated with the disease of interest and is a subset of the candidate node set. For this repository, the root nodes genes can be found in the folder "RootNodeSets". There are a total of 19 diseases represented.

In addition to these data files, the code for the KNGP algorithm is also available in the folder “KNGPCode”.


