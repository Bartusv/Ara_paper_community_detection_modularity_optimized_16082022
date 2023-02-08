# Ara_paper_community_detection_modularity_optimized_16082022

Automated Community detection for "GENIE3" regulatory networks, Bart Verwaaijen

The script performs 3 rounds of Random Walk community detection, with modularity optimisation for step lengths 1:20. 
All communities are saved as .html interactive graph with target genes in blue and transcription factors in yellow, 
edges are directed and line thickness corresponds to the GENIE3 value which is implicated as edge weight. Additionally
GO-term enrichment of the sub communities against the rest of the network is performed. Further graphs are created 
from the clustered communities were all members of a community are contracted into a single node, the edge thickness 
between these community nodes corresponds to the sum of all edges between these communities in the original network. 
Lastly a overview table is generated containing a description of each community with the top scoring BP class GO-term, 
modularity scores and the number of target and transcription factor genes included. A second table contains all genes 
with their community allocation and in case of Arabidopsis the TAIR annotation of these genes. 


Code for GO-term Enrichment with TopGO was adapted from Donat Wulf and Andrea Bräutigam
