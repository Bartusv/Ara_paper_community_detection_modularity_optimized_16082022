###### Automated Community detection for "GENIE3" regulatory networks, Bart Verwaaijen####### 
#The script performs 3 rounds of Random Walk community detection, with modularity optimisation for step lengths 1:20. All communities are saved as .html interactive graph with target genes in blue and transcription factors in yellow, edges are directed and line thickness corresponds to the GENIE3 value which is implicated as edge weight. 
#Additionally GO-term enrichment of the sub communities against the rest of the network is performed. Further graphs are created from the clustered communities were all members of a community are contracted into a single node, the edge thickness between these community nodes corresponds 
#to the sum of all edges between these communities in the original network. Lastly a overview table is generated containing a description of each community with the top scoring BP class GO-term, modularity scores and the number of target and transcription factor genes included. A second
# table contains all genes with their community allocation and in case of Arabidopsis the TAIR annotation of these genes. 


####Set your libpath and working directory 
.libPaths( c( "/mnt/Rlibrary/" , .libPaths() ) ) ###only use if non default libPaths is used
setwd("/mnt/Network_vis/Ara_paper_communities_modularity_optimized_16082022/")


#### (GENIE3) Target matrix example, 

# A tibble: 6 × 2,400
#target_id  AT1G01010 AT1G01030 AT1G01060 AT1G01160 AT1G01250 AT1G01260 AT1G01350 AT1G01380 AT1G01520  AT1G01530 AT1G01640 AT1G01720 AT1G01780 AT1G01920 AT1G02030 AT1G02065 AT1G02110   AT1G02210 AT1G02220 AT1G02230 AT1G02250 AT1G02340 AT1G02580 AT1G03010 AT1G03040 AT1G03150 AT1G03350  AT1G03490 AT1G03650 AT1G03750 AT1G03790 AT1G03800 AT1G03840 AT1G03970 AT1G04050 AT1G04100 AT1G04240
#<chr>          <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>       <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
#1 ATMG01410 0.0000370  0.0000115 0.000153   0.000167 0.0000164 0.0000149 0.0000461 0.0000899 0.0000150 0.00000611 0.0000382 0.0000817 0.0000315 0.0000261 0.0000112 0.000287  0.000587  0.00468     0.0000264 0.000112  0.000100  0.0000353 0.000106  0.0000214 0.0000463 0.000571  0.000434  0.000377   0.0000376 0.0000734 0.000238  0.000439  0.0000722 0.0000404 0.0000155 0.0000995 0.000138 
#2 ATMG01400 0.00000526 0.0000144 0.0000859  0.000464 0.000106  0.0000158 0.000111  0.000293  0.0000246 0.0000650  0.000178  0.0000346 0.0000166 0.0000438 0.0000131 0.0000205 0.0000692 0.00306     0.0000221 0.0000154 0.0000950 0.0000694 0.000960  0.0000541 0.0000893 0.000649  0.000272  0.000403   0.0000101 0.0000381 0.000213  0.000696  0.0000444 0.000906  0.0000973 0.0000161 0.0000415
#3 ATMG01390 0.0000755  0.000115  0.0000675  0.000528 0.0000965 0.000659  0.000105  0.0000426 0.0000679 0.000189   0.000389  0.000293  0.0000668 0.0000354 0.0000207 0.000626  0.000157  0.00000107  0.0000890 0.0000773 0.0000261 0.0000636 0.000234  0.000959  0.00104   0.0000923 0.000804  0.000457   0.000266  0.000203  0.0000701 0.000339  0.0000460 0.000337  0.0000539 0.00127   0.000489 
#4 ATMG01380 0.000364   0.000159  0.000491   0.000268 0.0000339 0.000355  0.000519  0.000348  0.0000875 0.000122   0.000125  0.000106  0.000217  0.0000755 0.000744  0.0000599 0.000538  0.000000203 0.000797  0.000398  0.0000708 0.000313  0.00101   0.000132  0.000592  0.000496  0.000157  0.0000103  0.000141  0.000570  0.000307  0.00137   0.000505  0.000206  0.0000520 0.0000590 0.000323 
#5 ATMG01370 0.0000882  0.000199  0.000140   0.000506 0.000107  0.000324  0.000264  0.000179  0.000128  0.0000707  0.000116  0.000412  0.000252  0.000764  0.000225  0.000360  0.000178  0.000212    0.000213  0.000114  0.000165  0.0000906 0.0000912 0.000367  0.000291  0.000161  0.000275  0.000683   0.000364  0.000237  0.000253  0.000215  0.000477  0.000309  0.000439  0.000349  0.000788 
#6 ATMG01360 0.0000838  0.000384  0.000129   0.000899 0.0000264 0.0000547 0.0000229 0.0000237 0.0000445 0.00000670 0.000630  0.000309  0.000171  0.0000435 0.0000185 0.0000175 0.000969  0.000000120 0.000965  0.000126  0.0000292 0.00106   0.0000695 0.0000304 0.000233  0.0000935 0.0000858 0.00000728 0.000228  0.000174  0.0000265 0.0000382 0.000191  0.000184  0.000292  0.0000221 0.000101 


####TopGO .map file example, tab separated file with comma separated GO-terms

# A tibble: 28,128 × 2
#AT1G01010 `GO:0003677,GO:0005634,GO:0005634,GO:0006355,GO:0007275,GO:0016021,GO:0003700,GO:0006355`                                                                                                                                                                                                                                                                                                
#<chr>     <chr>                                                                                                                                                                                                                                                                                                                                                                                    
#1 AT1G01020 GO:0003674,GO:0005739,GO:0005739,GO:0005739,GO:0005739,GO:0009507,GO:0009507,GO:0016020,GO:0016021,GO:0005783,GO:0006665,GO:0016125                                                                                                                                                                                                                                                      
#2 AT1G01030 GO:0003677,GO:0003677,GO:0005634,GO:0005634,GO:0005634,GO:0005634,GO:0003700,GO:0006355,GO:0009908,GO:0048366,GO:1901371                                                                                                                                                                                                                                                                 
#3 AT1G01040 GO:0003677,GO:0003723,GO:0004386,GO:0004525,GO:0004525,GO:0005524,GO:0005524,GO:0005737,GO:0005737,GO:0008026,GO:0010098,GO:0016075,GO:0046872,GO:0000911,GO:0006396,GO:0003725,GO:0005515,GO:0005515,GO:0005515,GO:0005515,GO:0005634,GO:0035279,GO:0031053,GO:0005515,GO:0010267,GO:0035196,GO:0009616,GO:0005634,GO:0010445,GO:0010599,GO:0005515,GO:2000034,GO:0048317,GO:0005515,GO…
#4 AT1G01046 GO:0003674,GO:0005575,GO:0008150     

####TAIR anno example 
# A tibble: 41,671 × 5
#Model_name  Type           Short_description                                 Curator_summary                       Computational_d…
#<chr>       <chr>          <chr>                                             <chr>                                 <chr>           
#1 AT1G01010.1 protein_coding NAC domain containing protein 1                   NA                                  NAC domain cont…
#2 AT1G01020.1 protein_coding Arv1-like protein                                 NA                                  ARV1; CONTAINS …
#3 AT1G01020.2 protein_coding Arv1-like protein                                 NA                                  ARV1; CONTAINS …
#4 AT1G01030.1 protein_coding AP2/B3-like transcriptional factor family protein NA                                  NGATHA3 (NGA3);…




# create empty dataframes for the script, if for some reason the skript needs to be run again remember to newly create these as they need to be empty
Membership_list <- matrix(nrow=0, ncol = 5)
Community_summary<- matrix(nrow=0, ncol = 16)
colnames(Community_summary)<- 1:16
colnames(Membership_list) <- 1:5
Community_summary_temp <- Community_summary
Community_summary_temp

#### install the needed libraries 
#For tidyverse/svglite first intall in shell terminal "sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev r-cran-systemfonts"
#install.packages("tidyverse")
#install.packages("igraph")
#install.packages("reshape2")
#install.packages("visNetwork")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("topGO")
#install.packages("doParallel")
#install.packages("svglite")

### load libraries 
library(visNetwork)
library(tidyverse)
library(igraph)
library(reshape2)
library(topGO)
source("topGO_alt.GenTable.R") #needed for low p-values written by D.Wulf
library(doParallel)

#####Load input data, change to your own  
target.matrix <- readRDS(file = "Ath_v1_1_7_5_mio_reads.Rds") #Target matrix from GENIE3
geneID2GO_Crein <- readMappings("gen2Go_At_2020401.map") #mapfile for TopGO
Tair_anno<- read_tsv(file = "TAIR10_functional_descriptions.txt") #TAIR annotation file 

####Tune-able parameters 
set.seed(42) #### only used to create the initial seed in the script 
#steps_1 <- 5  #### the steps parameter for each round of random walks determines the number of steps the random walkers take and influences the communities created, max number of useful steps is determined by "depth"of the network 
#steps_2 <- 3
#steps_3 <- 2

#Clean the TAIR annotation of splice variant annotations when needed, usually it is
Tair_anno$Model_name <- gsub("\\..*","",Tair_anno$Model_name)
Tair_anno <- Tair_anno[match(unique(Tair_anno$Model_name), Tair_anno$Model_name),]




# filter the target matrix by strength, filter was removed for this run!! 


target.matrix[target.matrix<0.005]<-NA
test.matrix_1 <- target.matrix


#######Warning!!! Target_ids in test.matrix_1 should be in first column instead of rownames, otherwise fix
#test.matrix_col_1 <- row.names(test.matrix_1)
#test.matrix_1 <- cbind(test.matrix_1, test.matrix_col_1)

####Create edge table for igraph 
test.matrix_long_1 <- melt(test.matrix_1)
Edge_table <- test.matrix_long_1 %>% filter(!is.na(value))
Edge_table <- data.frame(Edge_table[,2], Edge_table[,1], Edge_table[,3], stringsAsFactors = FALSE)
colnames(Edge_table) <- c("TF_factor", "Candidate_gene", "Edge_strenght")

####Create nodes table for igraph 
TFs <- as.matrix(colnames(target.matrix[,-1]))

####remove free nodes, by looking for nodes without edges in the edge table 

tempx <- rbind(as.matrix(Edge_table$Candidate_gene), as.matrix(Edge_table$TF_factor))
tablex <- as.data.frame(table(tempx))
CGs <- as.matrix(tablex$tempx)
x <- matrix("GC", ncol= 1, nrow = length(CGs))
submatrix2 <- cbind(CGs, x)
colnames(submatrix2) <- c("Gene_id", "TF_or_CG")
tibble2 <- as_tibble(submatrix2)
Main_matrix <- mutate(tibble2, TF_or_CG = ifelse(Gene_id %in% as.matrix(TFs), "TF", "CG"))

#Create igraph data
igraph_edges <- Edge_table
colnames(igraph_edges) <- c("from", "to", "weight")
igraph_nodes <- Main_matrix
colnames(igraph_nodes) <- c("id", "TF_or_CG")

#### create igraph network object
net <- graph_from_data_frame(d=igraph_edges, vertices = igraph_nodes, directed=T)
saveRDS(net, file = "Main_network_file.rds")



##Walktrap Community detection level 1
weights <- as.matrix(Edge_table[,3])
###modification 


 if (length(list.files(path = "./", pattern = "cluster_list.rds")) == 1) {
   
  cluster_list <- readRDS(file = "cluster_list.rds")
   
 } else {
   
   cl <- makePSOCKcluster(20)
   #clusterEvalQ(cl, .libPaths("/mnt/Rlibrary/"))
   registerDoParallel(cl)
   cluster_list <- foreach(steps=1:20) %dopar% igraph::cluster_walktrap(net, steps = steps)
   stopCluster(cl)
   
 }
 
 Modularities <- matrix(ncol=2, nrow = 0)
 
 
 for (steps in 1:20) {
   moduls <- data.frame(steps, modularity(cluster_list[[steps]]))
   Modularities <- rbind(Modularities, moduls)
   rm(moduls)
 } 
 
 colnames(Modularities) <- c("Walktrap_steps", "Modularity")
 c <- cluster_list[[Modularities[which.max(Modularities$Modularity),1]]]


saveRDS(cluster_list, file = "cluster_list.rds")

#c <- cluster_walktrap(net, steps = steps_1)


community_1 <-  ggplot(Modularities, aes(y=Modularity, x= Walktrap_steps)) + geom_line() + ggtitle("Community_Main") + 
geom_vline(xintercept = Modularities[which.max(Modularities$Modularity),1], color = "red", linetype="dashed")

ggsave("community_main.svg", community_1)
write.csv(file=paste0("Modularities_communities_main.csv"), Modularities)



##old code#### c <- cluster_walktrap(net, weights = weights, steps = steps_1) ##play with this value

##### See how many members per Community and Modularity 
table(c$membership)
modularity(c)

##### load memberhsip list 
vb <- data.frame(c$names, c$membership,0,0, vertex_attr(net,"TF_or_CG") )
colnames(vb) <- 1:5
Membership_list <- rbind(Membership_list, vb)
rm(vb)

####pull out first level sub graphs 

for (y in 1 : max(c$membership)) {
  #y <- 1
  Sub_community = which(c$membership == y)
  Keep = V(net)[Sub_community]
  net2 <- induced_subgraph(net, Keep)
  col <- vertex_attr(net2, "TF_or_CG")
  net2 <- set.vertex.attribute(net2, name = "group",value = col)
  
  test.visn <- toVisNetworkData(net2)
  
  test.visn$edges$value <- test.visn$edges$weight
  
  if (count(test.visn$edges)==0) {
    df2 <- matrix(nrow=1, ncol = 16)
    colnames(df2) <- 1:16
    Community_summary_temp <- rbind(Community_summary_temp, df2)
    rm(df2)
    
    # test expression
    next 
    
  }
  
  
  graph <- visNetwork(test.visn$nodes, test.visn$edges, height = 1080, width = 1920) %>%
    visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow"))) %>%
    visGroups(groupname = "TF", color = list(background = "#f0b51f", border = "#cf7113")) %>%
    visGroups(groupname = "CG", color = list( background = "#00aedb" , border = "#15679e")) %>%
    visNodes(borderWidth = 2)
  
  saveRDS(net2, file = paste0("Community_network_file_", y, ".rds"))
  
  
  ##Walktrap Community detection level 2
  cl <- makePSOCKcluster(20)
  #clusterEvalQ(cl, .libPaths("/mnt/Rlibrary/"))
  registerDoParallel(cl)
  cluster_list <- foreach(steps=1:20) %dopar% igraph::cluster_walktrap(net2, steps = steps)
  #stopCluster(cl)
  
  
  
  Modularities <- matrix(ncol=2, nrow = 0)
  
  
  for (steps in 1:20) {
  
  
    moduls <- data.frame(steps, modularity(cluster_list[[steps]]))
    Modularities <- rbind(Modularities, moduls)
    rm(moduls)
  } 
  
  colnames(Modularities) <- c("Walktrap_steps", "Modularity")
  c_sub <- cluster_list[[Modularities[which.max(Modularities$Modularity),1]]]
  
  #c_sub <- cluster_walktrap(net2, steps = steps_2)
  
  
  
  community_1 <-  ggplot(Modularities, aes(y=Modularity, x= Walktrap_steps)) + geom_line() + ggtitle(paste0("Community_", y)) + 
  geom_vline(xintercept = Modularities[which.max(Modularities$Modularity),1], color = "red", linetype="dashed")
  
  ggsave(filename = paste0("community_", y,".svg"), community_1)
  write.csv(file=paste0("Modularities_communities_", y, ".csv"), Modularities)
  
  
  
  
  #modularity(c_sub)
  #table(c_sub$membership)
  
  net2_nodes <- as.data.frame(cbind(vertex_attr(net2,"name"), vertex_attr(net2,"TF_or_CG")))
  colnames(net2_nodes) <- c("Ath_id", "TF_or_CG")
  net2_nodes <- merge(net2_nodes, Tair_anno, by.x = "Ath_id", by.y = "Model_name" )
  length(unique(net2_nodes$Ath_id))
  Tair_names <- unique(Tair_anno$Model_name) %in% unique(net2_nodes$Ath_id) 
  full_list_TRUE_FALSE = factor(as.integer(Tair_names))
  
  names <- unique(Tair_anno$Model_name)
  names(full_list_TRUE_FALSE)<- names
  
  for ( j in c("BP","MF", "CC")){
    #j = c("BP")
    controlvstreatment_up_GO <- new("topGOdata",
                                    description = "controlvstreatment_up",
                                    ontology = j,
                                    allGenes = full_list_TRUE_FALSE,
                                    nodeSize = 10,
                                    annot = annFUN.gene2GO,
                                    gene2GO = geneID2GO_Crein)
    
    #Statistische Test (Fischer)
    Fisher_controlvstreatment_up <- runTest(controlvstreatment_up_GO, 
                                            algorithm =
                                              "classic", statistic = "fisher")
    
    #this extracts the GO terms, you can change that by changing topNodes
    TableGO_Bor_cluster_18 <- topGO_altGenTable(controlvstreatment_up_GO,
                                                classicFisher = Fisher_controlvstreatment_up,
                                                topNodes= length(Fisher_controlvstreatment_up@score), 
                                                rawpvalue = 1e-1500)
    TableGO_Bor_cluster_18$q_value <- p.adjust(TableGO_Bor_cluster_18$classicFisher, method='BY')
    
    write.table(TableGO_Bor_cluster_18, file = paste0("TopGO_network_community_",y,"_",0,"_",0,"_",j,".csv"))
    level <- 1
    df2 <-    data.frame(y,0, 0,modularity(c),"NA","NA", table(net2_nodes$TF_or_CG)[1],table(net2_nodes$TF_or_CG)[2], j,TableGO_Bor_cluster_18[1,])
    colnames(df2) <- 1:16
    Community_summary <- rbind(Community_summary,df2)
    write_csv(Community_summary, file = "Community_overview_temp.csv")
    
    rm(TableGO_Bor_cluster_18)
    if (j == "BP") {
      
      visSave(graph, file= paste0("Network_community_",y,"_",0,"_",0,".html"), selfcontained = TRUE, background = "white")
      
      row.names(df2) <- y
      Community_summary_temp <- rbind(Community_summary_temp, df2)
      rm(df2)
    }
    
  } #To start loop until here use "}}" use "}" to start the whole loop
  
  Community_summary_temp_1 <- matrix(nrow=0, ncol = 16)
  colnames(Community_summary_temp_1)<- 1:16
  
  v <- data.frame(c_sub$names, y ,c_sub$membership,0, vertex_attr(net2,"TF_or_CG") )
  colnames(v) <- 1:5
  Membership_list <- rbind(Membership_list, v)
  rm(v)     
  
  for (k in 1 : max(c_sub$membership)) {
    #k <- 3
    Small_sub = which(c_sub$membership == k)
    Keep_sub = V(net2)[Small_sub]
    net3 <- induced_subgraph(net2, Keep_sub)
    col <- vertex_attr(net3, "TF_or_CG")
    net3 <- set.vertex.attribute(net3, name = "group",value = col)
    
    test.visnb <- toVisNetworkData(net3)
    test.visnb$edges$value <- test.visnb$edges$weight
    
    if (count(test.visnb$edges)==0)
    {
      df2 <- matrix(nrow=1, ncol = 16)
      colnames(df2) <- 1:16
      Community_summary_temp_1 <- rbind(Community_summary_temp_1, df2)
      rm(df2)
      # test expression
      next 
      
    }
    
    graph <- visNetwork(test.visnb$nodes, test.visnb$edges, height = 1080, width = 1920) %>%
      visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
      visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow"))) %>%
      visGroups(groupname = "TF", color = list(background = "#f0b51f", border = "#cf7113")) %>%
      visGroups(groupname = "CG", color = list( background = "#00aedb" , border = "#15679e")) %>%
      visNodes(borderWidth = 2)
    
    saveRDS(net3, file = paste0("Community_network_file_", y,"_",k, ".rds"))
    
    ##Walktrap Community detection level 3
    #cl <- makePSOCKcluster(20)
    registerDoParallel(cl)
    cluster_list <- foreach(steps=1:20) %dopar% igraph::cluster_walktrap(net3, steps = steps)
    #stopCluster(cl)



    Modularities <- matrix(ncol=2, nrow = 0)


    for (steps in 1:20) {


     moduls <- data.frame(steps, modularity(cluster_list[[steps]]))
      Modularities <- rbind(Modularities, moduls)
      rm(moduls)
    }

    colnames(Modularities) <- c("Walktrap_steps", "Modularity")
    c_subb <- cluster_list[[Modularities[which.max(Modularities$Modularity),1]]]

    #c_subb <- cluster_walktrap(net3, steps = steps_3)
    
    
    community_1 <-  ggplot(Modularities, aes(y=Modularity, x= Walktrap_steps)) + geom_line() + ggtitle(paste0("Community_", y, "_",k)) + 
    geom_vline(xintercept = Modularities[which.max(Modularities$Modularity),1], color = "red", linetype="dashed")
    
     ggsave(filename = paste0("community_", y,"_",k,".svg"), community_1)
     write.csv(file=paste0("Modularities_communities_", y,"_",k, ".csv"), Modularities)
    
    
    
    
    
    
    
    net3_nodes <- as.data.frame(cbind(vertex_attr(net3,"name"), vertex_attr(net3,"TF_or_CG")))
    colnames(net3_nodes) <- c("Ath_id", "TF_or_CG")
    net3_nodes <- merge(net3_nodes, Tair_anno, by.x = "Ath_id", by.y = "Model_name" )
    length(unique(net3_nodes$Ath_id))
    Tair_names <- unique(Tair_anno$Model_name) %in% unique(net3_nodes$Ath_id) 
    full_list_TRUE_FALSE = factor(as.integer(Tair_names))
    
    names <- unique(Tair_anno$Model_name)
    names(full_list_TRUE_FALSE)<- names
    
    for ( j in c("BP","MF", "CC")){
      #j = c("BP")
      controlvstreatment_up_GO <- new("topGOdata",
                                      description = "controlvstreatment_up",
                                      ontology = j,
                                      allGenes = full_list_TRUE_FALSE,
                                      nodeSize = 10,
                                      annot = annFUN.gene2GO,
                                      gene2GO = geneID2GO_Crein)
      
      #Statistische Test (Fischer)
      
      Fisher_controlvstreatment_up <- runTest(controlvstreatment_up_GO, 
                                              algorithm =
                                                "classic", statistic = "fisher")
      
      #this extracts the top100 GO terms, you can change that by changing topNodes
      
      TableGO_Bor_cluster_18 <- topGO_altGenTable(controlvstreatment_up_GO,
                                                  classicFisher = Fisher_controlvstreatment_up,
                                                  topNodes= length(Fisher_controlvstreatment_up@score), 
                                                  rawpvalue = 1e-1500)
      TableGO_Bor_cluster_18$q_value <- p.adjust(TableGO_Bor_cluster_18$classicFisher, method='BY')
      
      write.table(TableGO_Bor_cluster_18, file = paste0("TopGO_network_community_",y,"_",k,"_",0,"_",j,".csv"))
      level <- 2
      df2 <-    data.frame(y,k, 0,modularity(c),modularity(c_sub),"NA", table(net3_nodes$TF_or_CG)[1],table(net3_nodes$TF_or_CG)[2], j,TableGO_Bor_cluster_18[1,])
      colnames(df2) <- 1:16
      Community_summary<- rbind(Community_summary,df2)
      
      rm(TableGO_Bor_cluster_18)
      if (j == "BP") {
        
        visSave(graph, file= paste0("Network_community_",y,"_",k,"_",0,".html"), selfcontained = TRUE, background = "white")
        
        row.names(df2) <- k
        Community_summary_temp_1 <- rbind(Community_summary_temp_1, df2)
        rm(df2)
      }
      
    }### To start loop until here "}}}" use "}" to start the whole loop
    
    Community_summary_temp_2 <- matrix(nrow=0, ncol = 16)
    colnames(Community_summary_temp_2)<- 1:16
    
    v <- data.frame(c_subb$names, y ,k, c_subb$membership, vertex_attr(net3,"TF_or_CG") )
    colnames(v) <- 1:5
    Membership_list <- rbind(Membership_list, v)
    rm(v)     
    
    for (i in 1 : max(c_subb$membership)) {
      #i <- 9
      Small_sub = which(c_subb$membership == i)
      Keep_sub = V(net3)[Small_sub]
      net4 <- induced_subgraph(net3, Keep_sub)
      col <- vertex_attr(net4, "TF_or_CG")
      col1 <- vertex_attr(net4, "name")
      Node_ids <- cbind(col1, col)
      colnames(Node_ids) <- c("Ath_id", "TF_or_CG") 
      net4 <- set.vertex.attribute(net4, name = "group",value = col)
      
      
      
      ## convert to VisNetwork-list
      test.visnbb <- toVisNetworkData(net4)
      ## copy column "weight" to new column "value" in list "edges"
      test.visnbb$edges$value <- test.visnbb$edges$weight
      
      if (count(test.visnbb$edges)==0)
      {
        df2 <- matrix(nrow=1, ncol = 16)
        colnames(df2) <- 1:16
        Community_summary_temp_2 <- rbind(Community_summary_temp_2, df2)
        rm(df2)
        # test expression
        next 
        
      }
      
      graph <- visNetwork(test.visnbb$nodes, test.visnbb$edges, height = 1080, width = 1920) %>%
        visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
        visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow"))) %>%
        visGroups(groupname = "TF", color = list(background = "#f0b51f", border = "#cf7113")) %>%
        visGroups(groupname = "CG", color = list( background = "#00aedb" , border = "#15679e")) %>%
        visNodes(borderWidth = 2)
      
      saveRDS(net4, file = paste0("Community_network_file_", y,"_",k, "_",i, ".rds"))
      
      xv <- cbind(vertex_attr(net4,"name"), vertex_attr(net4,"TF_or_CG"))
      
      net4_nodes <- as.data.frame(cbind(vertex_attr(net4,"name"), vertex_attr(net4,"TF_or_CG")))
      colnames(net4_nodes) <- c("Ath_id", "TF_or_CG")
      net4_nodes <- merge(net4_nodes, Tair_anno, by.x = "Ath_id", by.y = "Model_name" )
      length(unique(net4_nodes$Ath_id))
      
      ### TOPGO
      
      #prepare TRUE_FALSE matrix
      Tair_names <- unique(Tair_anno$Model_name) %in% unique(net4_nodes$Ath_id) 
      full_list_TRUE_FALSE = factor(as.integer(Tair_names))
      names <- unique(Tair_anno$Model_name)
      names(full_list_TRUE_FALSE)<- names
      
      for ( j in c("BP","MF", "CC")){
        #j = c("BP")
        controlvstreatment_up_GO <- new("topGOdata",
                                        description = "controlvstreatment_up",
                                        ontology = j,
                                        allGenes = full_list_TRUE_FALSE,
                                        nodeSize = 10,
                                        annot = annFUN.gene2GO,
                                        gene2GO = geneID2GO_Crein)
        
        #Statistische Test (Fischer)
        Fisher_controlvstreatment_up <- runTest(controlvstreatment_up_GO, 
                                                algorithm =
                                                  "classic", statistic = "fisher")
        
        #this extracts the top100 GO terms, you can change that by changing topNodes
        TableGO_Bor_cluster_18 <- topGO_altGenTable(controlvstreatment_up_GO,
                                                    classicFisher = Fisher_controlvstreatment_up,
                                                    topNodes= length(Fisher_controlvstreatment_up@score), 
                                                    rawpvalue = 1e-1500)
        TableGO_Bor_cluster_18$q_value <- p.adjust(TableGO_Bor_cluster_18$classicFisher, method='BY')
        
        write.table(TableGO_Bor_cluster_18, file = paste0("TopGO_network_community_",y,"_",k,"_",i,"_",j,".csv"))
        level <- 3
        df2 <-    data.frame(y,k, i,modularity(c),modularity(c_sub),modularity(c_subb), table(net4_nodes$TF_or_CG)[1],table(net4_nodes$TF_or_CG)[2], j,TableGO_Bor_cluster_18[1,])
        colnames(df2) <- 1:16
        Community_summary<- rbind(Community_summary, df2)
        
        rm(TableGO_Bor_cluster_18)
        
        if (j == "BP") {
          
          visSave(graph, file= paste0("Network_community_",y,"_",k,"_",i,".html"), selfcontained = TRUE, background = "white")
          
          row.names(df2) <- i
          Community_summary_temp_2 <- rbind(Community_summary_temp_2, df2)
          rm(df2)
          
        }
        
      }}
    
    ####Plot the contracted communities as new network round 3
    x <- contract(net3, membership(c_subb), vertex.attr.comb=toString)
    x2 <- simplify(x)
    col <- vertex_attr(x2, "TF_or_CG")
    x2 <- set.vertex.attribute(x2, name = "group",value = col)
    
    #Convert to Visnetwork() object 
    test.visn <- toVisNetworkData(x2)
    
    if (count(test.visn$edges)==0)
    {
      
      next 
      
    }
    
    ## copy column "weight" to new column "value" in list "edges"
    test.visn$edges$value <- test.visn$edges$weight
    test.visn$nodes$label <- Community_summary_temp_2[,11]
    
    #Remove Communities with just a single member 
    remove <- as.matrix(table(c_subb$membership)[(table(c_subb$membership))==1])
    
    if(length(remove)==0){
      
      
      
      Community_Graph <- visNetwork(test.visn$nodes, test.visn$edges, height = 1080, width = 1920 ) %>%
        visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))
    } else {
      Community_Graph <- visNetwork(test.visn$nodes[-(as.double(row.names(remove))),], test.visn$edges, height = 1080, width = 1920 ) %>%
        visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))  
    }
    visSave(Community_Graph, file= paste0("Community_graph_",y,"_",k, ".html"), selfcontained = TRUE, background = "white") 
    
  }
  
  ####Plot the contracted communities as new network round 2
  x <- contract(net2, membership(c_sub), vertex.attr.comb=toString)
  x2 <- simplify(x)
  col <- vertex_attr(x2, "TF_or_CG")
  x2 <- set.vertex.attribute(x2, name = "group",value = col)
  
  #Convert to Visnetwork() object 
  test.visn <- toVisNetworkData(x2)
  
  if (count(test.visn$edges)==0)
  {
    
    next 
    
  }
  
  ## copy column "weight" to new column "value" in list "edges"
  test.visn$edges$value <- test.visn$edges$weight
  test.visn$nodes$label <- Community_summary_temp_1[,11]
  
  #Remove Communities with just a single member 
  remove <- as.matrix(table(c_sub$membership)[(table(c_sub$membership))==1])
  
  if (length(remove)==0) {
    
    Community_Graph <- visNetwork(test.visn$nodes, test.visn$edges, height = 1080, width = 1920 ) %>%
      visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))
    
  } else {
    
    Community_Graph <- visNetwork(test.visn$nodes[-(as.double(row.names(remove))),], test.visn$edges, height = 1080, width = 1920 ) %>%
      visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))
    
  }  
  visSave(Community_Graph, file= paste0("Community_graph_",y, ".html"), selfcontained = TRUE, background = "white") 
}

#### Save the Community_summary file
colnames(Community_summary) <- c("Level 1",	"Level 2", "Level 3",	"modularity 1", "modularity 2",	"modularity 3", "Target genes",	"TFs",	"GO-class",	"Term",	"Description",	"Annotated",	"Significant",	"Expected",	"p",	"p-adj")
rownames(Community_summary)<-NULL
write.csv(Community_summary, file = paste0("Communities_overview_", date(), ".csv"))

##### Clean up, attach TAIR anno and save Communities_memberships
All_memberships <- merge(Membership_list, Tair_anno, by.x = "1", by.y = "Model_name" )
All_memberships <- as_tibble(All_memberships)
All_memberships <- All_memberships %>% group_by(All_memberships$`1`) %>% top_n(1, `2`+`3`+`4`)
colnames(All_memberships)[1:5] <- c("Gene_id", "level 1", "level 2", "level 3", "TF_or_CG" )
All_memberships <- mutate(All_memberships, TF_or_CG = ifelse(Gene_id %in% as.matrix(TFs), "TF", "CG"))
All_memberships <- arrange(All_memberships, `level 1`, `level 2`,`level 3`,`TF_or_CG`)
All_memberships <- All_memberships[,-10]
write.csv(All_memberships , file = paste0("Communities_memberships_", date(), ".csv"))

####Plot the contracted communities as new network level 1
x <- contract(net, membership(c), vertex.attr.comb=toString)
x2 <- simplify(x)
col <- vertex_attr(x2, "TF_or_CG")
x2 <- set.vertex.attribute(x2, name = "group",value = col)

#Convert to Visnetwork() object 
test.visn <- toVisNetworkData(x2)

## copy column "weight" to new column "value" in list "edges"
test.visn$edges$value <- test.visn$edges$weight
test.visn$nodes$label <- Community_summary_temp[,11]

#Remove Communities with just a single member 
remove <- as.matrix(table(c$membership)[(table(c$membership))==1])

if(length(remove)==0){
  
  Community_Graph <- visNetwork(test.visn$nodes, test.visn$edges, height = 1080, width = 1920) %>%
    visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))
} else {
  Community_Graph <- visNetwork(test.visn$nodes[-(as.double(row.names(remove))),], test.visn$edges, height = 1080, width = 1920) %>%
    visIgraphLayout(randomSeed = 50, "layout_with_fr" ) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visEdges(arrows = list(to = list(enabled = TRUE, type = "arrow")))
  
}
visSave(Community_Graph, file= paste0("Main_Community_graph.html"), selfcontained = TRUE, background = "white")

