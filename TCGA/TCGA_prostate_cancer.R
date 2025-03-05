### Identifying the biomarker using TCGA datasets from https://www.youtube.com/watch?v=wiyl3Uv39mE&t=1219s
library("TCGAbiolinks")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("ggplot2")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("GenomicDataCommons")
library("GenomeInfoDbData")
library("keras")
library("tensorflow")
library("dplyr")

###
GDCprojects <- getGDCprojects() ### Get the list of all the projects or you can go to https://portal.gdc.cancer.gov/
TCGAbiolinks:::getProjectSummary("TCGA-PRAD")

query_TCGA <- GDCquery(
    project = "TCGA-PRAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    # workflow.type, GDC workflow type
    # legacy = FALSE,
    access = "open",
    # platform = " IlluminaGA_RNASeqV2",
    # file.type, To be used in the legacy database for some platforms, to define which file types to be used.
    # barcode,
    data.format = "TSV",
    experimental.strategy = "RNA-Seq",
    # sample.type
)

### This will download the data parallely with each chunk will download 100 files
GDCdownload(query = query_TCGA, method = "api", files.per.chunk = 100)

### This will summarized the dataset but it takes a lot of time to perform this.
tcga_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
class(tcga_data)

setwd("/diazlab/data3/.abhinav/.learning/TCGA_ML/")
saveRDS(tcga_data, file = "saveRDS/TCGA.rds")
# tcga_data <- readRDS("saveRDS/TCGA.rds")

### copying it into another variable so if it makes any changes we can identify that
sedf <- tcga_data
sedf <- sedf[, sedf@colData@listData$sample_type != "Metastatic"] ## Since it is only single sample
sedf <- sedf[, 1:550]

geneslist <- sedf@rowRanges$gene_id
samplelist <- sedf@colData@listData$sample

expr <- sedf@assays@data@listData$tpm_unstrand

rownames(expr) <- geneslist
colnames(expr) <- sedf@colData@listData$sample

drop <- apply(expr, 1, min) > 2000
expr_filtered <- expr[!drop, ]

drop2 <- apply(expr_filtered, 1, max) < 100
expr_filtered <- expr_filtered[!drop2, ]

dim(expr_filtered)

rm(expr)
# rm(tcga_data)
gc()

#### Dividing the data into train and test
### Machine learning in Python

train_data <- unlist(t(expr_filtered))
train_data <- log2(train_data+0.5)
#train_data <- train_data %>% normalize() %>% copydata()
dim(train_data) <- dim(t(expr_filtered))

hist(train_data[1,])
hist(train_data[2,])
hist(train_data[3,])
dim(train_data)

pca.gse <- PCA(train_data,graph = FALSE)
fviz_pca_ind(pca.gse, geom = "point", col.ind = sedf@colData@listData[["secondary_gleason_grade"]])

table(sedf@colData@listData[["secondary_gleason_grade"]])

```{r}
train_label <- sedf@colData@listData[["secondary_gleason_grade"]]
train_label <- train_label %>% as.factor() %>% as.numeric()
train_label <- train_label - 1
dim(train_label) <- c(dim(expr_filtered)[2], 1)

train_label <- to_categorical(train_label, num_classes = 3)

pheatmap(train_label,cluster_cols = FALSE,cluster_rows = FALSE)
```


# Construcing Neural Network 
rm(model)
dim(train_data)[2]


NNarray <- c(4096,2048, 1024, 512) 

model <- keras_model_sequential() %>%
  layer_dense(
    units = NNarray[1],
    activation = "relu",
    input_shape = dim(train_data)[2]
  ) %>%
  layer_dense(units = NNarray[2], activation = "relu") %>%
  layer_dense(units = NNarray[3], activation = "relu") %>%
  layer_dense(units = NNarray[4], activation = "relu") %>%
  layer_dense(units = dim(train_label)[2], activation = "sigmoid")

model %>% compile(optimizer = 'sgd', loss = "binary_crossentropy", metrics = c('accuracy'))

history <- model %>%
  fit(
    x = train_data,
    y = train_label,
    epochs = 200,
    use_multiprocessing = TRUE,
    batch_size = dim(train_data)[1]/25,
    #validation_split = 0.1
  )
```


```{r}
#save_model_hdf5(model,'C:/savepath/savemodel1.hdf5')

model <- keras::load_model_hdf5('C:/savepath/savemodel1.hdf5')

evaluate(model, train_data, train_label)

```

# Extraction of Weights and Bias 
```{r}
weight <- as.matrix(model$layers[[1]]$weights[[1]])
bias   <- as.matrix(model$layers[[1]]$weights[[2]])

expr_filtered <- as.matrix(expr_filtered)

dim(weight)
dim(expr_filtered)

rownames(weight) <- rownames(expr_filtered)

hist(weight[56,])
hist(bias)

#pheatmap(weight,cluster_rows = FALSE,cluster_cols = FALSE)
```
# Isolation of Gene of Interests (GOI)
```{r}
GOI <- c()

input_data  <- train_data

samplerate <- sample(1:nrow(input_data))[1:100]

for(j in samplerate) {
#for(j in 1:5) {
  sample1 <- train_data[j, ]
  
  total_weights <- weight * sample1
  total_weights_bias <- colSums(total_weights) + bias
  total_weights_bias <- as.data.frame(total_weights_bias)
  total_weights_bias <- cbind(seq(1, NNarray[1]), total_weights_bias)
  
  top_nodes <- as.data.frame(total_weights_bias[total_weights_bias[,2] > 1, ])
  
  goodnodes <- top_nodes$V1
  genes_to_goodnodes <- as.data.frame(total_weights[, goodnodes])
  
  rm(sample1)
  rm(total_weights)
  rm(total_weights_bias)
  rm(top_nodes)
  
  if(ncol(genes_to_goodnodes)>0){
    for(i in 1:ncol(genes_to_goodnodes)){
      temp <- genes_to_goodnodes[,i]
      names(temp) <- rownames(genes_to_goodnodes)
      temp <- sort(temp,decreasing = TRUE)
      temp <- temp[1:10]
      temp <- names(temp)
      
      GOI <- c(GOI, temp)
      rm(temp)
      }
    }
  
  print(paste("j=",j))
}

GOI_table <- as.data.frame(table(GOI))
GOI_table <- GOI_table %>% arrange(desc(Freq))
GOI_table <- GOI_table[GOI_table$Freq>1,]
GOI_list <- unique(GOI)
write.csv(GOI_table, file = "../Output/GOI_LIST.csv")
```

# Converting EnsembleID to EntrezID (GOI)
```{r}
gc()
genelist <- sub("[.][0-9]*", "", GOI_list)

library(clusterProfiler)
library(org.Hs.eg.db)

new_genelist <- bitr(
  genelist,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

nrow(new_genelist)
head(new_genelist)

```

# Converting EnsembleID to EntrezID (All Genes)
```{r}
fullgenelist <- rownames(expr_filtered)
fullgenelist <- sub("[.][0-9]*", "", fullgenelist)

new_fullgenelist <- bitr(
  fullgenelist,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

head(new_fullgenelist)
nrow(new_fullgenelist)

```
# Gene Set Enrichment using enrichGO() from Cluster profiler - Biological Process 
```{r}
gene <- new_genelist$ENTREZID
fullgene <- new_fullgenelist$ENTREZID

DOSEgenelist <- data(geneList, package = "DOSE")

library(org.Hs.eg.db)

ggo_BP <- groupGO(
  gene     = gene,
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  level    = 3,
  readable = TRUE
)

ego_BP <- enrichGO(
  gene       = gene,
  universe      = fullgene,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_BP, 10)
goplot(ego_BP)

ego_BP_df <- ego_BP@result
write.csv(ego_BP_df, file= "../Output/ego_BP_df.csv")
```



# Gene Set Enrichment using enrichGO() from Cluster profiler - Molecular Processes 
```{r}
ggo_MF <- groupGO(
  gene     = gene,
  OrgDb    = org.Hs.eg.db,
  ont      = "MF",
  level    = 3,
  readable = TRUE
)


ego_MF <- enrichGO(
  gene       = gene,
  universe      = fullgene,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_MF, 20)
goplot(ego_MF)

ego_MF_df <- ego_MF@result
write.csv(ego_MF_df, file= "../Output/ego_MF_df.csv")

```
# Gene Set Enrichment using enrichGO() from Cluster profiler - Cellular Categories  
```{r}


ggo_MF <- groupGO(
  gene     = gene,
  OrgDb    = org.Hs.eg.db,
  ont      = "CC",
  level    = 3,
  readable = TRUE
)


ego_CC <- enrichGO(
  gene       = gene,
  universe      = fullgene,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_CC, 20)
goplot(ego_CC)

ego_cc_df <- ego_CC@result
write.csv(ego_cc_df, file= "../Output/ego_cc_df.csv")
```

