#Install Packages
if (!require("BioocManager", quietly = T))
        install.packages("BiocManager")
BiocManager::install("AnnotationHub")

install.packages("devtools")

devtools::install_url("http://cran.r-project.org/src/contrib/Archive/**dplyr**/**dbplyr**1.3.0.tar.gz")
install.packages("igraph", type = "binary")
BiocManager::install("ensembldb")
library(AnnotationHub)
library(ensembldb)
#Connect to AnnotationHub
ah <- AnnotationHub()

ahDb <- query(ah,
              pattern = c("Mus musculus", "EnsDb"),
              ignore.case = T)

# Check versions of databases available
ahDb %>% 
        mcols()

# Acquire the latest annotation files
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")    

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

View(annotations)  

# Extract IDs for mitochondrial genes
mt <- annotations %>%
        dplyr::filter(seq_name == "MT") %>%
        dplyr::pull(gene_name)

# Number of UMIs assigned to mitochondrial genes
metadata$mtUMI <- Matrix::colSums(counts[which(rownames(counts) %in% mt),], na.rm = T)
