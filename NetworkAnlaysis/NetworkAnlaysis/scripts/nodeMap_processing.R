dir = 'OneDrive/OneDrive - Università degli Studi di Milano/ppmi/GRN/metacore_export/ppmi/MalesAndFemales/limma/females/FDR0.1/'
dir = '../metacore_export/ppmi/MalesAndFemales/limma/PDvsHC/Pre-processing/'

# Read the node file
node_file <- paste0(dir,"nodeMap.txt")
node <- read.delim(node_file, header = FALSE)
node <- node[-c(1,2,3),]
head(node)
dim(node)[1]

node <- unique(node)
head(node)
dim(node)[1]

# Read the geneslist file
hmc_geneslist <- paste0(dir, "hmc_genelist.txt")
hmc <- read.delim(hmc_geneslist, header = FALSE)
head(hmc)
dim(hmc)

# Intersect node genes symbols with the gene list
int <- intersect(node$V2, hmc$V1); cat(paste0('Intersection: ', length(int)))
node <- node[node$V2 %in% int,]
dim(node)[1]

node <- node[!duplicated(node$V2), ]
dim(node)[1]

# Salva il nuovo file senza duplicati
output_file <- "nodeMap_post.txt"
write.table(node, file = paste0(dir,output_file), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# Stampa un messaggio di completamento
cat("Duplicati rimossi. Nuovo file senza duplicati salvato in:","\n", paste0(dir,output_file), "\n")
