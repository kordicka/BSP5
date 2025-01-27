
#dir = 'OneDrive/OneDrive - Università degli Studi di Milano/ppmi/GRN/metacore_export/ppmi/MalesAndFemales/limma/females/FDR0.1/'
dir = '../metacore_export/ppmi/MalesAndFemales/limma/PDvsHC/Pre-processing/'

# Read the input file
input_file <- paste0(dir,"NetworkInteractions.txt")
data <- read.delim(input_file,col.names = c('Network Object "FROM"', 'Effect', 'Network Object "TO"'))
data <- data[-c(1,2),]
cat('check data: ')
head(data)
# Perform substitutions
data$Effect <- gsub("Unspecified", "--?", data$Effect)
data$Effect <- gsub("Inhibition", "--|", data$Effect)
data$Effect <- gsub("Activation", "-->", data$Effect)
data <- data[!(data$Effect == 'Technical'), ]
cat('Numer of interaction: ')
table(data$Effect)

# Stampa il dataframe risultante
cat('Data after format check: ')
head(data)

# Write the modified data to a new file
output_file <- paste0(dir,"NetworkInteraction_post.txt")
write.table(data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)

# Print a message indicating completion
cat("Substitutions completed. Output written to :", output_file, "\n")

# Leggi il file appena scritto
read_data <- read.delim(output_file, col.names = c('Network.Object..FROM.', 'Effect', 'Network.Object..TO.'))

# Controllo della corretta separazione dei campi nel file appena letto
# if (all(colnames(read_data) == colnames(data)) &&
#     all(read_data$`Network.Object..FROM.` == data$`Network.Object..FROM.`) &&
#     all(read_data$Effect == data$Effect) &&
#     all(read_data$`Network.Object..TO.` == data$`Network.Object..TO.`)) {
#   cat("File read successfully with correct field separation.\n")
# } else {
#   cat("Error: Field separation in the read file does not match.\n")
# }

