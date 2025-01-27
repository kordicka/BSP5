#
# Instructions to:
#
# 1) map differential expression statistics from limma to
# gene regulatory networks in GeneGo MetaCore
#
# 2) run R code to process the nodemap- and interaction-data
# downloaded from GeneGo MetaCores to create the input files
# for network perturbation analysis
#
# 3) run Jave code for the network perturbation analysis
#


#
# 1)
#
# To prepare your input data for the processing below, you need to run the following steps:
#
# 1.) Conduct a limma case-control differential expression analysis and write the output of the topTable-function to an Excel-file using the openxlsx-package and the write.xlsx-function (output: case_vs_controls_ranking.xlsx)
# 2.) Open case_vs_controls_ranking.xlsx in Excel or another spreadsheet program and copy all gene symbols with FDR <0.05 (or if there are only few, use the nominal p-value cut-off: p < 0.05)
# 3.) Login in GeneGo MetaCore (https://portal.genego.com/cgi/data_manager.cgi), goto "Build Network" - Paste the copied gene list into the text box - select "Gene symbol (official only)", deactivate: Group of proteins, complex of proteins, group of complexes; then submit
# 4.) Network algorithm: Select "Direct interactions"
# 5.) Prefilters: Select "Species: human"
# 6.) Object types: first select all, then remove DNA, RNA, Compound, Inorganic ion, Predicted metabolite
# 7.) Interaction types - select: B, TR, IE, cRT, Rg (Transcription regulation, Influence on expression, co-regulation of transcription (bottom: effects: only activating/inhibiting interactions)
# 8.) File - Export: nodemap (gene symbol/object mappings, select: homo sapiens for the nodemap: twice --> also for "Through") - save as: "case_vs_controls_network_nodes.xlsx"
# 9.) File - Export: interactions, Excel format - save as: "case_vs_controls_network_interactions.xlsx"
# 10.) For "case_vs_controls_network_nodes.xlsx": Open in Excel, need to remove hyperlinks first: mark all cells in columns C and D: copy them, then paste them again using "Special" - and "paste as values", then save the file under the same name
# 11.) Continue in R with the code below
#


#
# 2)
#
# Set your current working directory containing the downloaded nodemap and interaction data from GeneGo MetaCore
# workingdir = '/Users/downloads/'  # MAC example
# workingdir = 'C:/Users/downloads/' # #Windows example
# setwd(workingdir) 
#

# Required library for reading Excel files
require('xlsx')

# Reading network object and network interaction data from Excel files
netobj = read.xlsx("case_vs_controls_network_nodes.xlsx", sheetIndex = 1, startRow = 3)
net = read.xlsx("case_vs_controls_network_interactions.xlsx", sheetIndex = 1, startRow = 3)



# Function to preprocess gene list data from GeneGo nodemap file
preprocess_genelist = function(netobj, net, nodefile="nodemap.txt", intfile="interactions.txt", genefile="genelist.txt") {
    # Creating nodemap
    nodemap = cbind(as.matrix(netobj$Network.Object.Name), as.matrix(netobj$Gene.Symbol))
    nodemap = unique(nodemap)
    write.table(nodemap, file=nodefile, sep="\t", quote=F, col.names=F, row.names=F)
    
    # Processing interactions
    network = cbind(as.matrix(nodemap[match(net[,2],nodemap[,1]),2]), as.matrix(nodemap[match(net[,4],nodemap[,1]),2]), as.matrix(net[,6]))
    cat("",file=intfile,append=F)
    for(i in 1:nrow(network)){
        write(paste(network[i,1],"\t",ifelse(network[i,3]=="Unspecified","--?",ifelse(network[i,3]=="Activation","-->","--|")),"\t",network[i,2], sep=""),file=intfile,append=T)
    }
    
    # Generating gene list
    write.table(unique(as.matrix(netobj$Gene.Symbol)), file=genefile, quote=F, col.names=F, row.names=F)
    
    return()
}


# Function for processing expression data
expprocess = function(datmat, outfile="expression.txt") {
    genesymbol = as.matrix(datmat[,1])
    logfcs = datmat[,2]
    write.table(cbind(genesymbol, ifelse(logfcs < 0,0,1)), sep="\t", quote=F, col.names=F, row.names=F, file=outfile)
    return()
}

# Executing preprocess_genelist function with specified file paths
preprocess_genelist(netobj, net, paste(workingdir,"case_vs_controls_nodemap.txt",sep=""), paste(workingdir,"case_vs_controls_interactions.txt",sep=""), paste(workingdir,"case_vs_controls_genelist.txt",sep=""))

# Reading gene expression data (this is the output from a limma case-control analysis using the topTable-function and written to an Excel-file in xlsx-format)
expdat_case_vs_controls = read.xlsx('case_vs_controls_ranking.xlsx',sheetIndex = 1, startRow = 1)

# Prepare the expression data for the network perturbation analysis
expprocess(cbind(expdat_case_vs_controls[,c(1,3)]), paste(workingdir,"case_vs_controls_expression.txt",sep=""))


#
# 3)
#
# Network perturbation analysis (using Java and the jar-files Preprocessor.jar, DifferentialNetworkAnalysis.jar, CommonNetworkGenerator.jar, ComputeCycles.jar, SteadyStateCalculator.jar, PerturbagenListGenerator.jar and BruteForcePerturbationsUpdated.jar
#
#
#
## Pre-process input data
# java -jar Preprocessor.jar . case_vs_controls_nodemap.txt case_vs_controls_interactions.txt case_vs_controls_genelist.txt
## rename adjacency.txt to adjacency_case_vs_controls.txt
# mv adjacency.txt adjacency_case_vs_controls.txt
#
#
## run differential network analysis
# java -jar DifferentialNetworkAnalysis.jar case_vs_controls_expression.txt adjacency_case_vs_controls.txt GAResult.txt 0 true 1000 50 .
#
#
## rename NetworkPhenotype1.txt to sars_cov2_variant2_NetworkPhenotype1.txt
## rename NetworkPhenotype2.txt to sars_cov2_variant2_NetworkPhenotype2.txt
# mv NetworkPhenotype1.txt case_vs_controls_NetworkPhenotype1.txt
# mv NetworkPhenotype2.txt case_vs_controls_NetworkPhenotype2.txt
#
##
## create Cytoscape visualizations out of these networks
## import into Cytoscape:
## 1) File - import - import _Network_ from file: case_vs_controls_NetworkPhenotype1.txt
## 2) Advanced options: uncheck "use first row as column names", delimiter: Tabs
## 3) Column1 = source node; Column3 = target node, Column2 = interaction type
## 4) Create hierarchical network layout using "yfiles Hierarchical Layout"
## 5) Style - Edge: Target Arrow Shape - Column "interaction" - Discrete Mapping : ->, -|
## 6) Map gene expression data onto Node colors (import: import from table, export to tab-delimited text-file from case_vs_control_ranking.xlsx, node - fill color = log2foldchange, mapping = Continuous)
## 7) Save as: case_vs_controls_NetworkPhenotype1.cys (respectively, case_vs_controls_NetworkPhenotype2.cys, for the 2nd phenotype)
##
#
#
## Perturbation analysis
# java -jar CommonNetworkGenerator.jar case_vs_controls_NetworkPhenotype1.txt case_vs_controls_NetworkPhenotype2.txt case_vs_controls_CommonNetworkGenerator_Output.txt
# java -jar DifferentialNetworkGenerator.jar case_vs_controls_NetworkPhenotype1.txt case_vs_controls_NetworkPhenotype2.txt case_vs_controls_DifferentialNetworkGenerator_Output.txt
# java -jar ComputeCycles.jar case_vs_controls_CommonNetworkGenerator_Output.txt case_vs_controls_expression.txt case_vs_controls_pos.txt case_vs_controls_neg.txt
#
# cat case_vs_controls_pos.txt
# cat case_vs_controls_neg.txt
## check that there are positive or negative cycles (if not, the perturbation analysis may not be feasible, e.g. if the network is too small --> consider to relax p-value cut-off)
#
#
# java -jar SteadyStateCalculator.jar case_vs_controls_expression.txt case_vs_controls_NetworkPhenotype1.txt 1 case_vs_controls_SteadyStateCalculatorN1.txt
# java -jar SteadyStateCalculator.jar case_vs_controls_expression.txt case_vs_controls_NetworkPhenotype2.txt 2 case_vs_controls_SteadyStateCalculatorN2.txt
#
# java -jar PerturbagenListGenerator.jar case_vs_controls_pos.txt case_vs_controls_neg.txt case_vs_controls_DifferentialNetworkGenerator_Output.txt case_vs_controls_SteadyStateCalculatorN1.txt case_vs_controls_PerturbagenListGeneratorN1.txt
# java -jar PerturbagenListGenerator.jar case_vs_controls_pos.txt case_vs_controls_neg.txt case_vs_controls_DifferentialNetworkGenerator_Output.txt case_vs_controls_SteadyStateCalculatorN2.txt case_vs_controls_PerturbagenListGeneratorN2.txt
#
# java -jar BruteForcePerturbationsUpdated.jar case_vs_controls_expression.txt case_vs_controls_NetworkPhenotype1.txt 1 case_vs_controls_PerturbagenListGeneratorN1.txt 1 500000 case_vs_controls_BruteForcePerturbationsUpdatedN1_1.txt
# java -jar BruteForcePerturbationsUpdated.jar case_vs_controls_expression.txt case_vs_controls_NetworkPhenotype1.txt 1 case_vs_controls_PerturbagenListGeneratorN1.txt 2 500000 case_vs_controls_BruteForcePerturbationsUpdatedN1_2.txt
# java -jar BruteForcePerturbationsUpdated.jar case_vs_controls_expression.txt case_vs_controls_NetworkPhenotype1.txt 1 case_vs_controls_PerturbagenListGeneratorN1.txt 3 500000 case_vs_controls_BruteForcePerturbationsUpdatedN1_3.txt
#
# sort -rn case_vs_controls_BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
# mv Z.txt case_vs_controls_BruteForcePerturbationsUpdatedN1_1.txt
#
#
# sort -rn case_vs_controls_BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
# mv Z.txt case_vs_controls_BruteForcePerturbationsUpdatedN1_2.txt
#
#
# sort -rn case_vs_controls_BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
# mv Z.txt case_vs_controls_BruteForcePerturbationsUpdatedN1_3.txt
#