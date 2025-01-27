library(visNetwork)
library(dplyr)
plot_grns <- function(DF)
{
  ids <-  unique(c(DF$From, DF$To))
  
  source_Nodesize <- data.frame(plyr::count(DF$From))
  source_Nodesize$freq <- source_Nodesize$freq + 1
  
  to_Nodesize <- data.frame(x = unique(DF$To[!(DF$To %in% DF$From)]))
  to_Nodesize$freq <- 1
  nodes_size <- rbind(source_Nodesize, to_Nodesize)
  colnames(nodes_size)[2] <- "value"
  
  nodes <- data.frame(id = ids, label = ids,
                      shape = ifelse(ids %in% DF$From & ids %in% DF$To, "circle", "bigdot"), #change the shapes to circle and dot to get the v1
                      shadow = TRUE,
                      color = ifelse(ids %in% DF$From & !(ids %in% DF$To), "darkblue", "lightblue"),
                      font.color = "white")
  
  nodes <- merge(nodes, nodes_size, by.x = "id", by.y = "x")
  edges <- data.frame(from = DF$From, to = DF$To, arrows = "to",
                      color = case_when(
                        DF$effect == "-->" ~ "darkgreen",
                        DF$effect == "--?" ~ "darkgrey",
                        DF$effect == "--|" ~ "darkred",
                      ))
  
  net <- visNetwork(nodes, edges, width = "100%", height = 900)%>% 
    visNodes(shadow = list(enabled = TRUE)) %>% 
    visEdges(smooth = list(enabled = TRUE)) %>% 
    visInteraction(navigationButtons = TRUE)  %>% 
    visOptions(highlightNearest = list(enabled = TRUE, algorithm = "hierarchical"), nodesIdSelection = TRUE)
  
  net
}
interactions <- read.csv("4classes/consensus_analysis/NetworkPhenotype2.txt", header = F, sep = "\t")
colnames(interactions) <- c("From", "effect", "To")
network <- plot_grns(interactions)
