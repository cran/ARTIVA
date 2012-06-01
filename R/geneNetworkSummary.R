geneNetworkSummary <- function(ARTIVAres, Threshold){

 # Only interactions with a posterior probability higher than the specified threshold are written
 ARTIVAsubRes = ARTIVAres[ARTIVAres$PostProb >= Threshold,]

 # Interactions are ordered 
 ARTIVAsubRes = ARTIVAsubRes[order(ARTIVAsubRes$PostProb, decreasing = T),]
    
 # To get the signs of the interactions
 InteractionSigns = rep("NA", nrow(ARTIVAsubRes))
 InteractionSigns[ARTIVAsubRes$CoeffMean >= 0] = "+"
 InteractionSigns[ARTIVAsubRes$CoeffMean < 0]  = "-"

 # Information to be written
 ResTable =  data.frame(cbind(ARTIVAsubRes$Parent, ARTIVAsubRes$Target, ARTIVAsubRes$PostProb,
                                 ARTIVAsubRes$CPstart, ARTIVAsubRes$CPend, InteractionSigns))  
 colnames(ResTable) <- c("Parent", "Target", "PostProb", "CPstart", 
                            "CPend", "InteractionSign")
 #Calculate the number of pages (30 lines per page)
 nbpages <- nrow(ResTable)/30
 
 i <- 1
 # Browse the dataframe
 while(i <= nrow(ResTable))
   {
     if((i+29)/30 <= nbpages)   # It does not exceed the number of pages, so write 30 lines more
       textplot(ResTable[i:(i+29),], cex = 0.6, cmar = 1.1, rmar = 1, show.rownames = F, show.colnames = T,
                halign = "left",
                valign = "top")
     else   # It exceeds the number of pages, so write the end of the dataframe
       textplot(ResTable[i:nrow(ResTable),], cex = 0.6, cmar = 1.1, rmar = 1, show.rownames = F, show.colnames = T,
                halign = "left",
                valign = "top")
     i <- i+30
   }
 
# End of function
}

