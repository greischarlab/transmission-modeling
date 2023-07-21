phylo.dist <- read.csv("Phylogenetic Distances From Norwalk.csv")
delays <- read.csv("Average Delays.csv")

orderedDist <- phylo.dist[order(phylo.dist$Virus), ]
orderedDelay <- delays[order(delays$Virus), ]

norwalk <- grepl("Norwalk Virus", orderedDelay$Virus)
orderedDelay=orderedDelay[!norwalk, ]

png("H1.png", width=900, height=800)
plot(orderedDist$Phylo.dist, orderedDelay$delay.symp.peak, 
     xlab="Phylogenetic Distance from Norwalk virus", 
     ylab="Delay Between Symptom Onset and Peak Viral Load",
     col="purple", xlim=c(0.2, 0.8), cex=1.5, cex.lab=1.3)
abline(lm(orderedDelay$delay.symp.peak ~ orderedDist$Phylo.dist), 
       col="purple", lty=2)
text(orderedDist$Phylo.dist + .12, orderedDelay$delay.symp.peak, 
     labels=orderedDelay$Virus, cex=1.3)
dev.off()



