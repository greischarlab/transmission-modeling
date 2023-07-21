delays <- read.csv("Average Delays.csv")

separatedDelay <- delays[order(delays$type), ]
separatedDelay$color <- "slateblue"
separatedDelay$color[separatedDelay$type=="respiratory"] <- "violet"

png("SympSheddingDelay.png", width=1000, height=600)
barplot(separatedDelay$delay.symp.shed, names=separatedDelay$Virus,
        main="Delay Between Symptom Onset and Shedding Onset (Days)", 
        col=separatedDelay$color, cex.names=1.3, cex.main=1.3, cex=1.3,
        ylim=c(-0.3, 1))
abline(h=0, lwd=3, lty=2)
legend(4, 0.5, legend=c("Digestive", "Respiratory"), 
       lty=1, lwd=2, col=c("slateblue", "violet"), cex=1.3)
dev.off()