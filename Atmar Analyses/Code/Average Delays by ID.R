setwd("C:/Users/ktopa/OneDrive/Documents/Research Project/Summer URE")
delaysbyid <- read.csv("Delays by ID Norovirus.csv")
vl.atsymponset <- read.csv("Viral Load at Symptom Onset.csv")
patientparams <- read.csv("Atmar Patient Parameters local minimum.csv")

# Exclude participants with no clear peak in vl trajectories
delaysbyid <- delaysbyid[delaysbyid$id != 731 &
                          delaysbyid$id != 723 &
                          delaysbyid$id != 717,]

patientparams <- patientparams[patientparams$patient != 731 &
                                 patientparams$patient != 723 &
                                 patientparams$patient != 717,]

# Get average delays per id
avg.delays <- data.frame(dose = rev(unique(delaysbyid$dose)), 
                         average.delay = 
                           rep(NA, length(unique(delaysbyid$dose))))


for (dose in avg.delays$dose) {
  avg.delays[avg.delays$dose==dose,]$average.delay <- 
    mean(delaysbyid[delaysbyid$dose==dose,]$delay)
}

# Showing that p2 is correlated with delay

p2.linregdata <- data.frame(p2=patientparams[order(patientparams$patient),]$p2,
                            delay=delaysbyid[order(delaysbyid$id),]$delay)

lm.fit <- lm(delay ~ p2, data=p2.linregdata)

pval <- summary(lm.fit)$coefficients[2,4]

# Correlation coefficient is fairly close to -1
cor.coefficient <- cor(p2.linregdata$delay, p2.linregdata$p2)

# I will label the axis with the appropriate infectious doses (doing this to
# space out the IDs evenly)
id <- delaysbyid$dose
id[id==4.8] <- 1
id[id==48] <- 2
id[id==4800] <- 3

delays <- delaysbyid$delay
mean.delays <- avg.delays$average.delay
mean.delays.id <- avg.delays$infectious.dose

colors <- c("#D60270","#9B4F96","#0038A8")
cex <- 8
panellabelratio <- .98333

{
  postscript(file="Atmar Presymptomatic Transmission.eps", paper="special", 
             height=60, width=24, horizontal=FALSE)
  par(mar=c(3, 2, 4, 2), oma=c(20, 35, 2, 2),
      mgp=c(0, 7, 0), mai=c(2, 0, 1.5, 0))
  plot(jitter(id), delays,
          col=colors[factor(id)], xlim=c(.5, 3.5),
          ylim=c(-5, 140), cex.main=cex, cex.lab=cex, cex.axis=cex, cex=cex + 2, 
          font.main=2, font.lab=2, pch=ifelse(delays < 0, 16, 1),
          xaxt="n",
          ylab = "",
          xlab = "")
  axis(1, at=seq(1,3), labels=c("4.8", "48", "4800"), lwd=0, line=0, cex.axis=cex)
  mtext(text="Delay Between Symptom Onset and \n Time of Peak Pathogen Load (Hours)",
        side=2, line=22, cex=cex * (3/4))
  mtext(text="Infectious Dose", side=1, line=18, cex=cex * (3/4))
  
  #Average delays per infectious dose
  segments(0.7, mean.delays[1], 1.3, col=colors[1], lwd=cex + 2)
  segments(1.7, mean.delays[2], 2.3, col=colors[2], lwd=cex + 2)
  segments(2.7, mean.delays[3], 3.3, col=colors[3], lwd=cex + 2)
  
  #Pre-symptomatic vs post-symptomatic line
  abline(h=0, lty=3, lwd=cex + 2)
  text(panellabelratio * 3 + .5, panellabelratio * 145 - 5, labels="A", cex=cex, font=2)
  
  # Viral load vs. symptom onset
  
  plot(vl ~ jitter(symp.onset), data=vl.atsymponset, col=colors[factor(dose)],
       cex=cex+2, cex.lab=cex, cex.axis=cex,
       pch=ifelse(vl.atsymponset$presymp, 16, 1), ylim=c(7, 11),
       ylab="", xlab="")
  
  text(panellabelratio * 2 + 1, panellabelratio * 4 + 7, labels="B", cex=cex, font=2)
  mtext(text="Viral Load At Symptom Onset",
        side=2, line=21, cex=cex * (3/4))
  mtext(text="Day of Symptom Onset", side=1, line=18, cex=cex * (3/4))
  
  # Growth Rate vs. Delay
  
  plot(patientparams[order(patientparams$patient),]$p2,
       delaysbyid[order(delaysbyid$id),]$delay, 
       col=colors[factor(delaysbyid[order(delaysbyid$id),]$dose)],
       pch=ifelse(delaysbyid[order(delaysbyid$id),]$delay < 0, 16, 1),
       cex=cex+2, cex.lab=cex+2, cex.axis=cex, 
       ylab="", xlab="", main="Observed Data From Norovirus Study", yaxt="n",
       cex.main=cex-2)
  
  axis(side=2, at=seq(-10, 140, by=10), cex.axis=cex, labels=c("-10",
                                                 "0",
                                                 "10",
                                                 "20",
                                                 "30",
                                                 "40",
                                                 "50",
                                                 "60",
                                                 "70",
                                                 "80",
                                                 "90",
                                                 "100",
                                                 "110",
                                                 "120",
                                                 "130",
                                                 "140"))
  
  abline(lm.fit, col="red")
  
  abline(h=0, lty=3, lwd=cex + 2)
  
  text(panellabelratio * 145 - 5, panellabelratio * 5 - 4, labels="C", cex=cex, font=2)
  mtext(text="Viral Growth Rate",
        side=1, line=20, cex=cex * (3/4))
  mtext(text="Delay Between Symptom Onset \n And Peak Viral Load", side=2, 
        line=22, cex=cex * (3/4))
  
  legend(-3, 110, legend=c("Pre-symptomatic", "Post-symptomatic", 
                            "Least Squares Regression Line"), pch=c(19, 1, NA),
                            lty=c(NA, NA, 1), col=c("black", "black", "red"),
         cex=cex * (2/3))
  dev.off()
}
