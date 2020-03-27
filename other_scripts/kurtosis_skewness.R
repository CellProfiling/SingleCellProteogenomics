

#library(moments)
#skewness(CCNB1$Mean_ab_Cyto)
#skewness(GAPDH$Mean_ab_Cyto)
#skewness(DUSP18$Intensity_MeanIntensity_ResizedAb)
#library(agrmt)

#ajus(GAPDH$Mean_ab_Cyto, tolerance=0.1)
#ajus(CCNB1$Mean_ab_Cyto, tolerance=0.7, variant="modified")
#ajus(DUSP18$Intensity_MeanIntensity_ResizedAb, tolerance=0.7, variant="modified")
#library(diptest)
#dip.test(DUSP18$Intensity_MeanIntensity_ResizedAb)
#dip.test(GAPDH$Mean_ab_Cyto)
#dip.test(CCNB1$Mean_ab_Cyto)

#kurtosis(GAPDH$Mean_ab_Cyto)
#kurtosis(CCNB1$Mean_ab_Cyto)
#kurtosis(DUSP18$Intensity_MeanIntensity_ResizedAb)

#hist(CCNB1$Mean_ab_Cyto, prob=TRUE)
#lines(density(CCNB1$Mean_ab_Cyto))
#hist(GAPDH$Mean_ab_Cyto, prob=TRUE)
#lines(density(GAPDH$Mean_ab_Cyto), lwd=2)
#hist(DUSP18$Intensity_MeanIntensity_ResizedAb,prob=TRUE)
#lines(density(DUSP18$Intensity_MeanIntensity_ResizedAb, adjust=2), lwd=2)
#pdf("gp.pdf")
#gp<-ggplot2.violinplot(data=CCNB1$Mean_ab_Cyto, scale=0.1, addDot=TRUE,dotSize=1.7,removePanelGrid=TRUE,removePanelBorder=TRUE,axisLine=c(0.5, "solid", "black"), backgroundColor="white", xtitle="Dose (mg)", ytitle="length", mainTitle="CCNB1",  dotPosition="jitter", jitter=0.2)
#print(gp)
#dev.off()

workingdir<-setwd("/Users/diana.telessemian/Desktop/gaussian_per_plate/norm_volin_plot/")
well_plate<- unique(data$well_plate)
for (j in well_plate) {
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  datanuccyto<- data[grep(j,data$well_plate),]
  
  library(vioplot)
  require(viopoints)
  x <- rnorm(10)
  y <- rnorm(10)
  pdf(paste(workingdir,paste( "vio_plot_cyto",j,".pdf", sep=""),sep="/"),height=3, width=8 )
  #pdf("vio_plot_cyto.pdf", height=5, width=8) 
  par(mar=c(3,7,2,5))
  par(mfrow=c(1,2))
  #violin_plot<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Mean_ab_Cyto)), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Mean Intensity", las=1)
  #viopoints(datanuccyto$Mean_ab_Cyto, main="", col="aquamarine3", add=TRUE)
  #vioplot(datanuccyto$Mean_ab_Cyto ,col="transparent", add=TRUE)
  #hist<- hist(datanuccyto$Mean_ab_Cyto,prob=TRUE, main="", las=1)
  #lines(density(datanuccyto$Mean_ab_Cyto, adjust=2), lwd=2, main="")
  violin_plot_norm<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto))), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Normalized Mean Intensity",xlab="", las=1)
  viopoints(datanuccyto$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto), main="", col="aquamarine3", add=TRUE)
  vioplot(datanuccyto$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto) ,col="transparent", add=TRUE)
  hist_norm<-hist(datanuccyto$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto),prob=TRUE, main="", las=1, ylab="Normalized Density")
  lines(density(datanuccyto$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto), adjust=2), lwd=2, main="")
  dev.off()
  
  
  pdf(paste(workingdir,paste( "vio_plot_nuc",j,".pdf", sep=""),sep="/"),height=3, width=8 )
  #pdf("vio_plot_nuc.pdf", height=3, width=8) 
  par(mar=c(3,7,2,5))
  par(mfrow=c(1,2))
  #violin_plot_nuc<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Intensity_MeanIntensity_ResizedAb)), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Mean Intensity", las=1)
  #viopoints(datanuccyto$Intensity_MeanIntensity_ResizedAb, main="", col="aquamarine3", add=TRUE)
  #vioplot(datanuccyto$Intensity_MeanIntensity_ResizedAb ,col="transparent", add=TRUE)
  #hist_nuc<- hist(datanuccyto$Intensity_MeanIntensity_ResizedAb,prob=TRUE, main="", las=1)
  #  lines(density(datanuccyto$Intensity_MeanIntensity_ResizedAb, adjust=2), lwd=2, main="")
  violin_plot_nuc_norm<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb))), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Normalized Mean Intensity",xlab="", las=1)
  viopoints(datanuccyto$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb), main="", col="aquamarine3", add=TRUE)
  vioplot(datanuccyto$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb) ,col="transparent", add=TRUE)
  hist_nuc_norm<-hist(datanuccyto$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb),prob=TRUE, main="", las=1, ylab="Normalized Density")
  lines(density(datanuccyto$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb), adjust=2), lwd=2, main="")
  dev.off()
  
  pdf(paste(workingdir,paste( "vio_plot_cell",j,".pdf", sep=""),sep="/"),height=3, width=8 )
  #pdf("vio_plot_cell.pdf", height=5, width=8) 
  par(mar=c(3,7,2,5))
  par(mfrow=c(1,2))
  #violin_plot_cell<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Mean_ab_cell)), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Mean Intensity", las=1)
  # viopoints(datanuccyto$Mean_ab_cell, main="", col="aquamarine3", add=TRUE)
  #vioplot(datanuccyto$Mean_ab_cell ,col="transparent", add=TRUE)
  #hist_cell<- hist(datanuccyto$Mean_ab_cell,prob=TRUE, main="", las=1)
  # lines(density(datanuccyto$Mean_ab_cell, adjust=2), lwd=2, main="")
  violin_plot_cell_norm<-plot(x,y, xlim=c(0,2), ylim=c(0,max(datanuccyto$Mean_ab_cell/max(datanuccyto$Mean_ab_cell))), col="transparent", frame.plot=FALSE,xaxt="n.axis", ylab="Normalized Mean Intensity",xlab="", las=1)
  viopoints(datanuccyto$Mean_ab_cell/max(datanuccyto$Mean_ab_cell), main="", col="aquamarine3", add=TRUE)
  vioplot(datanuccyto$Mean_ab_cell/max(datanuccyto$Mean_ab_cell) ,col="transparent", add=TRUE)
  hist_cell_norm<-hist(datanuccyto$Mean_ab_cell/max(datanuccyto$Mean_ab_cell),prob=TRUE, main="", las=1, ylab="Normalized Density")
  lines(density(datanuccyto$Mean_ab_cell/max(datanuccyto$Mean_ab_cell), adjust=2), lwd=2, main="")
  dev.off()
}



workingdir<-setwd("/Users/diana.telessemian/Desktop/gaussian_per_plate/kurtosis_skewness/")
well_plate<- unique(data$well_plate)
for (j in well_plate) {
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  datanuccyto<- data[grep(j,data$well_plate),]
  
  
  stdev_nuc<-sd(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  stdev_cyto<-sd(datanuccyto$Mean_ab_Cyto)
  stdev_cell<-sd(datanuccyto$Mean_ab_cell)
  median_nuc<-median(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  median_cyto<-median(datanuccyto$Mean_ab_Cyto)
  median_cell<-median(datanuccyto$Mean_ab_cell)
  mean_nuc<-mean(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  mean_cyto<-mean(datanuccyto$Mean_ab_Cyto)
  mean_cell<-mean(datanuccyto$Mean_ab_cell)
  capture.output(median_nuc,file=paste(workingdir,paste( "median_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(median_cyto,file=paste(workingdir,paste( "median_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(median_cell,file=paste(workingdir,paste( "median_cell_",j,".csv", sep=""),sep="/"))
  capture.output(stdev_nuc,file=paste(workingdir,paste( "stdev_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(stdev_cyto,file=paste(workingdir,paste( "stdev_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(stdev_cell,file=paste(workingdir,paste( "stdev_cell_",j,".csv", sep=""),sep="/"))
  capture.output(mean_nuc,file=paste(workingdir,paste( "mean_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(mean_cyto,file=paste(workingdir,paste( "mean_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(mean_cell,file=paste(workingdir,paste( "mean_cell_",j,".csv", sep=""),sep="/"))
}
workingdir<-setwd("/Users/diana.telessemian/Desktop/gaussian_per_plate/kurtosis_skewness/")
well_plate<- unique(data$well_plate)
for (j in well_plate) {
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  datanuccyto<- data[grep(j,data$well_plate),]
  require(diptest)
  require(moments)
  Kurtosis_nuc<-kurtosis(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  Kurtosis_cyto<-kurtosis(datanuccyto$Mean_ab_Cyto)
  Kurtosis_cell<-kurtosis(datanuccyto$Mean_ab_cell)
  skewness_nuc<-skewness(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  skewness_cyto<-skewness(datanuccyto$Mean_ab_Cyto)
  skewness_cell<-skewness(datanuccyto$Mean_ab_cell)
  dip_nuc<-dip.test(datanuccyto$Intensity_MeanIntensity_ResizedAb)$p.value
  dip_cyto<-dip.test(datanuccyto$Mean_ab_Cyto)$p.value
  dip_cell<-dip.test(datanuccyto$Mean_ab_cell)$p.value
  capture.output(Kurtosis_nuc,file=paste(workingdir,paste( "Kurtosis_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(Kurtosis_cyto,file=paste(workingdir,paste( "Kurtosis_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(Kurtosis_cell,file=paste(workingdir,paste( "Kurtosis_cell_",j,".csv", sep=""),sep="/"))
  capture.output(skewness_nuc,file=paste(workingdir,paste( "skewness_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(skewness_cyto,file=paste(workingdir,paste( "skewness_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(skewness_cell,file=paste(workingdir,paste( "skewness_cell_",j,".csv", sep=""),sep="/"))
  capture.output(dip_nuc,file=paste(workingdir,paste( "dip_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(dip_cyto,file=paste(workingdir,paste( "dip_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(dip_cell,file=paste(workingdir,paste( "dip_cell_",j,".csv", sep=""),sep="/"))
  
  
  ######Normality test
  require(nortest)
  normality_nuc<-ad.test(datanuccyto$Intensity_MeanIntensity_ResizedAb)$p.value
  normality_cyto<-ad.test(datanuccyto$Mean_ab_Cyto)$p.value
  normality_cell<-ad.test(datanuccyto$Mean_ab_cell)$p.value
  capture.output(normality_nuc,file=paste(workingdir,paste( "normality_nuc",j,".csv", sep=""),sep="/"))
  capture.output(normality_cyto,file=paste(workingdir,paste( "normality_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(normality_cell,file=paste(workingdir,paste( "normality_cell_",j,".csv", sep=""),sep="/"))
  
  
}

workingdir<-setwd("/Users/dianatelessemian/Desktop/kruskal_fv/clustering")
well_plate<- unique(data$well_plate)
for (j in well_plate) {
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  datanuccyto<- data[grep(j,data$well_plate),]
  
  
  
  nbofcells_expressing_nuc<-nrow(datanuccyto[datanuccyto$Intensity_MeanIntensity_ResizedAb > mean(datanuccyto$Intensity_MeanIntensity_ResizedAb),])
  nbofcells_expressing_cyto<-nrow(datanuccyto[datanuccyto$Mean_ab_Cyto > mean(datanuccyto$Mean_ab_Cyto),])
  nbofcells_expressing_cell<-nrow(datanuccyto[datanuccyto$Mean_ab_cell > mean(datanuccyto$Mean_ab_cell),])
  nbofcells<-length(datanuccyto$Intensity_MeanIntensity_ResizedAb)
  capture.output(nbofcells_expressing_nuc,file=paste(workingdir,paste( "nbofcells_expressing_nuc_",j,".csv", sep=""),sep="/"))
  capture.output(nbofcells_expressing_cyto,file=paste(workingdir,paste( "nbofcells_expressing_cyto_",j,".csv", sep=""),sep="/"))
  capture.output(nbofcells_expressing_cell,file=paste(workingdir,paste( "nbofcells_expressing_cell_",j,".csv", sep=""),sep="/"))
  capture.output(nbofcells,file=paste(workingdir,paste( "nbofcells_",j,".csv", sep=""),sep="/"))
}


#######clustering

library(xlsx)
ccd_notccd_features<-read.xlsx("/Users/diana.telessemian/Desktop/gaussian_per_plate/features_all.xlsx",header=T ,sheetIndex =1 , stringsAsFactors = FALSE)
###mean_phases_var_norm<-read.xlsx("/Users/dianatelessemian/Desktop/Kruskal_fv/clustering/features/mean_G1_G1S_SG2nuc_cyto_cell.xlsx",header=T ,sheetIndex =1 , stringsAsFactors = FALSE)
##ccd_notccd_features_new<-merge(ccd_notccd_features,mean_phases_var_norm, by="well_plate" )
ccd_notccd_features_compartments<-read.xlsx("/Users/diana.telessemian/Desktop/gaussian_per_plate/features_all.xlsx",header=T ,sheetIndex =3 , stringsAsFactors = FALSE)

ccd_features<-subset(ccd_notccd_features,ccd_notccd_features$correlation=="ccd")
ccd_features_compartments<-subset(ccd_notccd_features_compartments,ccd_notccd_features_compartments$corr_FDR_explainedvar=="ccd")
notccd_features<-subset(ccd_notccd_features,ccd_notccd_features$correlation=="notccd")
notccd_features_compartments<-subset(ccd_notccd_features_compartments,ccd_notccd_features_compartments$corr_FDR_explainedvar=="notccd")

ccd_mean_stdev<-data.frame(ccd_features[,c(75,81)])
#mean_stdev<-data.frame(ccd_features[c(1:397),c(94,100)])
mean<-as.numeric(ccd_features$mean_cell.y)
stdev<-as.numeric(ccd_features$X.stdev_cell)
results<-kmeans(ccd_mean_stdev,3)
results
results$size
table(ccd_features$NA..1,results$cluster)
plot(mean,stdev, col=results$cluster)
plot(log(mean),log(stdev), col=results$cluster)
plot(mean,stdev, col=factor(ccd_features$NA..1))


ccd_compartments_mean_stdev<-data.frame(ccd_features_compartments[,c(49,51)])
#mean_stdev<-data.frame(ccd_features[c(1:397),c(94,100)])
mean<-as.numeric(ccd_features_compartments$mean_compartment.y)
stdev<-as.numeric(ccd_features_compartments$X.stdev_compartment)
results<-kmeans(ccd_mean_stdev,3)
results
results$size
table(ccd_features$NA..1,results$cluster)
plot(mean,stdev, col=results$cluster)
plot(log(mean),log(stdev), col=results$cluster)
plot(mean,stdev, col=factor(ccd_features$NA..1))



p_ccd<-subset(kurtosis_skewness, kurtosis_skewness$Gene=="CCNB1"| kurtosis_skewness$Gene=="GATA6"| kurtosis_skewness$Gene=="NECAP2"| kurtosis_skewness$Gene=="AP3D1")
p_notccd<-subset(kurtosis_skewness, kurtosis_skewness$Gene=="SC5D"| kurtosis_skewness$Gene=="GAPDH"| kurtosis_skewness$Gene=="TRIM32"| kurtosis_skewness$Gene=="INCA1")

mycolors = c('gray','darkblue')
kurtosis_skewness$col<-mycolors[kurtosis_skewness$corr_FDR_explainedvar]
kurtosis_skewness<-data.frame(ccd_notccd_features_compartments[,c(2,4,5,38,40,41)])
kurtosis_skewness$X.Kurtosis_compartment<-as.numeric(kurtosis_skewness$Kurtosis_comp)
kurtosis_skewness$skewness_compartment<-as.numeric(kurtosis_skewness$skewness_comp)
#kurtosis_skewness$X.stdev_compartment<-as.numeric(kurtosis_skewness$X.stdev_compartment)
kurtosis_skewness$correlation<-factor(kurtosis_skewness$corr_FDR_explainedvar)
plot(kurtosis_skewness$skewness_compartment,kurtosis_skewness$X.Kurtosis_compartment,  col=mycolors[kurtosis_skewness$correlation],frame.plot=FALSE, ylab="Kurtosis",xlab="skewness", las=1 )
points(p_ccd$skewness_compartment, p_ccd$X.Kurtosis_compartment, pch=16, cex=2, col="gray28")
points(p_notccd$skewness_compartment, p_notccd$X.Kurtosis_compartment, pch="*", cex=3, col="darkblue")
text(0.13,16,"GATA6", col="gray28", cex=0.8)
text(0.95,24,"CCNB1", col="gray28", cex=0.8)
text(3.9345,42,"NECAP2", col="gray28", cex=0.8)
text(12,202,"AP3D1", col="gray28", cex=0.8)
text(0.13,31,"SC5D", col="darkblue", cex=0.8)
text(1.910,38,"GAPDH", col="darkblue", cex=0.8)
text(3.9345,50,"TRIM32", col="darkblue", cex=0.8)
text(10,145.005,"INCA1", col="darkblue", cex=0.8)
legend("bottomright",legend=levels(kurtosis_skewness$correlation), fill =mycolors, cex=0.8)
#####kurtosis skewness stdev
col=mycolors[kurtosis_skewness$correlation]
require(scatterplot3d)
scatterplot3d(x = kurtosis_skewness$skewness_compartment, # x axis
              y = kurtosis_skewness$X.Kurtosis_compartment,  # y axis
              z = kurtosis_skewness$X.stdev_compartment, color=col)         # z axis
legend("topleft",legend=levels(kurtosis_skewness$correlation), fill =mycolors, cex=0.8, box.lty=0)

p_ccd<-subset(kurtosis_skewness, kurtosis_skewness$Gene=="CCNB1"| kurtosis_skewness$Gene=="GATA6"| kurtosis_skewness$Gene=="NECAP2")       
mycolors = c('royal blue','grey52')
kurtosis_skewness$col<-mycolors[kurtosis_skewness$correlation]
kurtosis_skewness<-data.frame(ccd_notccd_features_compartments[,c(1,4,10,51,52,53)])
kurtosis_skewness$X.Kurtosis_compartment<-as.numeric(kurtosis_skewness$X.Kurtosis_compartment)
kurtosis_skewness$skewness_compartment<-as.numeric(kurtosis_skewness$skewness_compartment)
#kurtosis_skewness$X.stdev_compartment<-as.numeric(kurtosis_skewness$X.stdev_compartment)
kurtosis_skewness$corr_FDR_explainedvar<-factor(kurtosis_skewness$corr_FDR_explainedvar)
par(mar=c(5,5,4,4))
plot(kurtosis_skewness$skewness_compartment,kurtosis_skewness$X.Kurtosis_compartment,  col=mycolors[kurtosis_skewness$corr_FDR_explainedvar],frame.plot=FALSE, ylab="Kurtosis",xlab="skewness", las=1, xlim=c(0,15) )
points(p_ccd$skewness_compartment, p_ccd$X.Kurtosis_compartment, pch="*", cex=2, col="darkblue")
text(0.25,16,"GATA6", col="royal blue", cex=0.8)
text(0.955,24,"CCNB1", col="royal blue", cex=0.8)
text(3.9345,45,"NECAP2", col="royal blue", cex=0.8)
legend(0,150,legend=c("CCD","Non CCD"), fill =mycolors, cex=0.8, box.lty=0)

set.seed(1)
mycolors = c('"lightblue"','lightpink',"palegreen3")
mycolors = c("slateblue1","orchid3","tan1")
Kmean_kurtosis_skewness<-data.frame(kurtosis_skewness[,c(7,8)])
kurtosis<-as.numeric(Kmean_kurtosis_skewness$X.Kurtosis_compartment)
skewness<-as.numeric(Kmean_kurtosis_skewness$skewness_compartment)
Kmean_kurtosis_skewness<-data.frame(kurtosis, skewness)
results<-kmeans(Kmean_kurtosis_skewness,2)
results$cluster<-factor(results$cluster)
plot(skewness, kurtosis, col=mycolors[results$cluster],frame.plot=FALSE, ylab="Kurtosis",xlab="Skewness", las=1,xlim=c(0,15))
results$cluster<-factor(results$cluster)
legend("bottomright",legend=levels(results$cluster), fill =mycolors, cex=0.8)
plot(skewness, kurtosis, col=kurtosis_skewness$corr_FDR_explainedvar,frame.plot=FALSE, ylab="Kurtosis",xlab="Skewness", las=1 )
table(kurtosis_skewness$correlation,results$cluster)
require(ColorPalette)
monoPalette(3, name = c("slateblue1", "orchid3", "tan1"))
col=colorRampPalette(c("slateblue1", "orchid3", "tan1"))

mycolors = c("slateblue1","orchid3","tan1")

skewness_kurtosis_clustered<-data.frame(kurtosis_skewness, results$cluster)
skewness_kurtosis_clustered_ccd<-subset(skewness_kurtosis_clustered,skewness_kurtosis_clustered$correlation=="ccd")
skewness_kurtosis_clustered_nonccd<-subset(skewness_kurtosis_clustered,skewness_kurtosis_clustered$correlation=="notccd")
par(mar=c(5,5,2,8))
plot(skewness_kurtosis_clustered_ccd$skewness_compartment,skewness_kurtosis_clustered_ccd$X.Kurtosis_compartment, col=c("slateblue1","orchid3","tan1")[as.numeric(skewness_kurtosis_clustered_ccd$results.cluster)],pch=1,cex=0.8,frame.plot=FALSE, ylab="Kurtosis",xlab="Skewness", las=1,xlim=c(0,15))
points(skewness_kurtosis_clustered_nonccd$skewness_compartment,skewness_kurtosis_clustered_nonccd$X.Kurtosis_compartment, col=c("slateblue1","orchid3","tan1")[as.numeric(skewness_kurtosis_clustered_nonccd$results.cluster)],pch="*", cex=1.2)
legend(13,35,legend=c("1","2","3"), fill =c("tan1","slateblue1","orchid3"), cex=0.8, box.lty=0)
legend(9,30,legend=c("CCD","Non CCD"), pch=c(1,8), col=c("black", "black"), cex=0.8, box.lty=0)
points(p_ccd$skewness_compartment, p_ccd$X.Kurtosis_compartment, pch=1, cex=0.8, col="gray")
text(0.27,16,"GATA6", col="tan1", cex=0.8)
text(0.955,24,"CCNB1", col="tan1", cex=0.8)
text(3.9345,45,"NECAP2", col="slateblue1", cex=0.8)


#plot(skewness_kurtosis_clustered_nonccd$skewness_compartment,skewness_kurtosis_clustered_nonccd$X.Kurtosis_compartment, col=c("slateblue1","orchid3","tan1")[as.numeric(skewness_kurtosis_clustered_nonccd$results.cluster)],pch="*",frame.plot=FALSE, ylab="Kurtosis",xlab="Skewness", las=1,xlim=c(0,15))
#require(cluster)
#clusplot(Kmean_kurtosis_skewness,results$cluste, color=TRUE, shade=FALSE, labels=2, lines=0 )

#####elbow method to choose number of clusters
k.max <- 15
data <- Kmean_kurtosis_skewness
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


re<-cbind(results$cluster, kurtosis_skewness)
