workingdir<-setwd("/Users/diana.telessemian/Desktop/gaussian_per_plate/FUCCI_staining_review/")


#######################
######all plates together  
nuc6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6717<-cbind(nuc6717,well=rep(nuc6717$ImageNumber))
nuc6717<-cbind(nuc6717,plate=6717)
Image_data_6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
nuc6717$well<- Image_data_6717$Metadata_Well[match(nuc6717$well, Image_data_6717$Group_Index)]
nuc6717$site<-nuc6717$ImageNumber
nuc6717$site<- Image_data_6717$Metadata_Site[match(nuc6717$site, Image_data_6717$Group_Index)]
nuc6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6718<-cbind(nuc6718,well=rep(nuc6718$ImageNumber))
nuc6718<-cbind(nuc6718,plate=6718)
Image_data_6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
nuc6718$well<- Image_data_6718$Metadata_Well[match(nuc6718$well, Image_data_6718$Group_Index)]
nuc6718$site<-nuc6718$ImageNumber
nuc6718$site<- Image_data_6718$Metadata_Site[match(nuc6718$site, Image_data_6718$Group_Index)]
nuc6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6719<-cbind(nuc6719,well=rep(nuc6719$ImageNumber))
nuc6719<-cbind(nuc6719,plate=6719)
Image_data_6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6719$well<- Image_data_6719$Metadata_Well[match(nuc6719$well, Image_data_6719$Group_Index)]
nuc6719$site<- nuc6719$ImageNumber
nuc6719$site<- Image_data_6719$Metadata_Site[match(nuc6719$site, Image_data_6719$Group_Index)]
nuc6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6720<-cbind(nuc6720,well=rep(nuc6720$ImageNumber))
nuc6720<-cbind(nuc6720,plate=6720)
Image_data_6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
nuc6720$well<- Image_data_6720$Metadata_Well[match(nuc6720$well, Image_data_6720$Group_Index)]
nuc6720$site<-nuc6720$ImageNumber
nuc6720$site<- Image_data_6720$Metadata_Site[match(nuc6720$site, Image_data_6720$Group_Index)]
nuc6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6721<-cbind(nuc6721,well=rep(nuc6721$ImageNumber))
nuc6721<-cbind(nuc6721,plate=6721)
Image_data_6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6721$well<- Image_data_6721$Metadata_Well[match(nuc6721$well, Image_data_6721$Group_Index)]
nuc6721$site<-nuc6721$ImageNumber
nuc6721$site<- Image_data_6721$Metadata_Site[match(nuc6721$site, Image_data_6721$Group_Index)]
nuc6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6722<-cbind(nuc6722,well=rep(nuc6722$ImageNumber))
nuc6722<-cbind(nuc6722,plate=6722)
Image_data_6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6722$well<- Image_data_6722$Metadata_Well[match(nuc6722$well, Image_data_6722$Group_Index)]
nuc6722$site<-nuc6722$ImageNumber
nuc6722$site<- Image_data_6722$Metadata_Site[match(nuc6722$site, Image_data_6722$Group_Index)]
nuc6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6724<-cbind(nuc6724,well=rep(nuc6724$ImageNumber))
nuc6724<-cbind(nuc6724,plate=6724)
Image_data_6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6724$well<- Image_data_6724$Metadata_Well[match(nuc6724$well, Image_data_6724$Group_Index)]
nuc6724$site<-nuc6724$ImageNumber
nuc6724$site<- Image_data_6724$Metadata_Site[match(nuc6724$site, Image_data_6724$Group_Index)]
nuc6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6725<-cbind(nuc6725,well=rep(nuc6725$ImageNumber))
nuc6725<-cbind(nuc6725,plate=6725)
Image_data_6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6725$well<- Image_data_6725$Metadata_Well[match(nuc6725$well, Image_data_6725$Group_Index)]
nuc6725$site<-nuc6725$ImageNumber
nuc6725$site<- Image_data_6725$Metadata_Site[match(nuc6725$site, Image_data_6725$Group_Index)]
nuc6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6731<-cbind(nuc6731,well=rep(nuc6731$ImageNumber))
nuc6731<-cbind(nuc6731,plate=6731)
Image_data_6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6731$well<- Image_data_6731$Metadata_Well[match(nuc6731$well, Image_data_6731$Group_Index)]
nuc6731$site<-nuc6731$ImageNumber
nuc6731$site<- Image_data_6731$Metadata_Site[match(nuc6731$site, Image_data_6731$Group_Index)]
nuc6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6734<-cbind(nuc6734,well=rep(nuc6734$ImageNumber))
nuc6734<-cbind(nuc6734,plate=6734)
Image_data_6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6734$well<- Image_data_6734$Metadata_Well[match(nuc6734$well, Image_data_6734$Group_Index)]
nuc6734$site<-nuc6734$ImageNumber
nuc6734$site<- Image_data_6734$Metadata_Site[match(nuc6734$site, Image_data_6734$Group_Index)]
nuc6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6735<-cbind(nuc6735,well=rep(nuc6735$ImageNumber))
nuc6735<-cbind(nuc6735,plate=6735)
Image_data_6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6735$well<- Image_data_6735$Metadata_Well[match(nuc6735$well, Image_data_6735$Group_Index)]
nuc6735$site<-nuc6735$ImageNumber
nuc6735$site<- Image_data_6735$Metadata_Site[match(nuc6735$site, Image_data_6735$Group_Index)]
nuc6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6736<-cbind(nuc6736,well=rep(nuc6736$ImageNumber))
nuc6736<-cbind(nuc6736,plate=6736)
Image_data_6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6736$well<- Image_data_6736$Metadata_Well[match(nuc6736$well, Image_data_6736$Group_Index)]
nuc6736$site<-nuc6736$ImageNumber
nuc6736$site<- Image_data_6736$Metadata_Site[match(nuc6736$site, Image_data_6736$Group_Index)]
nuc6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Nuclei.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,23,24,25,30,31,32,33,58,59,60,61,62,63,64,65)]
nuc6745<-cbind(nuc6745,well=rep(nuc6745$ImageNumber))
nuc6745<-cbind(nuc6745,plate=6745)
Image_data_6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
nuc6745$well<- Image_data_6745$Metadata_Well[match(nuc6745$well, Image_data_6745$Group_Index)]
nuc6745$site<-nuc6745$ImageNumber
nuc6745$site<- Image_data_6745$Metadata_Site[match(nuc6745$site, Image_data_6745$Group_Index)]
nuc_all<-rbind(nuc6717,nuc6718,nuc6719,nuc6720,nuc6721,nuc6722,nuc6724,nuc6725,nuc6731,nuc6734,nuc6735,nuc6736,nuc6745)
G1nuc <-nuc_all[nuc_all$Children_FilteredG1_Count=="1", ]
G1Snuc<-nuc_all[nuc_all$Children_FilteredG1S_Count=="1", ]
SG2nuc<-nuc_all[nuc_all$Children_FilteredSG2_Count=="1", ]
nuc_all[nuc_all$Children_FilteredG1S_Count=="1",5 ]<-"2"
nuc_all[nuc_all$Children_FilteredSG2_Count=="1",7 ]<-"3"
G1num<-as.numeric(nuc_all$Children_FilteredG1_Count)
G1Snum<-as.numeric(nuc_all$Children_FilteredG1S_Count)
SG2num<-as.numeric(nuc_all$Children_FilteredSG2_Count)
CCnuc_all <- within(nuc_all, CC <- paste(G1num+G1Snum+SG2num, sep=''))
CCdatanuc_phases<-CCnuc_all[CCnuc_all$CC=="1",23 ] <-"G1"
CCdatanuc_phases<-CCnuc_all[CCnuc_all$CC=="2",23 ] <-"G1S"
CCdatanuc_phases<-CCnuc_all[CCnuc_all$CC=="3",23 ] <-"SG2"
write.table(CCnuc_all, file="CCnuc_all.csv", sep = ",",qmethod = "double")
#########cytoplasm
cyto6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6717<-cbind(cyto6717,well=rep(cyto6717$ImageNumber))
cyto6717<-cbind(cyto6717,plate=6717)
Image_data_6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cyto6717$well<- Image_data_6717$Metadata_Well[match(cyto6717$well, Image_data_6717$Group_Index)]
cyto6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6718<-cbind(cyto6718,well=rep(cyto6718$ImageNumber))
cyto6718<-cbind(cyto6718,plate=6718)
Image_data_6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cyto6718$well<- Image_data_6718$Metadata_Well[match(cyto6718$well, Image_data_6718$Group_Index)]
cyto6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6719<-cbind(cyto6719,well=rep(cyto6719$ImageNumber))
cyto6719<-cbind(cyto6719,plate=6719)
Image_data_6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6719$well<- Image_data_6719$Metadata_Well[match(cyto6719$well, Image_data_6719$Group_Index)]
cyto6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6720<-cbind(cyto6720,well=rep(cyto6720$ImageNumber))
cyto6720<-cbind(cyto6720,plate=6720)
Image_data_6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cyto6720$well<- Image_data_6720$Metadata_Well[match(cyto6720$well, Image_data_6720$Group_Index)]
cyto6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6721<-cbind(cyto6721,well=rep(cyto6721$ImageNumber))
cyto6721<-cbind(cyto6721,plate=6721)
Image_data_6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6721$well<- Image_data_6721$Metadata_Well[match(cyto6721$well, Image_data_6721$Group_Index)]
cyto6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6722<-cbind(cyto6722,well=rep(cyto6722$ImageNumber))
cyto6722<-cbind(cyto6722,plate=6722)
Image_data_6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6722$well<- Image_data_6722$Metadata_Well[match(cyto6722$well, Image_data_6722$Group_Index)]
cyto6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6724<-cbind(cyto6724,well=rep(cyto6724$ImageNumber))
cyto6724<-cbind(cyto6724,plate=6724)
Image_data_6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6724$well<- Image_data_6724$Metadata_Well[match(cyto6724$well, Image_data_6724$Group_Index)]
cyto6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6725<-cbind(cyto6725,well=rep(cyto6725$ImageNumber))
cyto6725<-cbind(cyto6725,plate=6725)
Image_data_6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6725$well<- Image_data_6725$Metadata_Well[match(cyto6725$well, Image_data_6725$Group_Index)]
cyto6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6731<-cbind(cyto6731,well=rep(cyto6731$ImageNumber))
cyto6731<-cbind(cyto6731,plate=6731)
Image_data_6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6731$well<- Image_data_6731$Metadata_Well[match(cyto6731$well, Image_data_6731$Group_Index)]
cyto6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6734<-cbind(cyto6734,well=rep(cyto6734$ImageNumber))
cyto6734<-cbind(cyto6734,plate=6734)
Image_data_6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6734$well<- Image_data_6734$Metadata_Well[match(cyto6734$well, Image_data_6734$Group_Index)]
cyto6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6735<-cbind(cyto6735,well=rep(cyto6735$ImageNumber))
cyto6735<-cbind(cyto6735,plate=6735)
Image_data_6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6735$well<- Image_data_6735$Metadata_Well[match(cyto6735$well, Image_data_6735$Group_Index)]
cyto6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6736<-cbind(cyto6736,well=rep(cyto6736$ImageNumber))
cyto6736<-cbind(cyto6736,plate=6736)
Image_data_6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6736$well<- Image_data_6736$Metadata_Well[match(cyto6736$well, Image_data_6736$Group_Index)]
cyto6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Cytoplasm.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,25,26,27,28,53,54,55,56,57,58,59,60)]
cyto6745<-cbind(cyto6745,well=rep(cyto6745$ImageNumber))
cyto6745<-cbind(cyto6745,plate=6745)
Image_data_6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cyto6745$well<- Image_data_6745$Metadata_Well[match(cyto6745$well, Image_data_6745$Group_Index)]

cyto_all<-rbind(cyto6717,cyto6718,cyto6719,cyto6720,cyto6721,cyto6722,cyto6724,cyto6725,cyto6731,cyto6734,cyto6735,cyto6736,cyto6745)
write.table(cyto_all, file="cyto_all_plates.csv", sep = ",",qmethod = "double")
names(cyto_all) <- c("ImageNumber", "ObjectNumber","Area_cyto","perimeter_cyto","Integrated_green_cyto","Integrated_red_cyto", "Integrated_ab_cyto","Integrated_mt_cyto" ,"Mean_Green_Fucci_Cyto","Mean_Red_Fucci_Cyto", "Mean_ab_Cyto", "Mean_mt_Cyto","Median_Green_Fucci_Cyto","Median_Red_Fucci_Cyto","Median_ab_Cyto","Median_mt_Cyto","well","plate" )





#########cellplasm
cell6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6717<-cbind(cell6717,well=rep(cell6717$ImageNumber))
cell6717<-cbind(cell6717,plate=6717)
Image_data_6717 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6717/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cell6717$well<- Image_data_6717$Metadata_Well[match(cell6717$well, Image_data_6717$Group_Index)]
cell6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6718<-cbind(cell6718,well=rep(cell6718$ImageNumber))
cell6718<-cbind(cell6718,plate=6718)
Image_data_6718 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6718/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cell6718$well<- Image_data_6718$Metadata_Well[match(cell6718$well, Image_data_6718$Group_Index)]
cell6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6719<-cbind(cell6719,well=rep(cell6719$ImageNumber))
cell6719<-cbind(cell6719,plate=6719)
Image_data_6719 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6719/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6719$well<- Image_data_6719$Metadata_Well[match(cell6719$well, Image_data_6719$Group_Index)]
cell6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6720<-cbind(cell6720,well=rep(cell6720$ImageNumber))
cell6720<-cbind(cell6720,plate=6720)
Image_data_6720 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6720/measurement/Image.csv", header=T, stringsAsFactors = FALSE)[,c(46,47,537,538)]
cell6720$well<- Image_data_6720$Metadata_Well[match(cell6720$well, Image_data_6720$Group_Index)]
cell6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6721<-cbind(cell6721,well=rep(cell6721$ImageNumber))
cell6721<-cbind(cell6721,plate=6721)
Image_data_6721<- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6721/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6721$well<- Image_data_6721$Metadata_Well[match(cell6721$well, Image_data_6721$Group_Index)]
cell6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6722<-cbind(cell6722,well=rep(cell6722$ImageNumber))
cell6722<-cbind(cell6722,plate=6722)
Image_data_6722 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6722/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6722$well<- Image_data_6722$Metadata_Well[match(cell6722$well, Image_data_6722$Group_Index)]
cell6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6724<-cbind(cell6724,well=rep(cell6724$ImageNumber))
cell6724<-cbind(cell6724,plate=6724)
Image_data_6724 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6724/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6724$well<- Image_data_6724$Metadata_Well[match(cell6724$well, Image_data_6724$Group_Index)]
cell6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6725<-cbind(cell6725,well=rep(cell6725$ImageNumber))
cell6725<-cbind(cell6725,plate=6725)
Image_data_6725 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6725/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6725$well<- Image_data_6725$Metadata_Well[match(cell6725$well, Image_data_6725$Group_Index)]
cell6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6731<-cbind(cell6731,well=rep(cell6731$ImageNumber))
cell6731<-cbind(cell6731,plate=6731)
Image_data_6731 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6731/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6731$well<- Image_data_6731$Metadata_Well[match(cell6731$well, Image_data_6731$Group_Index)]
cell6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6734<-cbind(cell6734,well=rep(cell6734$ImageNumber))
cell6734<-cbind(cell6734,plate=6734)
Image_data_6734 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6734/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6734$well<- Image_data_6734$Metadata_Well[match(cell6734$well, Image_data_6734$Group_Index)]
cell6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6735<-cbind(cell6735,well=rep(cell6735$ImageNumber))
cell6735<-cbind(cell6735,plate=6735)
Image_data_6735 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6735/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6735$well<- Image_data_6735$Metadata_Well[match(cell6735$well, Image_data_6735$Group_Index)]
cell6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6736<-cbind(cell6736,well=rep(cell6736$ImageNumber))
cell6736<-cbind(cell6736,plate=6736)
Image_data_6736 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6736/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6736$well<- Image_data_6736$Metadata_Well[match(cell6736$well, Image_data_6736$Group_Index)]
cell6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Cells.csv", header=T, stringsAsFactors = FALSE)[,c(1,2,3,19,26,27,28,29,54,55,56,57,58,59,60,61)]
cell6745<-cbind(cell6745,well=rep(cell6745$ImageNumber))
cell6745<-cbind(cell6745,plate=6745)
Image_data_6745 <- read.csv ("/Volumes/DIANASDRIVE/Greiner SensoPLate 655891 96 40x 4 channels/Cell_profiler_output/Analysis/Plate_6745/Image.csv", header=T, stringsAsFactors = FALSE)[,c(52,53,547,549)]
cell6745$well<- Image_data_6745$Metadata_Well[match(cell6745$well, Image_data_6745$Group_Index)]


cell_all<-rbind(cell6717,cell6718,cell6719,cell6720,cell6721,cell6722,cell6724,cell6725,cell6731,cell6734,cell6735,cell6736,cell6745)
write.table(cell_all, file="cell_all_plates.csv", sep = ",",qmethod = "double")
names(cell_all) <-  c("ImageNumber", "ObjectNumber","Area_cell","perimeter_cell","Integrated_green_cell","Integrated_red_cell", "Integrated_ab_cell","Integrated_mt_cell","Mean_Green_Fucci_cell","Mean_Red_Fucci_cell", "Mean_ab_cell", "Mean_mt_cell","Median_Green_Fucci_cell","Median_Red_Fucci_cell","Median_ab_cell","Median_mt_cell","well","plate" )




########merging nucleus and cyto and cell data together 
nuc_cyto<-merge(cyto_all, CCnuc_all,by=c('ImageNumber','ObjectNumber',"plate","well"), all=TRUE)
nuc_cyto$well_plate <- paste(nuc_cyto$well,nuc_cyto$plate, sep="_")
nuc_cyto_all<-merge(nuc_cyto, cell_all,by=c('ImageNumber','ObjectNumber',"plate","well"), all=TRUE)

########remove out of focus images 
nuc_cyto_all$well_plate_imagenb <- paste(nuc_cyto_all$well_plate,nuc_cyto_all$ImageNumber, sep="_")
nuc_cyto_all_excl_outoffocus<-nuc_cyto_all[nuc_cyto_all$well_plate_imagenb!= "H07_55185977_543" 
                                           &nuc_cyto_all$well_plate_imagenb!='F10_55185977_418'
                                           &nuc_cyto_all$well_plate_imagenb!='D02_55185977_225'
                                           &nuc_cyto_all$well_plate_imagenb!='F09_55185977_412'
                                           &nuc_cyto_all$well_plate_imagenb!='A06_55185977_35'
                                           &nuc_cyto_all$well_plate_imagenb!='D05_55185977_246'
                                           &nuc_cyto_all$well_plate_imagenb!='B07_55185977_109'
                                           &nuc_cyto_all$well_plate_imagenb!='F02_55185977_372'
                                           &nuc_cyto_all$well_plate_imagenb!='F02_55185977_368'
                                           &nuc_cyto_all$well_plate_imagenb!='E01_55185977_294'
                                           &nuc_cyto_all$well_plate_imagenb!='A05_55185977_25'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55185977_49'
                                           &nuc_cyto_all$well_plate_imagenb!='F02_55195978_370'
                                           &nuc_cyto_all$well_plate_imagenb!='H04_55195978_523'
                                           &nuc_cyto_all$well_plate_imagenb!='C02_55195978_151'
                                           &nuc_cyto_all$well_plate_imagenb!='F12_55195978_427'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_55205980_202'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55205980_49'
                                           &nuc_cyto_all$well_plate_imagenb!='A10_55205980_58'
                                           &nuc_cyto_all$well_plate_imagenb!='E08_55205980_331'
                                           &nuc_cyto_all$well_plate_imagenb!='D04_55205980_240'
                                           &nuc_cyto_all$well_plate_imagenb!='D01_55215982_220'
                                           &nuc_cyto_all$well_plate_imagenb!='B09_55225983_124'
                                           &nuc_cyto_all$well_plate_imagenb!='B09_55225983_123'
                                           &nuc_cyto_all$well_plate_imagenb!='B09_55225983_121'
                                           &nuc_cyto_all$well_plate_imagenb!='B09_55225983_122'
                                           &nuc_cyto_all$well_plate_imagenb!='G09_55225983_481'
                                           &nuc_cyto_all$well_plate_imagenb!='D01_55225983_217'
                                           &nuc_cyto_all$well_plate_imagenb!='F10_55225983_415'
                                           &nuc_cyto_all$well_plate_imagenb!='B01_55225983_73'
                                           &nuc_cyto_all$well_plate_imagenb!='D04_55225983_238'
                                           &nuc_cyto_all$well_plate_imagenb!='G11_55225983_493'
                                           &nuc_cyto_all$well_plate_imagenb!='H04_55225983_526'
                                           &nuc_cyto_all$well_plate_imagenb!='C02_55225983_156'
                                           &nuc_cyto_all$well_plate_imagenb!='A03_55235979_15'
                                           &nuc_cyto_all$well_plate_imagenb!='A01_55235979_1'
                                           &nuc_cyto_all$well_plate_imagenb!='G07_55235979_470'
                                           &nuc_cyto_all$well_plate_imagenb!='A02_55235979_12'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_516'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_515'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_514'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_513'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_512'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55245981_511'
                                           &nuc_cyto_all$well_plate_imagenb!='C02_55245981_152'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55245981_52'
                                           &nuc_cyto_all$well_plate_imagenb!='D01_55245981_217'
                                           &nuc_cyto_all$well_plate_imagenb!='G06_55335984_465'
                                           &nuc_cyto_all$well_plate_imagenb!='F06_55335984_393'
                                           &nuc_cyto_all$well_plate_imagenb!='G12_55345985_503'
                                           &nuc_cyto_all$well_plate_imagenb!='E10_55345985_344'
                                           &nuc_cyto_all$well_plate_imagenb!='B04_55345985_95'
                                           &nuc_cyto_all$well_plate_imagenb!='B04_55345985_92'
                                           &nuc_cyto_all$well_plate_imagenb!='B04_55345985_96'
                                           &nuc_cyto_all$well_plate_imagenb!='F07_55355986_397'
                                           &nuc_cyto_all$well_plate_imagenb!='D04_55355986_235'
                                           &nuc_cyto_all$well_plate_imagenb!='A07_55365987_37'
                                           &nuc_cyto_all$well_plate_imagenb!='D10_55365987_274'
                                           &nuc_cyto_all$well_plate_imagenb!='A03_55365987_13'
                                           &nuc_cyto_all$well_plate_imagenb!='A01_55365987_2'
                                           &nuc_cyto_all$well_plate_imagenb!='B10_55365987_130'
                                           &nuc_cyto_all$well_plate_imagenb!='B08_55365987_118'
                                           &nuc_cyto_all$well_plate_imagenb!='F11_55375988_424'
                                           &nuc_cyto_all$well_plate_imagenb!='D03_55385989_230'
                                           &nuc_cyto_all$well_plate_imagenb!='D05_55385989_241'
                                           &nuc_cyto_all$well_plate_imagenb!='F05_55385989_385'
                                           &nuc_cyto_all$well_plate_imagenb!='D03_55385989_229'
                                           &nuc_cyto_all$well_plate_imagenb!='H03_55385989_518'
                                           &nuc_cyto_all$well_plate_imagenb!='G01_55385989_435'
                                           &nuc_cyto_all$well_plate_imagenb!='E12_55385989_355'
                                           &nuc_cyto_all$well_plate_imagenb!='C04_55385989_166'
                                           &nuc_cyto_all$well_plate_imagenb!='E12_55385989_358'
                                           &nuc_cyto_all$well_plate_imagenb!='D04_55395990_238'
                                           &nuc_cyto_all$well_plate_imagenb!='H04_55395990_528'
                                           &nuc_cyto_all$well_plate_imagenb!='F05_55395990_385'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55395990_514'
                                           &nuc_cyto_all$well_plate_imagenb!='E02_55395990_298'
                                           &nuc_cyto_all$well_plate_imagenb!='F06_55395990_392'
                                           &nuc_cyto_all$well_plate_imagenb!='F10_55395990_416'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55405991_50'
                                           &nuc_cyto_all$well_plate_imagenb!='B01_55405991_73'
                                           &nuc_cyto_all$well_plate_imagenb!='B05_55405991_97'
                                           &nuc_cyto_all$well_plate_imagenb!='B10_55405991_132'
                                           &nuc_cyto_all$well_plate_imagenb!='D11_55405991_277'
                                           &nuc_cyto_all$well_plate_imagenb!='F04_55405991_381'
                                           &nuc_cyto_all$well_plate_imagenb!='G07_55405991_471'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55195978_51'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55195978_52'
                                           &nuc_cyto_all$well_plate_imagenb!='A09_55195978_54'
                                           &nuc_cyto_all$well_plate_imagenb!='F12_55195978_432'
                                           &nuc_cyto_all$well_plate_imagenb!='C12_55195978_211'
                                           &nuc_cyto_all$well_plate_imagenb!='D02_55195978_223'
                                           &nuc_cyto_all$well_plate_imagenb!='A12_55205980_67'
                                           &nuc_cyto_all$well_plate_imagenb!='F12_55215982_432'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_55225983_516'
                                           &nuc_cyto_all$well_plate_imagenb!='F05_55225983_386'
                                           &nuc_cyto_all$well_plate_imagenb!='B12_55235979_140'
                                           &nuc_cyto_all$well_plate_imagenb!='A11_55245981_64'
                                           &nuc_cyto_all$well_plate_imagenb!='A11_55245981_65'
                                           &nuc_cyto_all$well_plate_imagenb!='A11_55245981_66'
                                           &nuc_cyto_all$well_plate_imagenb!='E06_55355986_319'
                                           &nuc_cyto_all$well_plate_imagenb!='E06_55355986_319'
                                           &nuc_cyto_all$well_plate_imagenb!="D05_55185977_243"
                                           &nuc_cyto_all$well_plate_imagenb!="D05_55185977_244"
                                           &nuc_cyto_all$well_plate_imagenb!="D05_55185977_246"
                                           &nuc_cyto_all$well_plate_imagenb!="B02_55195978_79"
                                           &nuc_cyto_all$well_plate_imagenb!="C02_55195978_151"
                                           &nuc_cyto_all$well_plate_imagenb!="A05_55195978_25"
                                           &nuc_cyto_all$well_plate_imagenb!="F06_55195978_391"
                                           &nuc_cyto_all$well_plate_imagenb!="F06_55195978_395"
                                           &nuc_cyto_all$well_plate_imagenb!="A07_55195978_40"
                                           &nuc_cyto_all$well_plate_imagenb!="A07_55195978_41"
                                           &nuc_cyto_all$well_plate_imagenb!="A09_55195978_39"
                                           &nuc_cyto_all$well_plate_imagenb!="A09_55195978_40"
                                           &nuc_cyto_all$well_plate_imagenb!="A09_55195978_42"
                                           &nuc_cyto_all$well_plate_imagenb!="E08_55235979_331"
                                           &nuc_cyto_all$well_plate_imagenb!="D02_55205980_224"
                                           &nuc_cyto_all$well_plate_imagenb!="D02_55205980_225"
                                           &nuc_cyto_all$well_plate_imagenb!="B12_55205980_139"
                                           &nuc_cyto_all$well_plate_imagenb!="C03_55245981_160"
                                           &nuc_cyto_all$well_plate_imagenb!="C03_55245981_161"
                                           &nuc_cyto_all$well_plate_imagenb!="C03_55245981_162"
                                           &nuc_cyto_all$well_plate_imagenb!="D06_55245981_250"
                                           &nuc_cyto_all$well_plate_imagenb!="D06_55245981_251"
                                           &nuc_cyto_all$well_plate_imagenb!="D06_55245981_252"
                                           &nuc_cyto_all$well_plate_imagenb!="C08_55245981_187"
                                           &nuc_cyto_all$well_plate_imagenb!="F12_55245981_430"
                                           &nuc_cyto_all$well_plate_imagenb!="F12_55245981_431"
                                           &nuc_cyto_all$well_plate_imagenb!="F12_55245981_432"
                                           &nuc_cyto_all$well_plate_imagenb!="B01_55215982_73"
                                           &nuc_cyto_all$well_plate_imagenb!="B01_55215982_74"
                                           &nuc_cyto_all$well_plate_imagenb!="D10_55215982_274"
                                           &nuc_cyto_all$well_plate_imagenb!="D10_55215982_275"
                                           &nuc_cyto_all$well_plate_imagenb!="D10_55215982_276"
                                           &nuc_cyto_all$well_plate_imagenb!="C12_55215982_211"
                                           &nuc_cyto_all$well_plate_imagenb!="G09_55225983_481"
                                           &nuc_cyto_all$well_plate_imagenb!="H08_55335984_549"
                                           &nuc_cyto_all$well_plate_imagenb!="H08_55335984_550"
                                           &nuc_cyto_all$well_plate_imagenb!="H08_55335984_552"
                                           &nuc_cyto_all$well_plate_imagenb!="B09_55335984_123"
                                           &nuc_cyto_all$well_plate_imagenb!="C06_55345985_175"
                                           &nuc_cyto_all$well_plate_imagenb!="B10_55345985_131"
                                           &nuc_cyto_all$well_plate_imagenb!="B10_55345985_132"
                                           &nuc_cyto_all$well_plate_imagenb!="A01_55365987_2"
                                           &nuc_cyto_all$well_plate_imagenb!="G06_55365987_463"
                                           &nuc_cyto_all$well_plate_imagenb!="C01_55375988_148"
                                           &nuc_cyto_all$well_plate_imagenb!="C01_55375988_149"
                                           &nuc_cyto_all$well_plate_imagenb!="C01_55375988_150"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55375988_517"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55375988_520"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55375988_521"
                                           &nuc_cyto_all$well_plate_imagenb!="D02_55385989_223"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55385989_518"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55385989_521"
                                           &nuc_cyto_all$well_plate_imagenb!="H03_55385989_522"
                                           &nuc_cyto_all$well_plate_imagenb!="E04_55385989_307"
                                           &nuc_cyto_all$well_plate_imagenb!="E04_55385989_312"
                                           &nuc_cyto_all$well_plate_imagenb!="E12_55385989_355"
                                           &nuc_cyto_all$well_plate_imagenb!="E12_55385989_358"
                                           &nuc_cyto_all$well_plate_imagenb!="D03_55395990_232"
                                           &nuc_cyto_all$well_plate_imagenb!="D03_55395990_233"
                                           &nuc_cyto_all$well_plate_imagenb!="D04_55395990_238"
                                           &nuc_cyto_all$well_plate_imagenb!="B07_55395990_112"
                                           &nuc_cyto_all$well_plate_imagenb!="G08_55395990_478"
                                           &nuc_cyto_all$well_plate_imagenb!="G08_55395990_479"
                                           &nuc_cyto_all$well_plate_imagenb!="A10_55405991_56"
                                           &nuc_cyto_all$well_plate_imagenb!="C10_55405991_202"
                                           &nuc_cyto_all$well_plate_imagenb!="A11_55405991_64"
                                           &nuc_cyto_all$well_plate_imagenb!="D11_55405991_277"
                                           &nuc_cyto_all$well_plate_imagenb!="A12_55405991_69"
                                           &nuc_cyto_all$well_plate_imagenb!="A12_55405991_72"
                                           &nuc_cyto_all$well_plate_imagenb!='A04_75836284_21'
                                           &nuc_cyto_all$well_plate_imagenb!='A06_75836284_35'
                                           &nuc_cyto_all$well_plate_imagenb!='B10_75836284_128'
                                           &nuc_cyto_all$well_plate_imagenb!='B10_75836284_131'
                                           &nuc_cyto_all$well_plate_imagenb!='B10_75836284_132'
                                           &nuc_cyto_all$well_plate_imagenb!='B11_75836284_133'
                                           &nuc_cyto_all$well_plate_imagenb!='B11_75836284_135'
                                           &nuc_cyto_all$well_plate_imagenb!='B11_75836284_136'
                                           &nuc_cyto_all$well_plate_imagenb!='C04_75836284_166'
                                           &nuc_cyto_all$well_plate_imagenb!='C08_75836284_187'
                                           &nuc_cyto_all$well_plate_imagenb!='C08_75836284_188'
                                           &nuc_cyto_all$well_plate_imagenb!='C08_75836284_191'
                                           &nuc_cyto_all$well_plate_imagenb!='D07_75836284_254'
                                           &nuc_cyto_all$well_plate_imagenb!='D09_75836284_265'
                                           &nuc_cyto_all$well_plate_imagenb!='D09_75836284_266'
                                           &nuc_cyto_all$well_plate_imagenb!='D09_75836284_268'
                                           &nuc_cyto_all$well_plate_imagenb!='D12_75836284_284'
                                           &nuc_cyto_all$well_plate_imagenb!='D12_75836284_285'
                                           &nuc_cyto_all$well_plate_imagenb!='D12_75836284_287'
                                           &nuc_cyto_all$well_plate_imagenb!='E01_75836284_291'
                                           &nuc_cyto_all$well_plate_imagenb!='E01_75836284_293'
                                           &nuc_cyto_all$well_plate_imagenb!='E03_75836284_305'
                                           &nuc_cyto_all$well_plate_imagenb!='E04_75836284_307'
                                           &nuc_cyto_all$well_plate_imagenb!='E07_75836284_325'
                                           &nuc_cyto_all$well_plate_imagenb!='E07_75836284_328'
                                           &nuc_cyto_all$well_plate_imagenb!='E08_75836284_331'
                                           &nuc_cyto_all$well_plate_imagenb!='E09_75836284_340'
                                           &nuc_cyto_all$well_plate_imagenb!='E11_75836284_353'
                                           &nuc_cyto_all$well_plate_imagenb!='F05_75836284_389'
                                           &nuc_cyto_all$well_plate_imagenb!='F05_75836284_390'
                                           &nuc_cyto_all$well_plate_imagenb!='F07_75836284_398'
                                           &nuc_cyto_all$well_plate_imagenb!='F08_75836284_407'
                                           &nuc_cyto_all$well_plate_imagenb!='F08_75836284_408'
                                           &nuc_cyto_all$well_plate_imagenb!='G04_75836284_451'
                                           &nuc_cyto_all$well_plate_imagenb!='G04_75836284_453'
                                           &nuc_cyto_all$well_plate_imagenb!='G12_75836284_502'
                                           &nuc_cyto_all$well_plate_imagenb!='H02_75836284_513'
                                           &nuc_cyto_all$well_plate_imagenb!='H05_75836284_531'
                                           &nuc_cyto_all$well_plate_imagenb!='H08_75836284_547'
                                           &nuc_cyto_all$well_plate_imagenb!='H08_75836284_552'
                                           &nuc_cyto_all$well_plate_imagenb!='H10_75836284_560'
                                           &nuc_cyto_all$well_plate_imagenb!='A02_75846286_7'
                                           &nuc_cyto_all$well_plate_imagenb!='A02_75846286_8'
                                           &nuc_cyto_all$well_plate_imagenb!='A02_75846286_11'
                                           &nuc_cyto_all$well_plate_imagenb!='A05_75846286_26'
                                           &nuc_cyto_all$well_plate_imagenb!='A05_75846286_28'
                                           &nuc_cyto_all$well_plate_imagenb!='B02_75846286_53'
                                           &nuc_cyto_all$well_plate_imagenb!='B07_75846286_82'
                                           &nuc_cyto_all$well_plate_imagenb!='C03_75846286_99'
                                           &nuc_cyto_all$well_plate_imagenb!='C03_75846286_100'
                                           &nuc_cyto_all$well_plate_imagenb!='C06_75846286_115'
                                           &nuc_cyto_all$well_plate_imagenb!='C06_75846286_118'
                                           &nuc_cyto_all$well_plate_imagenb!='C06_75846286_119'
                                           &nuc_cyto_all$well_plate_imagenb!='D02_75846286_128'
                                           &nuc_cyto_all$well_plate_imagenb!='D02_75846286_130'
                                           &nuc_cyto_all$well_plate_imagenb!='D02_75846286_131'
                                           &nuc_cyto_all$well_plate_imagenb!='E02_75846286_168'
                                           &nuc_cyto_all$well_plate_imagenb!='E03_75846286_170'
                                           &nuc_cyto_all$well_plate_imagenb!='E03_75846286_171'
                                           &nuc_cyto_all$well_plate_imagenb!='E03_75846286_172'
                                           &nuc_cyto_all$well_plate_imagenb!='E03_75846286_173'
                                           &nuc_cyto_all$well_plate_imagenb!='E04_75846286_177'
                                           &nuc_cyto_all$well_plate_imagenb!='E04_75846286_178'
                                           &nuc_cyto_all$well_plate_imagenb!='G03_75846286_245'
                                           &nuc_cyto_all$well_plate_imagenb!='H04_75846286_283'
                                           &nuc_cyto_all$well_plate_imagenb!='H04_75846286_284'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_193'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_194'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_195'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_196'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_197'
                                           &nuc_cyto_all$well_plate_imagenb!='C09_75836284_198'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_199'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_200'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_201'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_202'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_203'
                                           &nuc_cyto_all$well_plate_imagenb!='C10_75836284_204'
                                           ,]



nuc_cyto_all_excl_outoffocus$ratio<-nuc_cyto_all_excl_outoffocus$Area_cell/nuc_cyto_all_excl_outoffocus$AreaShape_Area



mitotic_green<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$ratio<2.2,]
a<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$ratio>2.2,]
mitotic_cytokinesis<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$AreaShape_Area<1550,]
b<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$AreaShape_Area>=1550,]
mitotic_red_green<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$ratio<2.2 |nuc_cyto_all_excl_outoffocus$AreaShape_Area<1550,] 
nuc_cyto_all_excl_outoffocus<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$ratio>2.2 & nuc_cyto_all_excl_outoffocus$AreaShape_Area>=1550,] 

#######gaussian per plate

gauss_6717<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6717",]
gauss_6718<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6718",]
gauss_6719<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6719",]
gauss_6720<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6720",]
gauss_6721<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6721",]
gauss_6722<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6722",]
gauss_6724<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6724",]
gauss_6725<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6725",]
gauss_6731<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6731",]
gauss_6734<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6734",]
gauss_6735<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6735",]
gauss_6736<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6736",]
gauss_6745<-nuc_cyto_all_excl_outoffocus[nuc_cyto_all_excl_outoffocus$plate=="6745",]


df_gauss_6717<-data.frame (as.numeric(log(gauss_6717$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6717$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6718<-data.frame (as.numeric(log(gauss_6718$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6718$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6719<-data.frame (as.numeric(log(gauss_6719$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6719$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6720<-data.frame (as.numeric(log(gauss_6720$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6720$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6721<-data.frame (as.numeric(log(gauss_6721$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6721$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6722<-data.frame (as.numeric(log(gauss_6722$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6722$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6724<-data.frame (as.numeric(log(gauss_6724$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6724$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6725<-data.frame (as.numeric(log(gauss_6725$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6725$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6731<-data.frame (as.numeric(log(gauss_6731$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6731$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6734<-data.frame (as.numeric(log(gauss_6734$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6734$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6735<-data.frame (as.numeric(log(gauss_6735$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6735$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6736<-data.frame (as.numeric(log(gauss_6736$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6736$Intensity_MeanIntensity_CorrResizedRedFUCCI))
df_gauss_6745<-data.frame (as.numeric(log(gauss_6745$Intensity_MeanIntensity_CorrResizedGreenFUCCI)),log(gauss_6745$Intensity_MeanIntensity_CorrResizedRedFUCCI))



#### Gaussian clustering for each plate 
require(mclust)


M <- 1e3

#xyMclust <- Mclust(df, G=3,initialization=list(subset=sample(1:nrow(df), size=M)))
set.seed(1)
xyMclust6717 <- Mclust(df_gauss_6717, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6717), size=M)))
set.seed(1)
xyMclust6718 <- Mclust(df_gauss_6718, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6718), size=M)))
set.seed(1)
xyMclust6719 <- Mclust(df_gauss_6719, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6719), size=M)))
set.seed(1)
xyMclust6720 <- Mclust(df_gauss_6720, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6720), size=M)))
set.seed(1)
xyMclust6721 <- Mclust(df_gauss_6721, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6721), size=M)))
set.seed(1)
xyMclust6722 <- Mclust(df_gauss_6722, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6722), size=M)))
set.seed(1)
xyMclust6724 <- Mclust(df_gauss_6724, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6724), size=M)))
set.seed(1)
xyMclust6725 <- Mclust(df_gauss_6725, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6725), size=M)))
set.seed(1)
xyMclust6731 <- Mclust(df_gauss_6731, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6731), size=M)))
set.seed(1)
xyMclust6734 <- Mclust(df_gauss_6734, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6734), size=M)))
set.seed(1)
xyMclust6735 <- Mclust(df_gauss_6735, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6735), size=M)))
set.seed(1)
xyMclust6736 <- Mclust(df_gauss_6736, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6736), size=M)))
set.seed(1)
xyMclust6745 <- Mclust(df_gauss_6745, G=3,modelNames="EEE" ,initialization=list(subset=sample(1:nrow(df_gauss_6745), size=M)))



predicted_clusters_6717<-predict(xyMclust6717,df_gauss_6717)
predicted_clus_6717<-data.frame(predicted_clusters_6717)
predicted_clust_6717<-as.matrix(predicted_clus_6717)
summary_clusters_6717<-summary(xyMclust6717, parameters=TRUE)
predicted_clus_6717<- within(predicted_clus_6717, Z1prob<- paste(pnorm(predicted_clus_6717$z.1), sep=''))
predicted_clus_6717 <- within(predicted_clus_6717, Z2prob <- paste(pnorm(predicted_clus_6717$z.2), sep=''))
predicted_clus_6717 <- within(predicted_clus_6717, Z3prob <- paste(pnorm(predicted_clus_6717$z.3), sep=''))
nuc_predicted_prob_6717<-cbind(gauss_6717, predicted_clus_6717)
nuc_predicted_prob_6717<- within(nuc_predicted_prob_6717, log_green <- paste(log(nuc_predicted_prob_6717$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6717 <- within(nuc_predicted_prob_6717, log_red <- paste(log(nuc_predicted_prob_6717$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6717$d<-nuc_predicted_prob_6717$classification
nuc_predicted_prob_6717[nuc_predicted_prob_6717$d=="1",64 ]<-"SG2"
nuc_predicted_prob_6717[nuc_predicted_prob_6717$d=="2",64 ]<-"G1S"
nuc_predicted_prob_6717[nuc_predicted_prob_6717$d=="3",64 ]<-"G1"
nuc_predicted_prob_6717$d->nuc_predicted_prob_6717$d1
nuc_predicted_prob_6717$a<-nuc_predicted_prob_6717$z.1>0.9 | nuc_predicted_prob_6717$z.2>0.9 | nuc_predicted_prob_6717$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6717) ){
  if(nuc_predicted_prob_6717$a[i]=="FALSE")
    nuc_predicted_prob_6717$d[i]<-"5"}
summary_clust_6717<-as.matrix(summary_clusters_6717)
write.table(summary_clust_6717, file="summary_clusters_6717.csv", sep = ",",qmethod = "double")       

predicted_clusters_6718<-predict(xyMclust6718,df_gauss_6718)
predicted_clus_6718<-data.frame(predicted_clusters_6718)
predicted_clust_6718<-as.matrix(predicted_clus_6718)
summary_clusters_6718<-summary(xyMclust6718, parameters=TRUE)
predicted_clus_6718<- within(predicted_clus_6718, Z1prob<- paste(pnorm(predicted_clus_6718$z.1), sep=''))
predicted_clus_6718 <- within(predicted_clus_6718, Z2prob <- paste(pnorm(predicted_clus_6718$z.2), sep=''))
predicted_clus_6718 <- within(predicted_clus_6718, Z3prob <- paste(pnorm(predicted_clus_6718$z.3), sep=''))
nuc_predicted_prob_6718<-cbind(gauss_6718, predicted_clus_6718)
nuc_predicted_prob_6718<- within(nuc_predicted_prob_6718, log_green <- paste(log(nuc_predicted_prob_6718$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6718 <- within(nuc_predicted_prob_6718, log_red <- paste(log(nuc_predicted_prob_6718$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6718$d<-nuc_predicted_prob_6718$classification
nuc_predicted_prob_6718[nuc_predicted_prob_6718$d=="1",64 ]<-"G1S"
nuc_predicted_prob_6718[nuc_predicted_prob_6718$d=="2",64 ]<-"G1"
nuc_predicted_prob_6718[nuc_predicted_prob_6718$d=="3",64 ]<-"SG2"
nuc_predicted_prob_6718$d->nuc_predicted_prob_6718$d1
nuc_predicted_prob_6718$a<-nuc_predicted_prob_6718$z.1>0.9 | nuc_predicted_prob_6718$z.2>0.9 | nuc_predicted_prob_6718$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6718) ){
  if(nuc_predicted_prob_6718$a[i]=="FALSE")
    nuc_predicted_prob_6718$d[i]<-"5"}
summary_clust_6718<-as.matrix(summary_clusters_6718)
write.table(summary_clust_6718, file="summary_clusters_6718.csv", sep = ",",qmethod = "double")  

predicted_clusters_6719<-predict(xyMclust6719,df_gauss_6719)
predicted_clus_6719<-data.frame(predicted_clusters_6719)
predicted_clust_6719<-as.matrix(predicted_clus_6719)
summary_clusters_6719<-summary(xyMclust6719, parameters=TRUE)
predicted_clus_6719<- within(predicted_clus_6719, Z1prob<- paste(pnorm(predicted_clus_6719$z.1), sep=''))
predicted_clus_6719 <- within(predicted_clus_6719, Z2prob <- paste(pnorm(predicted_clus_6719$z.2), sep=''))
predicted_clus_6719 <- within(predicted_clus_6719, Z3prob <- paste(pnorm(predicted_clus_6719$z.3), sep=''))
nuc_predicted_prob_6719<-cbind(gauss_6719, predicted_clus_6719)
nuc_predicted_prob_6719<- within(nuc_predicted_prob_6719, log_green <- paste(log(nuc_predicted_prob_6719$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6719 <- within(nuc_predicted_prob_6719, log_red <- paste(log(nuc_predicted_prob_6719$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6719$d<-nuc_predicted_prob_6719$classification
nuc_predicted_prob_6719[nuc_predicted_prob_6719$d=="1",64 ]<-"G1"
nuc_predicted_prob_6719[nuc_predicted_prob_6719$d=="2",64 ]<-"G1S"
nuc_predicted_prob_6719[nuc_predicted_prob_6719$d=="3",64 ]<-"SG2"
nuc_predicted_prob_6719$d->nuc_predicted_prob_6719$d1
nuc_predicted_prob_6719$a<-nuc_predicted_prob_6719$z.1>0.9 | nuc_predicted_prob_6719$z.2>0.9 | nuc_predicted_prob_6719$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6719) ){
  if(nuc_predicted_prob_6719$a[i]=="FALSE")
    nuc_predicted_prob_6719$d[i]<-"5"}
summary_clust_6719<-as.matrix(summary_clusters_6719)
write.table(summary_clust_6719, file="summary_clusters_6719.csv", sep = ",",qmethod = "double")  


predicted_clusters_6720<-predict(xyMclust6720,df_gauss_6720)
predicted_clus_6720<-data.frame(predicted_clusters_6720)
predicted_clust_6720<-as.matrix(predicted_clus_6720)
summary_clusters_6720<-summary(xyMclust6720, parameters=TRUE)
predicted_clus_6720<- within(predicted_clus_6720, Z1prob<- paste(pnorm(predicted_clus_6720$z.1), sep=''))
predicted_clus_6720 <- within(predicted_clus_6720, Z2prob <- paste(pnorm(predicted_clus_6720$z.2), sep=''))
predicted_clus_6720 <- within(predicted_clus_6720, Z3prob <- paste(pnorm(predicted_clus_6720$z.3), sep=''))
nuc_predicted_prob_6720<-cbind(gauss_6720, predicted_clus_6720)
nuc_predicted_prob_6720<- within(nuc_predicted_prob_6720, log_green <- paste(log(nuc_predicted_prob_6720$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6720 <- within(nuc_predicted_prob_6720, log_red <- paste(log(nuc_predicted_prob_6720$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6720$d<-nuc_predicted_prob_6720$classification
nuc_predicted_prob_6720[nuc_predicted_prob_6720$d=="1",64 ]<-"G1S"
nuc_predicted_prob_6720[nuc_predicted_prob_6720$d=="2",64 ]<-"G1"
nuc_predicted_prob_6720[nuc_predicted_prob_6720$d=="3",64 ]<-"SG2"
nuc_predicted_prob_6720$d->nuc_predicted_prob_6720$d1
nuc_predicted_prob_6720$a<-nuc_predicted_prob_6720$z.1>0.9 | nuc_predicted_prob_6720$z.2>0.9 | nuc_predicted_prob_6720$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6720) ){
  if(nuc_predicted_prob_6720$a[i]=="FALSE")
    nuc_predicted_prob_6720$d[i]<-"5"}
summary_clust_6720<-as.matrix(summary_clusters_6720)
write.table(summary_clust_6720, file="summary_clusters_6720.csv", sep = ",",qmethod = "double")  



predicted_clusters_6721<-predict(xyMclust6721,df_gauss_6721)
predicted_clus_6721<-data.frame(predicted_clusters_6721)
predicted_clust_6721<-as.matrix(predicted_clus_6721)
summary_clusters_6721<-summary(xyMclust6721, parameters=TRUE)
predicted_clus_6721<- within(predicted_clus_6721, Z1prob<- paste(pnorm(predicted_clus_6721$z.1), sep=''))
predicted_clus_6721 <- within(predicted_clus_6721, Z2prob <- paste(pnorm(predicted_clus_6721$z.2), sep=''))
predicted_clus_6721 <- within(predicted_clus_6721, Z3prob <- paste(pnorm(predicted_clus_6721$z.3), sep=''))
nuc_predicted_prob_6721<-cbind(gauss_6721, predicted_clus_6721)
nuc_predicted_prob_6721<- within(nuc_predicted_prob_6721, log_green <- paste(log(nuc_predicted_prob_6721$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6721 <- within(nuc_predicted_prob_6721, log_red <- paste(log(nuc_predicted_prob_6721$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6721$d<-nuc_predicted_prob_6721$classification
nuc_predicted_prob_6721[nuc_predicted_prob_6721$d=="1",64 ]<-"G1"
nuc_predicted_prob_6721[nuc_predicted_prob_6721$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6721[nuc_predicted_prob_6721$d=="3",64 ]<-"G1S"
nuc_predicted_prob_6721$d->nuc_predicted_prob_6721$d1
nuc_predicted_prob_6721$a<-nuc_predicted_prob_6721$z.1>0.9 | nuc_predicted_prob_6721$z.2>0.9 | nuc_predicted_prob_6721$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6721) ){
  if(nuc_predicted_prob_6721$a[i]=="FALSE")
    nuc_predicted_prob_6721$d[i]<-"5"}
summary_clust_6721<-as.matrix(summary_clusters_6721)
write.table(summary_clust_6721, file="summary_clusters_6721.csv", sep = ",",qmethod = "double")  


predicted_clusters_6722<-predict(xyMclust6722,df_gauss_6722)
predicted_clus_6722<-data.frame(predicted_clusters_6722)
predicted_clust_6722<-as.matrix(predicted_clus_6722)
summary_clusters_6722<-summary(xyMclust6722, parameters=TRUE)
predicted_clus_6722<- within(predicted_clus_6722, Z1prob<- paste(pnorm(predicted_clus_6722$z.1), sep=''))
predicted_clus_6722 <- within(predicted_clus_6722, Z2prob <- paste(pnorm(predicted_clus_6722$z.2), sep=''))
predicted_clus_6722 <- within(predicted_clus_6722, Z3prob <- paste(pnorm(predicted_clus_6722$z.3), sep=''))
nuc_predicted_prob_6722<-cbind(gauss_6722, predicted_clus_6722)
nuc_predicted_prob_6722<- within(nuc_predicted_prob_6722, log_green <- paste(log(nuc_predicted_prob_6722$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6722 <- within(nuc_predicted_prob_6722, log_red <- paste(log(nuc_predicted_prob_6722$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6722$d<-nuc_predicted_prob_6722$classification
nuc_predicted_prob_6722[nuc_predicted_prob_6722$d=="1",64 ]<-"G1"
nuc_predicted_prob_6722[nuc_predicted_prob_6722$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6722[nuc_predicted_prob_6722$d=="3",64 ]<-"G1S"
nuc_predicted_prob_6722$d->nuc_predicted_prob_6722$d1
nuc_predicted_prob_6722$a<-nuc_predicted_prob_6722$z.1>0.9 | nuc_predicted_prob_6722$z.2>0.9 | nuc_predicted_prob_6722$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6722) ){
  if(nuc_predicted_prob_6722$a[i]=="FALSE")
    nuc_predicted_prob_6722$d[i]<-"5"}
summary_clust_6722<-as.matrix(summary_clusters_6722)
write.table(summary_clust_6722, file="summary_clusters_6722.csv", sep = ",",qmethod = "double")  



predicted_clusters_6724<-predict(xyMclust6724,df_gauss_6724)
predicted_clus_6724<-data.frame(predicted_clusters_6724)
predicted_clust_6724<-as.matrix(predicted_clus_6724)
summary_clusters_6724<-summary(xyMclust6724, parameters=TRUE)
predicted_clus_6724<- within(predicted_clus_6724, Z1prob<- paste(pnorm(predicted_clus_6724$z.1), sep=''))
predicted_clus_6724 <- within(predicted_clus_6724, Z2prob <- paste(pnorm(predicted_clus_6724$z.2), sep=''))
predicted_clus_6724 <- within(predicted_clus_6724, Z3prob <- paste(pnorm(predicted_clus_6724$z.3), sep=''))
nuc_predicted_prob_6724<-cbind(gauss_6724, predicted_clus_6724)
nuc_predicted_prob_6724<- within(nuc_predicted_prob_6724, log_green <- paste(log(nuc_predicted_prob_6724$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6724 <- within(nuc_predicted_prob_6724, log_red <- paste(log(nuc_predicted_prob_6724$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6724$d<-nuc_predicted_prob_6724$classification
nuc_predicted_prob_6724[nuc_predicted_prob_6724$d=="1",64 ]<-"G1"
nuc_predicted_prob_6724[nuc_predicted_prob_6724$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6724[nuc_predicted_prob_6724$d=="3",64 ]<-"G1S"
nuc_predicted_prob_6724$d->nuc_predicted_prob_6724$d1
nuc_predicted_prob_6724$a<-nuc_predicted_prob_6724$z.1>0.9 | nuc_predicted_prob_6724$z.2>0.9 | nuc_predicted_prob_6724$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6724) ){
  if(nuc_predicted_prob_6724$a[i]=="FALSE")
    nuc_predicted_prob_6724$d[i]<-"5"}
summary_clust_6724<-as.matrix(summary_clusters_6724)
write.table(summary_clust_6724, file="summary_clusters_6724.csv", sep = ",",qmethod = "double")  



predicted_clusters_6725<-predict(xyMclust6725,df_gauss_6725)
predicted_clus_6725<-data.frame(predicted_clusters_6725)
predicted_clust_6725<-as.matrix(predicted_clus_6725)
summary_clusters_6725<-summary(xyMclust6725, parameters=TRUE)
predicted_clus_6725<- within(predicted_clus_6725, Z1prob<- paste(pnorm(predicted_clus_6725$z.1), sep=''))
predicted_clus_6725 <- within(predicted_clus_6725, Z2prob <- paste(pnorm(predicted_clus_6725$z.2), sep=''))
predicted_clus_6725 <- within(predicted_clus_6725, Z3prob <- paste(pnorm(predicted_clus_6725$z.3), sep=''))
nuc_predicted_prob_6725<-cbind(gauss_6725, predicted_clus_6725)
nuc_predicted_prob_6725<- within(nuc_predicted_prob_6725, log_green <- paste(log(nuc_predicted_prob_6725$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6725 <- within(nuc_predicted_prob_6725, log_red <- paste(log(nuc_predicted_prob_6725$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6725$d<-nuc_predicted_prob_6725$classification
nuc_predicted_prob_6725[nuc_predicted_prob_6725$d=="1",64 ]<-"G1"
nuc_predicted_prob_6725[nuc_predicted_prob_6725$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6725[nuc_predicted_prob_6725$d=="3",64 ]<-"G1S"
nuc_predicted_prob_6725$d->nuc_predicted_prob_6725$d1
nuc_predicted_prob_6725$a<-nuc_predicted_prob_6725$z.1>0.9 | nuc_predicted_prob_6725$z.2>0.9 | nuc_predicted_prob_6725$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6725) ){
  if(nuc_predicted_prob_6725$a[i]=="FALSE")
    nuc_predicted_prob_6725$d[i]<-"5"}
summary_clust_6725<-as.matrix(summary_clusters_6725)
write.table(summary_clust_6725, file="summary_clusters_6725.csv", sep = ",",qmethod = "double")  



predicted_clusters_6731<-predict(xyMclust6731,df_gauss_6731)
predicted_clus_6731<-data.frame(predicted_clusters_6731)
predicted_clust_6731<-as.matrix(predicted_clus_6731)
summary_clusters_6731<-summary(xyMclust6731, parameters=TRUE)
predicted_clus_6731<- within(predicted_clus_6731, Z1prob<- paste(pnorm(predicted_clus_6731$z.1), sep=''))
predicted_clus_6731 <- within(predicted_clus_6731, Z2prob <- paste(pnorm(predicted_clus_6731$z.2), sep=''))
predicted_clus_6731 <- within(predicted_clus_6731, Z3prob <- paste(pnorm(predicted_clus_6731$z.3), sep=''))
nuc_predicted_prob_6731<-cbind(gauss_6731, predicted_clus_6731)
nuc_predicted_prob_6731<- within(nuc_predicted_prob_6731, log_green <- paste(log(nuc_predicted_prob_6731$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6731 <- within(nuc_predicted_prob_6731, log_red <- paste(log(nuc_predicted_prob_6731$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6731$d<-nuc_predicted_prob_6731$classification
nuc_predicted_prob_6731[nuc_predicted_prob_6731$d=="1",64 ]<-"SG2"
nuc_predicted_prob_6731[nuc_predicted_prob_6731$d=="2",64 ]<-"G1S"
nuc_predicted_prob_6731[nuc_predicted_prob_6731$d=="3",64 ]<-"G1"
nuc_predicted_prob_6731$d->nuc_predicted_prob_6731$d1
nuc_predicted_prob_6731$a<-nuc_predicted_prob_6731$z.1>0.9 | nuc_predicted_prob_6731$z.2>0.9 | nuc_predicted_prob_6731$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6731) ){
  if(nuc_predicted_prob_6731$a[i]=="FALSE")
    nuc_predicted_prob_6731$d[i]<-"5"}
summary_clust_6731<-as.matrix(summary_clusters_6731)
write.table(summary_clust_6731, file="summary_clusters_6731.csv", sep = ",",qmethod = "double")  



predicted_clusters_6734<-predict(xyMclust6734,df_gauss_6734)
predicted_clus_6734<-data.frame(predicted_clusters_6734)
predicted_clust_6734<-as.matrix(predicted_clus_6734)
summary_clusters_6734<-summary(xyMclust6734, parameters=TRUE)
predicted_clus_6734<- within(predicted_clus_6734, Z1prob<- paste(pnorm(predicted_clus_6734$z.1), sep=''))
predicted_clus_6734 <- within(predicted_clus_6734, Z2prob <- paste(pnorm(predicted_clus_6734$z.2), sep=''))
predicted_clus_6734 <- within(predicted_clus_6734, Z3prob <- paste(pnorm(predicted_clus_6734$z.3), sep=''))
nuc_predicted_prob_6734<-cbind(gauss_6734, predicted_clus_6734)
nuc_predicted_prob_6734<- within(nuc_predicted_prob_6734, log_green <- paste(log(nuc_predicted_prob_6734$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6734 <- within(nuc_predicted_prob_6734, log_red <- paste(log(nuc_predicted_prob_6734$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6734$d<-nuc_predicted_prob_6734$classification
nuc_predicted_prob_6734[nuc_predicted_prob_6734$d=="1",64 ]<-"SG2"
nuc_predicted_prob_6734[nuc_predicted_prob_6734$d=="2",64 ]<-"G1S"
nuc_predicted_prob_6734[nuc_predicted_prob_6734$d=="3",64 ]<-"G1"
nuc_predicted_prob_6734$d->nuc_predicted_prob_6734$d1
nuc_predicted_prob_6734$a<-nuc_predicted_prob_6734$z.1>0.9 | nuc_predicted_prob_6734$z.2>0.9 | nuc_predicted_prob_6734$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6734) ){
  if(nuc_predicted_prob_6734$a[i]=="FALSE")
    nuc_predicted_prob_6734$d[i]<-"5"}
summary_clust_6734<-as.matrix(summary_clusters_6734)
write.table(summary_clust_6734, file="summary_clusters_6734.csv", sep = ",",qmethod = "double")  


predicted_clusters_6735<-predict(xyMclust6735,df_gauss_6735)
predicted_clus_6735<-data.frame(predicted_clusters_6735)
predicted_clust_6735<-as.matrix(predicted_clus_6735)
summary_clusters_6735<-summary(xyMclust6735, parameters=TRUE)
predicted_clus_6735<- within(predicted_clus_6735, Z1prob<- paste(pnorm(predicted_clus_6735$z.1), sep=''))
predicted_clus_6735 <- within(predicted_clus_6735, Z2prob <- paste(pnorm(predicted_clus_6735$z.2), sep=''))
predicted_clus_6735 <- within(predicted_clus_6735, Z3prob <- paste(pnorm(predicted_clus_6735$z.3), sep=''))
nuc_predicted_prob_6735<-cbind(gauss_6735, predicted_clus_6735)
nuc_predicted_prob_6735<- within(nuc_predicted_prob_6735, log_green <- paste(log(nuc_predicted_prob_6735$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6735 <- within(nuc_predicted_prob_6735, log_red <- paste(log(nuc_predicted_prob_6735$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6735$d<-nuc_predicted_prob_6735$classification
nuc_predicted_prob_6735[nuc_predicted_prob_6735$d=="1",64 ]<-"G1S"
nuc_predicted_prob_6735[nuc_predicted_prob_6735$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6735[nuc_predicted_prob_6735$d=="3",64 ]<-"G1"
nuc_predicted_prob_6735$d->nuc_predicted_prob_6735$d1
nuc_predicted_prob_6735$a<-nuc_predicted_prob_6735$z.1>0.9 | nuc_predicted_prob_6735$z.2>0.9 | nuc_predicted_prob_6735$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6735) ){
  if(nuc_predicted_prob_6735$a[i]=="FALSE")
    nuc_predicted_prob_6735$d[i]<-"5"}
summary_clust_6735<-as.matrix(summary_clusters_6735)
write.table(summary_clust_6735, file="summary_clusters_6735.csv", sep = ",",qmethod = "double")  




predicted_clusters_6736<-predict(xyMclust6736,df_gauss_6736)
predicted_clus_6736<-data.frame(predicted_clusters_6736)
predicted_clust_6736<-as.matrix(predicted_clus_6736)
summary_clusters_6736<-summary(xyMclust6736, parameters=TRUE)
predicted_clus_6736<- within(predicted_clus_6736, Z1prob<- paste(pnorm(predicted_clus_6736$z.1), sep=''))
predicted_clus_6736 <- within(predicted_clus_6736, Z2prob <- paste(pnorm(predicted_clus_6736$z.2), sep=''))
predicted_clus_6736 <- within(predicted_clus_6736, Z3prob <- paste(pnorm(predicted_clus_6736$z.3), sep=''))
nuc_predicted_prob_6736<-cbind(gauss_6736, predicted_clus_6736)
nuc_predicted_prob_6736<- within(nuc_predicted_prob_6736, log_green <- paste(log(nuc_predicted_prob_6736$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6736 <- within(nuc_predicted_prob_6736, log_red <- paste(log(nuc_predicted_prob_6736$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6736$d<-nuc_predicted_prob_6736$classification
nuc_predicted_prob_6736[nuc_predicted_prob_6736$d=="1",64 ]<-"G1"
nuc_predicted_prob_6736[nuc_predicted_prob_6736$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6736[nuc_predicted_prob_6736$d=="3",64]<-"G1S"
nuc_predicted_prob_6736$d->nuc_predicted_prob_6736$d1
nuc_predicted_prob_6736$a<-nuc_predicted_prob_6736$z.1>0.9 | nuc_predicted_prob_6736$z.2>0.9 | nuc_predicted_prob_6736$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6736) ){
  if(nuc_predicted_prob_6736$a[i]=="FALSE")
    nuc_predicted_prob_6736$d[i]<-"5"}
summary_clust_6736<-as.matrix(summary_clusters_6736)
write.table(summary_clust_6736, file="summary_clusters_6736.csv", sep = ",",qmethod = "double")  


predicted_clusters_6745<-predict(xyMclust6745,df_gauss_6745)
predicted_clus_6745<-data.frame(predicted_clusters_6745)
predicted_clust_6745<-as.matrix(predicted_clus_6745)
summary_clusters_6745<-summary(xyMclust6745, parameters=TRUE)
predicted_clus_6745<- within(predicted_clus_6745, Z1prob<- paste(pnorm(predicted_clus_6745$z.1), sep=''))
predicted_clus_6745 <- within(predicted_clus_6745, Z2prob <- paste(pnorm(predicted_clus_6745$z.2), sep=''))
predicted_clus_6745 <- within(predicted_clus_6745, Z3prob <- paste(pnorm(predicted_clus_6745$z.3), sep=''))
nuc_predicted_prob_6745<-cbind(gauss_6745, predicted_clus_6745)
nuc_predicted_prob_6745<- within(nuc_predicted_prob_6745, log_green <- paste(log(nuc_predicted_prob_6745$Intensity_MeanIntensity_CorrResizedGreenFUCCI), sep=''))
nuc_predicted_prob_6745 <- within(nuc_predicted_prob_6745, log_red <- paste(log(nuc_predicted_prob_6745$Intensity_MeanIntensity_CorrResizedRedFUCCI), sep=''))
#sub<-subset(nuc_predicted_prob_phases, nuc_predicted_prob_phases$z.1>0.9 | nuc_predicted_prob_phases$z.2>0.9 | nuc_predicted_prob_phases$z.3>0.9)
nuc_predicted_prob_6745$d<-nuc_predicted_prob_6745$classification
nuc_predicted_prob_6745[nuc_predicted_prob_6745$d=="1",64 ]<-"G1S"
nuc_predicted_prob_6745[nuc_predicted_prob_6745$d=="2",64 ]<-"SG2"
nuc_predicted_prob_6745[nuc_predicted_prob_6745$d=="3",64 ]<-"G1"
nuc_predicted_prob_6745$d->nuc_predicted_prob_6745$d1
nuc_predicted_prob_6745$a<-nuc_predicted_prob_6745$z.1>0.9 | nuc_predicted_prob_6745$z.2>0.9 | nuc_predicted_prob_6745$z.3>0.9
for (i in 1:nrow(nuc_predicted_prob_6745) ){
  if(nuc_predicted_prob_6745$a[i]=="FALSE")
    nuc_predicted_prob_6745$d[i]<-"5"}
summary_clust_6745<-as.matrix(summary_clusters_6745)
write.table(summary_clust_6745, file="summary_clusters_6745.csv", sep = ",",qmethod = "double")  






#plot(xyMclust1,col=c("red","green","orange"),pch=c(46,15,14), cex=0.5,xlab="",ylab="")
#plot(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)", col=factor(nuc_cyto_all_excl_outoffocus$CC),pch=16, cex=0.4)
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6717",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6717, parameters=xyMclust6717$parameters,
             z=xyMclust6717$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6718",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6718, parameters=xyMclust6718$parameters,
             z=xyMclust6718$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6719",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6719, parameters=xyMclust6719$parameters,
             z=xyMclust6719$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6720",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6720, parameters=xyMclust6720$parameters,
             z=xyMclust6720$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6721",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6721, parameters=xyMclust6721$parameters,
             z=xyMclust6721$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6722",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6722, parameters=xyMclust6722$parameters,
             z=xyMclust6722$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6724",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6724, parameters=xyMclust6724$parameters,
             z=xyMclust6724$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6725",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6725, parameters=xyMclust6725$parameters,
             z=xyMclust6725$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6731",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6731, parameters=xyMclust6731$parameters,
             z=xyMclust6731$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6734",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6734, parameters=xyMclust6734$parameters,
             z=xyMclust6734$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6735",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6735 ,parameters=xyMclust6735$parameters,
             z=xyMclust6735$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6736",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6736, parameters=xyMclust6736$parameters,
             z=xyMclust6736$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()
pdf(paste(workingdir,paste( "log_Fucci_gaussian_clustering_6745",".pdf", sep=""),sep="/"))
par(pty="s")
mclust2Dplot(df_gauss_6745, parameters=xyMclust6745$parameters,
             z=xyMclust6745$z, what = "classification", identify = TRUE,col=c("green","red","orange"),CEX=0.2,symbols = c(16,16,16),xlab="Log(Mean Intensity Green)",ylab="Log(Mean Intensity Red)",xlim=c(-6,-1))
dev.off()

nuc_predicted_prob_phases<-rbind(nuc_predicted_prob_6717,nuc_predicted_prob_6718,nuc_predicted_prob_6719,nuc_predicted_prob_6720,nuc_predicted_prob_6721,nuc_predicted_prob_6722,nuc_predicted_prob_6724, nuc_predicted_prob_6725,nuc_predicted_prob_6731,nuc_predicted_prob_6734,nuc_predicted_prob_6735,nuc_predicted_prob_6736,nuc_predicted_prob_6745)

write.table(nuc_predicted_prob_phases, file = "nuc_predicted_prob_phases.csv", sep = ",",qmethod = "double")
data_nuc_cyto_merged<-nuc_predicted_prob_phases
######remove the ones that are on border of the clusters
data<-subset(data_nuc_cyto_merged, data_nuc_cyto_merged$d=="G1" | data_nuc_cyto_merged$d=="G1S" | data_nuc_cyto_merged$d=="SG2")



workingdir<-setwd("/Users/diana.telessemian/Desktop/gaussian_per_plate/plots")

well_plate<- unique(data$well_plate)
for (j in well_plate) {
  
  datanuccyto<- data[grep(j,data$well_plate),]
  
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  #datanuccyto<- data[data$well_plate=="G08_55375988",]
  #datanuccyto<- datanuccyto[datanuccyto$ImageNumber!="80",]
  G1<-datanuccyto[datanuccyto$d=="G1", ]
  G1S<-datanuccyto[datanuccyto$d=="G1S", ]
  SG2<-datanuccyto[datanuccyto$d=="SG2", ]
  
  G1_ab_int_nuc <- (G1$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb))
  G1S_ab_int_nuc <- (G1S$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb))
  SG2_ab_int_nuc <- (SG2$Intensity_MeanIntensity_ResizedAb/max(datanuccyto$Intensity_MeanIntensity_ResizedAb))
  ab_intenity_phases_nuc<-list(G1_ab_int_nuc ,G1S_ab_int_nuc ,SG2_ab_int_nuc)
  
  G1_ab_int_cyto <- (G1$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto))
  G1S_ab_int_cyto <- (G1S$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto))
  SG2_ab_int_cyto <- (SG2$Mean_ab_Cyto/max(datanuccyto$Mean_ab_Cyto))
  ab_intenity_phases_cyto<-list(G1_ab_int_cyto ,G1S_ab_int_cyto ,SG2_ab_int_cyto)
  
  
  G1_ab_int_cell <- (G1$Mean_ab_cell/max(datanuccyto$Mean_ab_cell))
  G1S_ab_int_cell <- (G1S$Mean_ab_cell/max(datanuccyto$Mean_ab_cell))
  SG2_ab_int_cell <- (SG2$Mean_ab_cell/max(datanuccyto$Mean_ab_cell))
  ab_intenity_phases_cell<-list(G1_ab_int_cell ,G1S_ab_int_cell ,SG2_ab_int_cell)
  
  require(PMCMR)
  datanuccyto$d<-as.factor(datanuccyto$d)
  kruskal_nuc<-kruskal.test(datanuccyto$Intensity_MeanIntensity_ResizedAb~ datanuccyto$d, data = datanuccyto)
  pvalue_eachclass_nuc<-posthoc.kruskal.nemenyi.test(x=datanuccyto$Intensity_MeanIntensity_ResizedAb, g=datanuccyto$d, method="Tukey") 
  kruskal_cyto<-kruskal.test(datanuccyto$Mean_ab_Cyto~ datanuccyto$d, data = datanuccyto)
  pvalue_eachclass_cyto<-posthoc.kruskal.nemenyi.test(x=datanuccyto$Mean_ab_Cyto, g=datanuccyto$d, method="Tukey")
  kruskal_cell<-kruskal.test(datanuccyto$Mean_ab_cell~ datanuccyto$d, data = datanuccyto)
  pvalue_eachclass_cell<-posthoc.kruskal.nemenyi.test(x=datanuccyto$Mean_ab_cell, g=datanuccyto$d, method="Tukey")
  
  pdf(paste(workingdir,paste( "boxplotnuc_",j,".pdf", sep=""),sep="/"))
  #pdf("boxplotnuc.pdf")
  par(pty="s")
  boxplotnuc <- boxplot(ab_intenity_phases_nuc , main="",ylab="Mean Intensity",las=1, xlab="", ylim=c(0,max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)), border=c(),names=c(),cex.axis=0.9,xaxt="n.axis",frame.plot=FALSE,cex.lab=1.5)
  segments(1,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),1.8)
  segments(1,0.88*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),1,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  segments(1.8,0.88*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),1.8,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  segments(2.2,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),3)
  segments(2.2,0.88*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),2.2,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  segments(3,0.88*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),3,0.9*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  segments(1,0.98*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),3)
  segments(3,0.96*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),3,0.98*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  segments(1,0.96*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)),1,0.98*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)))
  text(1.4,0.92*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)), signif(pvalue_eachclass_nuc$p.value[1,1], digits=2), col = "black", cex=1.2)
  text(2,1*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)), signif(pvalue_eachclass_nuc$p.value[2,1], digits=2), col = "black", cex=1.2)
  text(2.5,0.92*(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)+(max(datanuccyto$Intensity_MeanIntensity_ResizedAb/datanuccyto$Intensity_MeanIntensity_ResizedAb)/4)), signif(pvalue_eachclass_nuc$p.value[2,2], digits=2), col = "black", cex=1.2)
  text(2.1,1, signif(kruskal_nuc$p.value,digits=2), col = "black", cex=1.2)
  text(1.8,1, "P= " , col = "black", cex=1.2)
  text(1,0, "G1", col = "black", cex=1.2)
  text(2,0, "G1S" , col = "black", cex=1.2)
  text(3,0, "SG2", col = "black", cex=1.2)
  dev.off()
  
  pdf(paste(workingdir,paste( "boxplotcyto_",j,".pdf", sep=""),sep="/"))
  par(pty="s")
  boxplotcyto <- boxplot(ab_intenity_phases_cyto , main="",ylab="Mean Intensity",las=1, ylim=c(0,max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)), border=c(),names=c(),cex.axis=0.9,xaxt="n.axis",frame.plot=FALSE,cex.lab=1.5)
  segments(1,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),1.8)
  segments(1,0.88*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),1,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  segments(1.8,0.88*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),1.8,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  segments(2.2,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),3)
  segments(2.2,0.88*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),2.2,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  segments(3,0.88*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),3,0.9*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  segments(1,0.98*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),3)
  segments(3,0.96*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),3,0.98*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  segments(1,0.96*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),1,0.98*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)))
  text(1.4,0.92*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)), signif(pvalue_eachclass_cyto$p.value[1,1], digits=2), col = "black", cex=1.2)
  text(2,1*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),  signif(pvalue_eachclass_cyto$p.value[2,1], digits=2), col = "black", cex=1.2)
  text(2.5,0.92*(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)+(max(datanuccyto$Mean_ab_Cyto/datanuccyto$Mean_ab_Cyto)/4)),  signif(pvalue_eachclass_cyto$p.value[2,2], digits=2), col = "black", cex=1.2)
  text(2.1,1, signif(kruskal_cyto$p.value,digits=2), col = "black", cex=1.2)
  text(1.8,1, "P= " , col = "black", cex=1.2)
  text(1,0, "G1", col = "black", cex=1.2)
  text(2,0, "G1S" , col = "black", cex=1.2)
  text(3,0, "SG2", col = "black", cex=1.2)
  dev.off()
  
  pdf(paste(workingdir,paste( "boxplotcell_",j,".pdf", sep=""),sep="/"))
  par(pty="s")
  boxplotcell <- boxplot(ab_intenity_phases_cell , main="",ylab="Mean Intensity",las=1, xlab="", ylim=c(0,max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cel)/4)), border=c(),names=c(),cex.axis=0.9,xaxt="n.axis",frame.plot=FALSE, cex.lab=1.5)
  segments(1,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),1.8)
  segments(1,0.88*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),1,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  segments(1.8,0.88*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),1.8,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  segments(2.2,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),3)
  segments(2.2,0.88*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),2.2,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  segments(3,0.88*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),3,0.9*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  segments(1,0.98*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),3)
  segments(3,0.96*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),3,0.98*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  segments(1,0.96*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),1,0.98*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)))
  text(1.4,0.92*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)), signif(pvalue_eachclass_cell$p.value[1,1], digits=2), col = "black", cex=1.2)
  text(2,1*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),  signif(pvalue_eachclass_cell$p.value[2,1], digits=2), col = "black", cex=1.2)
  text(2.5,0.92*(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)+(max(datanuccyto$Mean_ab_cell/datanuccyto$Mean_ab_cell)/4)),  signif(pvalue_eachclass_cell$p.value[2,2], digits=2), col = "black", cex=1.2)
  text(2.1,1, signif(kruskal_cell$p.value,digits=2), col = "black", cex=1.2)
  text(1.8,1, "P= " , col = "black", cex=1.2)
  text(1,0, "G1", col = "black", cex=1.2)
  text(2,0, "G1S" , col = "black", cex=1.2)
  text(3,0, "SG2", col = "black", cex=1.2)
  dev.off()
}

require(ggplot2)
well_plate<- unique(data_nuc_cyto_merged$well_plate)

for (i in well_plate) {
  data_nuc_cyto_well<- data_nuc_cyto_merged[grep(i,data_nuc_cyto_merged$well_plate),]
  
  #datanuccyto<- data[data$well_plate=="E11_55345985",]
  
  #data_nuc_cyto_well<- data_nuc_cyto_merged[data_nuc_cyto_merged$well_plate=="H02_5540",]
  rbPal <- colorRampPalette(c('yellow','red','blue'))
  #rbPal<-colorRampPalette(topo.colors(10))
  data_nuc_cyto_well$Col <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Mean_ab_Cyto,breaks = 10))]
  cuts<-levels(cut(data_nuc_cyto_well$Mean_ab_Cyto,breaks = 10))
  cuts<-gsub(","," - ",cuts)
  cuts<-gsub("\\(","[",cuts)
  tiff(paste(workingdir,paste( "cc_plot_cyto_",i,".tiff", sep=""),sep="/"))
  par(pty="s")
  cc_plot_cyto<-hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
                       col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main ="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)", col=data_nuc_cyto_well$Col,pch=16, cex=0.7)
  #legend("top",cuts,col=rbPal(10),pch=16)
  dev.off()
  
  tiff(paste(workingdir,paste( "cc_plot_nuc_",i,".tiff", sep=""),sep="/"))
  cc_plot_nuc<-data_nuc_cyto_well$Colnuc <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb,breaks = 10))]
  cuts_nuc<-levels(cut(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb,breaks = 10))
  cuts_nuc<-gsub(","," - ",cuts_nuc)
  cuts_nuc<-gsub("\\(","[",cuts_nuc)
  hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
         col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)", col=data_nuc_cyto_well$Colnuc,pch=16, cex=0.7)
  #legend("top",cuts_nuc,col=rbPal(10),pch=16, cex=1)
  dev.off()
  
  
  tiff(paste(workingdir,paste( "cc_plot_cell_",i,".tiff", sep=""),sep="/"))
  cc_plot_cell<-data_nuc_cyto_well$Colcell <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Mean_ab_cell,breaks = 10))]
  cuts_cell<-levels(cut(data_nuc_cyto_well$Mean_ab_cell,breaks = 10))
  cuts_cell<-gsub(","," - ",cuts_cell)
  cuts_cell<-gsub("\\(","[",cuts_cell)
  hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
         col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean_Intensity_Red)",xlab="Log(Mean_Green)", col=data_nuc_cyto_well$Colcell,pch=16, cex=0.7)
  #legend("top",cuts_nuc,col=rbPal(10),pch=16, cex=1)
  dev.off()
  
}
workingdir<-setwd("/Users/diana.telessemian/Desktop/kruskal_fv/norm_cc_plots")
library(ggplot2)
require(gplots)
well_plate<- unique(data_nuc_cyto_merged$well_plate)

for (i in well_plate) {
  data_nuc_cyto_well<- data_nuc_cyto_merged[grep(i,data_nuc_cyto_merged$well_plate),]
  
  #datanuccyto<- data[data$well_plate=="H02_55405991",]
  
  #data_nuc_cyto_well<- data_nuc_cyto_merged[data_nuc_cyto_merged$well_plate=="H02_55405991",]
  rbPal <- colorRampPalette(c('yellow','red','blue'))
  #rbPal<-colorRampPalette(topo.colors(10))
  data_nuc_cyto_well$Col <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Mean_ab_Cyto/max(data_nuc_cyto_well$Mean_ab_Cyto),breaks = 10))]
  cuts<-levels(cut(data_nuc_cyto_well$Mean_ab_Cyto/max(data_nuc_cyto_well$Mean_ab_Cyto),breaks = 10))
  cuts<-gsub(","," - ",cuts)
  cuts<-gsub("\\(","[",cuts)
  tiff(paste(workingdir,paste( "cc_plot_cyto_",i,".tiff", sep=""),sep="/"))
  par(pty="s")
  cc_plot_cyto<-hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
                       col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main ="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)", col=data_nuc_cyto_well$Col,pch=16, cex=0.7)
  #legend("top",cuts,col=rbPal(10),pch=16)
  dev.off()
  
  tiff(paste(workingdir,paste( "cc_plot_nuc_",i,".tiff", sep=""),sep="/"))
  cc_plot_nuc<-data_nuc_cyto_well$Colnuc <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb/max(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb),breaks = 10))]
  cuts_nuc<-levels(cut(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb/max(data_nuc_cyto_well$Intensity_MeanIntensity_ResizedAb),breaks = 10))
  cuts_nuc<-gsub(","," - ",cuts_nuc)
  cuts_nuc<-gsub("\\(","[",cuts_nuc)
  hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
         col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)", col=data_nuc_cyto_well$Colnuc,pch=16, cex=0.7)
  #legend("top",cuts_nuc,col=rbPal(10),pch=16, cex=1)
  dev.off()
  
  
  tiff(paste(workingdir,paste( "cc_plot_cell_",i,".tiff", sep=""),sep="/"))
  cc_plot_cell<-data_nuc_cyto_well$Colcell <- rbPal(10)[as.numeric(cut(data_nuc_cyto_well$Mean_ab_cell/max(data_nuc_cyto_well$Mean_ab_cell),breaks = 10))]
  cuts_cell<-levels(cut(data_nuc_cyto_well$Mean_ab_cell/max(data_nuc_cyto_well$Mean_ab_cell),breaks = 10))
  cuts_cell<-gsub(","," - ",cuts_cell)
  cuts_cell<-gsub("\\(","[",cuts_cell)
  hist2d(log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(nuc_cyto_all_excl_outoffocus$Intensity_MeanIntensity_CorrResizedRedFUCCI),nbins=200, same.scale=FALSE, na.rm=TRUE, show=TRUE, FUN=base::length,ylab="Log(Mean Intensity Red)",xlab="Log(Mean Intensity Green)",
         col=c("white", gray.colors(12,start=0.7, end=1 )),cex.lab=1.5, frame.plot=FALSE, main="",xlim=c(-6,-1))
  points(log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedGreenFUCCI),log(data_nuc_cyto_well$Intensity_MeanIntensity_CorrResizedRedFUCCI),ylab="Log(Mean_Intensity_Red)",xlab="Log(Mean_Green)", col=data_nuc_cyto_well$Colcell,pch=16, cex=0.7)
  #legend("top",cuts_nuc,col=rbPal(10),pch=16, cex=1)
  dev.off()
  
}
















