################################################
#### Prt1 :Merge new files with old files ######
#### date1 = date of the old one ###############
#### date2 = new date ##########################
################################################
date1 = "12062020"
date2 = "19062020"
REF ="EPI_ISL_402125"
DIR = "~/Desktop/Genomics-COVID19/"

dataold <- read.table(paste(DIR,"summaryfile_",date1,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)
snpsold <- read.table(paste(DIR,"snps_summaryfile_",date1,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)
datanew <- read.table(paste(DIR,"fromrim/summaryfile_",date2,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)
snpsnew <- read.table(paste(DIR,"fromrim/snps_summaryfile_",date2,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)

library(lubridate)
getyearweek <- function(date_vect){
  year_week<-vector();
  datetemp<-unlist(strsplit(as.character(parse_date_time(date_vect, c("dmy", "ymd","y","ym")))," "))
  year_week<-(paste(year(datetemp),week(datetemp),sep="-"))
  return(year_week)
}
datanew$actualdate <- datanew$date
datanew$date <- getyearweek(datanew$date)
snpsnew$actualdate <- snpsnew$date
snpsnew$date <- getyearweek(snpsnew$date)

data <- rbind(dataold,datanew[-1,])
snps <- rbind(snpsold,snpsnew)

if (length(which(duplicated(data)))) {
  data <- data[-which(duplicated(data)),]
} 
if (length(which(duplicated(snps)))) {
  snps <- snps[-which(duplicated(snps)),]
}
## Saving new files
write.table(data,paste(DIR,"summaryfile_",date2,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)
write.table(snps,paste(DIR,"snps_summaryfile_",date2,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)
### Save annotation file
annot <- readLines(paste(DIR,"annotation/mutations_annotations_",date2,".vcf",sep=""))
annot <- annot[grep("NC_045512.2",annot,fixed = T)]
annot <- annot[-1]
pos <- unlist(lapply(annot,function(x){unlist(strsplit(x,"\t"))[2]}))
ref <- unlist(lapply(annot,function(x){unlist(strsplit(x,"\t"))[4]}))
qry <- unlist(lapply(annot,function(x){unlist(strsplit(x,"\t"))[5]}))
mut <- paste(pos,ref,">",qry,sep = "")
annotation <- unlist(lapply(annot,function(x){unlist(strsplit(unlist(strsplit(x,"\t"))[8],"\\|"))[2]}))
towrite <- data.frame(mutation=mut,ANN=annotation)
old <- read.table(paste(DIR,"annotation/mutations_annotations_parsed_",date1,".txt",sep=""),stringsAsFactors = F,header = T)
all <- rbind(old,towrite) 
all <- unique(all)
towrite=all
write.table(towrite,paste(DIR,"annotation/mutations_annotations_parsed_",date2,".txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")

tabmut <- table(snps$mut)
tabmut <- data.frame(tabmut)
tabmut$effect <- NA

get_snpEFFmutation <- function(x) {
  if (length(which(towrite$mutation==x))) {
    as.character(towrite$ANN[which(towrite$mutation==x)])[1]
  } else {
    "No annotation"
  }
}
tabmut$snpEFF <- unlist(lapply(as.character(tabmut$Var1),FUN = get_snpEFFmutation))
tabmut$effect<-NULL
tabmut <- tabmut[-which(is.na(tabmut$snpEFF)),]
write.table(tabmut[order(tabmut$Freq,decreasing = T),],paste(DIR,"annotation/mutations_annotations_parsed_withfreq",date2,".txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
date=date2
rm(dataold,datanew,snpsold,snpsnew,old,all,tabmut,towrite,annot,annotation,date1,date2)
###############################################
####### Part2: Create files for dashboard #####
###############################################
orfs <- read.table(paste(DIR,"orfs_covid19.txt",sep=""),stringsAsFactors = F,header = F)
br <- (sort(unique(c(0,as.vector(orfs$V1),as.vector(orfs$V2),30000))))
lab <- c("other","ORF1A","ORF1B","other","S","other","ORF3A","other","E","other","M","other","ORF6","other","ORF7A","ORF7B","other","ORF8","other","N","other","ORF10","other")
orfs$V4 <- orfs$V2-orfs$V1 + 1
s1 = 30000-sum(orfs$V4)
orfs<-rbind(orfs,c(1,1,"other",s1)) 
row.names(orfs) <- orfs$V3
heatmapfile <- NULL
for (country in unique(snps$country)) {
  temp = snps[which(snps$country==country),]
  res <- data.frame(Location=as.character(orfs$V3))
  row.names(res) <- res$Location
  
  if (length(unique(temp$id))<2) {
    id = unique(temp$id)
    res = as.data.frame(table(temp$orf[which(temp$id==id)]))
    res <- t(res)
    colnames(res) <- res[1,]
    if (length(colnames(res)) < 13) {
      restemp <- data.frame(Location=as.character(orfs$V3),Freq=rep(0,length(orfs$V3)))
      row.names(restemp) <- restemp$Location
      restemp[colnames(res),2] <- as.numeric(res[2,])
      res <- restemp
      res <- t(res)
      colnames(res) <- res[1,]
    } 
    res <- as.numeric(res[-1,orfs$V3])/as.numeric(orfs[orfs$V3,4])*1000
    res <- t(res)
    colnames(res) <- orfs$V3
    data2write <- data.frame(Country=country,ID=unique(temp$id),Date=unique(temp$date))
    data2write <- cbind(data2write,res)
  } else {
    data2write = data.frame(Country=character(0),ID=character(0),Date=character(0))
    for (k in 1:length(unique(temp$id))) {
      id = unique(temp$id)[k]
      tabtemp = as.data.frame(table(temp$orf[which(temp$id==id)]))
      colnames(tabtemp) = c("Location",id)
      res <- merge(res,tabtemp,all.x = T,all.y = F,by="Location")
      row.names(res) <- res$Location
      data2write <- rbind(data2write,data.frame(Country=country,ID=id,Date=temp$date[which(temp$id==id)][1]))
    }
    res <-t(res[orfs$V3,-c(1)]/as.numeric(orfs[orfs$V3,4])*1000)
    res[is.na(res)] <- 0
    data2write <- cbind(data2write,res)
  }
  heatmapfile <- rbind(heatmapfile,data2write)
}

write.csv(heatmapfile,paste(DIR,"files4dashboard/counts_per_genekbp_percountry_",date,"_",REF,".csv",sep=""),row.names = F,quote = F)

piechartfile <- orfs$V3
piechartfile <- c("Country",piechartfile)
for (country in unique(snps$country)) {
  temp = as.data.frame(table(snps$orf[which(snps$country==country)])/length(which(data$country==country)))
  row.names(temp) <- temp$Var1
  temp[orfs$V3,2] = temp[orfs$V3,2]/as.numeric(orfs$V4)*1000
  temp$Freq[which(is.na(temp$Freq))] <- 0
  temp$angle <- temp$Freq*360/sum(temp$Freq)
  towrite <- c(country,temp[orfs$V3,3])
  piechartfile <- rbind(piechartfile,towrite)
}  

colnames(piechartfile) <- piechartfile[1,]
piechartfile <- piechartfile[-1,]
write.csv(piechartfile,paste(DIR,"files4dashboard/averagecounts_per_genekbp_percountry_",date,"_",REF,".csv",sep=""),row.names = F,quote = F)

numbers <- as.data.frame(table(data$country))
numbers <- rbind(numbers,data.frame(Var1="Global",Freq=length(data$country)))
write.csv(numbers,paste(DIR,"files4dashboard/number_of_samples_",date,"_",REF,".csv",sep=""),row.names = F,quote = F)

mut <- read.table(paste(DIR,"annotation/mutations_annotations_parsed_withfreq",date,".txt",sep=""),header = T,sep="\t",stringsAsFactors = F)

geteffect <- function(x) {
  if (length(which(mut$Var1==x))) {res = mut$snpEFF[which(mut$Var1==x)]} else{res = "No annotation"}
  res
}

snps <- snps[-which(snps$Qry.sub=="."),]

snps$eff <- unlist(lapply(snps$mut,geteffect))
snps$eff[which(snps$eff=="frameshift_variant&stop_gained")] <- "Frameshift"
snps$eff[which(snps$eff=="frameshift_variant")] <- "Frameshift"
snps$eff[which(snps$eff=="missense_variant&splice_region_variant")] <- "Missense"
snps$eff[which(snps$eff=="downstream_gene_variant")] <- "Downstream gene"
snps$eff[which(snps$eff=="upstream_gene_variant")] <- "Upstream gene"
snps$eff[which(snps$eff=="missense_variant")] <- "Missense"
snps$eff[which(snps$eff=="splice_region_variant&synonymous_variant")] <- "Splice_region"
snps$eff[which(snps$eff=="splice_region_variant")] <- "Splice_region"
snps$eff[which(snps$eff=="splice_region_variant&stop_retained_variant")] <- "Splice_region"
snps$eff[which(snps$eff=="synonymous_variant")] <- "Synonymous"
snps$eff[which(snps$eff=="stop_lost&splice_region_variant")] <- "Stop_lost"
snps$eff[which(snps$eff=="stop_lost")] <- "Stop_lost"
snps$eff[which(snps$eff=="stop_lost")] <- "Stop_lost"
snps$eff[which(snps$eff=="stop_gained")] <- "Stop_gained"
snps$eff[which(snps$eff=="No annotation")] <- "No Annotation"
snps$eff[which(snps$eff=="start_lost")] <- "Start_lost"
snps$eff[which(snps$eff=="frameshift_variant&stop_lost&splice_region_variant")] <- "Frameshift"

tabeffect <- as.data.frame(table(snps$eff)/length(data$country))
tabeffect$country="Global"

piechartmut <- data.frame(Country=character(0),Consequence=character(0),Number=numeric(0))
piechartmut <- rbind(piechartmut,tabeffect[,c(3,1,2)])

for (country in unique(snps$country)) {
  temp = as.data.frame(table(snps$eff[which(snps$country==country)])/length(which(data$country==country)))
  temp$country=country
  piechartmut <- rbind(piechartmut,temp[,c(3,1,2)])
}  
write.csv(piechartmut,paste(DIR,"files4dashboard/average_consequence_percountry_",date,"_",REF,".csv",sep=""),row.names = F,quote = F)


piechartmut <- data.frame(Country=character(0),Consequence=character(0),Number=numeric(0),Gene=character(0))
listofgenes <- orfs$V3
for (gene in listofgenes) {
  snps2use = snps[which(snps$orf==gene),]
  tabeffect <- as.data.frame(table(snps2use$eff)/length(data$country))
  tabeffect$country="Global"
  
  piechartmut1 <- data.frame(Country=character(0),Consequence=character(0),Number=numeric(0))
  piechartmut1 <- rbind(piechartmut1,tabeffect[,c(3,1,2)])
  
  for (country in unique(snps2use$country)) {
    temp = as.data.frame(table(snps2use$eff[which(snps2use$country==country)])/length(which(data$country==country)))
    temp$country=country
    piechartmut1 <- rbind(piechartmut1,temp[,c(3,1,2)])
  }  
  piechartmut1$Gene = gene
  piechartmut <- rbind(piechartmut,piechartmut1)
}
write.csv(piechartmut,paste(DIR,"files4dashboard/average_consequence_percountry_",date,"_",REF,"_pergene.csv",sep=""),row.names = F,quote = F)
###############################################
####### Part3: Create table ###################
###############################################
snps <- snps[which(snps$mut%in%mut$Var1),]
df <- snps[,c(2,3,9,8,11)]
colnames(df) <- c("Country","Date","Mutation","Gene","Consequence")
write.csv(df,paste(DIR,"files4dashboard/covid19-genomics-dataset-",date,".csv",sep=""),col.names = T,row.names = F,quote =F)

