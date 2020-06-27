
datefile = paste(unlist(str_split(Sys.Date(), pattern = "-"))[c(3,2,1)], collapse = "")
REF ="EPI_ISL_402125"

dir.create("Stats/")
dir.create("VCFs")


names4files <- as.vector(read.csv(paste("FileLists/gisaid_files_list_",datefile,".csv",sep=""),stringsAsFactors = F))       
names4files <- names4files$x[which(names4files$x!=REF)] 
df = data.frame(id=REF,country="Wuhan",date="2019-12-31",length=29903,Ref.alignedbases=100,Qry.alignedbases=100,
                breakpoints=0,relocation=0,translocation=0,inversion=0,inertion=0,snps=0,indels=0)
write.table(df,paste("Stats/summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)

for (QRY in names4files) {
  cat(QRY,"\n")
  filename <- (paste("SubFiles/",datefile,"/",QRY,"_",datefile,".fasta",sep=""))
  lines <- readLines(filename,n=1,skip =0)
  country <- unlist(lapply(lines,function(x){unlist(strsplit(x,"/"))[2]}))
  date <- unlist(lapply(lines,function(x){unlist(strsplit(unlist(strsplit(x,"/"))[4],"\\|"))[3]}))  
  #put each set of nucmers in a folder per dat and add path after Analysis Files
  reportfile <- paste("AnalysisFiles/",datefile,"/dnadiff_nucmer_",REF,"_",QRY,"_",datefile,".report",sep="")
  #reportfile <- paste("dnadiff_nucmer_",REF,"_",QRY,"_",datefile,".report",sep="")
  if (file.exists(reportfile)) {
    test <- readLines(reportfile,skip = 0)
    length = as.numeric(unlist(strsplit(test[grep("TotalBases",test)]," "))[length(unlist(strsplit(test[grep("TotalBases",test)]," ")))])
    Ref.alignedbases = as.numeric(unlist(strsplit(unlist(strsplit(test[grep("AlignedBases",test)],"\\("))[2],"%)"))[1])
    Qry.alignedbases= as.numeric(unlist(strsplit(unlist(strsplit(test[grep("AlignedBases",test)],"\\("))[3],"%)"))[1])
    brkpts = as.numeric(unlist(strsplit(test[grep("Breakpoints",test)],"\\s "))[length(unlist(strsplit(test[grep("Breakpoints",test)],"\\s ")))])
    reloc = as.numeric(unlist(strsplit(test[grep("Relocations",test)],"\\s "))[length(unlist(strsplit(test[grep("Relocations",test)],"\\s ")))])
    transloc = as.numeric(unlist(strsplit(test[grep("Translocations",test)],"\\s "))[length(unlist(strsplit(test[grep("Translocations",test)],"\\s ")))])
    inv = as.numeric(unlist(strsplit(test[grep("Inversions",test)],"\\s "))[length(unlist(strsplit(test[grep("Inversions",test)],"\\s ")))])
    ins = as.numeric(unlist(strsplit(test[grep("Insertions",test)],"\\s "))[length(unlist(strsplit(test[grep("Insertions",test)],"\\s ")))])
    snps = as.numeric(unlist(strsplit(test[grep("TotalSNPs",test)],"\\s "))[length(unlist(strsplit(test[grep("TotalSNPs",test)],"\\s ")))])
    indels = as.numeric(unlist(strsplit(test[grep("TotalIndels",test)],"\\s "))[length(unlist(strsplit(test[grep("TotalIndels",test)],"\\s ")))])
    if(length<29000){
      next
    } else {
      df = data.frame(id=QRY,country=country,date=date,length=length,Ref.alignedbases=Ref.alignedbases,Qry.alignedbases=Qry.alignedbases,
           breakpoints=brkpts,relocation=reloc,translocation=transloc,inversion=inv,inertion=ins,snps=snps,indels=indels)
      write.table(df,paste("Stats/summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = F,append = T)
    }
  } else {
    cat("Skipping missing report file ", QRY, "\n")
  }
}


data <- read.table(paste("Stats/summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)
dim(data)
summary(data)
#### Manual curation
# data$country[which(data$country=="NetherlandsL")] <- "Netherlands"
# data$country[which(data$country=="New Zealand")] <- "NewZealand"
# data$country[which(data$country=="Northern Ireland")] <- "Ireland"
# data$country[which(data$country=="Wuhan-Hu-1")] <- "Wuhan"
# data$country[which(data$country=="NanChang")] <- "Nanchang"
# data$country[which(data$country=="DK")] <- "Denmark"
# data$country[which(data$country=="ITALY")] <- "Italy"
data$country<-replace(as.character(data$country), data$country == "Northern Ireland", "Ireland")
data$country<-replace(as.character(data$country), data$country == "New Zealand", "NewZealand")
data$country<-replace(as.character(data$country), data$country == "ITALY", "Italy")
data$country<-replace(as.character(data$country), data$country == "Wuhan-Hu-1", "Wuhan")
data$country<-replace(as.character(data$country), data$country == "DK", "Denmark")
data$country<-replace(as.character(data$country), data$country == "NetherlandsL", "Netherlands")
data$country<-replace(as.character(data$country), data$country == "NanChang", "Nanchang")




write.table(data,paste("Stats/summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)


for (QRY in data$id) {
  cat(QRY,"\n")
  country <- data$country[which(data$id == QRY)]
  date <- data$date[which(data$id == QRY)]
  if (data$snps[which(data$id==QRY)]) {
    reportfile <- paste("AnalysisFiles/", datefile,"/dnadiff_nucmer_",REF,"_",QRY,"_",datefile,".snps",sep="")
    snps <- read.table(reportfile,header = F,stringsAsFactors = F)
    n =dim(snps)[1]
    df = data.frame(id=rep(QRY,n),country=rep(country,n),date=rep(date,n),pos.REF=snps[,1],sub.REF=snps[,2],sub.QRY=snps[,3],pos.QRY=snps[,4])
    write.table(df,paste("Stats/snps_summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = F,append = T)
  } else {
    cat("Skipping missing report file ", QRY, "\n")
  }
}

snps <- read.table(paste("Stats/snps_summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",header = F,stringsAsFactors = F)
colnames(snps) <- c("id","country","date","Ref.pos","Ref.sub","Qry.sub","Qry.pos")



orfs <- read.table("orfs_covid19.txt",stringsAsFactors = F,header = F)
br <- c(orfs$V1,30000)
snps$orf <- cut(snps$Ref.pos, breaks=br,labels = orfs$V3,right = F) 
snps$mut <- paste(snps$Ref.pos,snps$Ref.sub,">",snps$Qry.sub,sep = "")
snps$mutnuc <- paste(snps$Ref.sub,">",snps$Qry.sub,sep = "")

#write.table(snps,paste("Stats/snps_summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)
orfs <- read.table("orfs_covid19.txt",stringsAsFactors = F,header = F)
br <- (sort(unique(c(0,as.vector(orfs$V1),as.vector(orfs$V2),30000))))
lab <- c("other","ORF1A","ORF1B","other","S","other","ORF3A","other","E","other","M","other","ORF6","other","ORF7A","ORF7B","other","ORF8","other","N","other","ORF10","other")
s1 = 30000-sum(orfs$V4)
orfs<-rbind(orfs,c(1,1,"other",s1)) 
row.names(orfs) <- orfs$V3
snps$orf <- cut(snps$Ref.pos, breaks=br,labels =lab,right = F) 
write.table(snps,paste("Stats/snps_summaryfile_",datefile,"_",REF,".txt",sep=""),sep="\t",quote = F,row.names = F,col.names = T)



#### Prepare vcf file
# get_mutation <- function(x) {
#   if (length(which(mutations$mut==x))) {
#     as.character(mutations$type[which(mutations$mut==x)])[1]
#   } else {
#     "NAN"
#   }
# }


#fix the TRUE in the snp summary file


date = paste(unlist(str_split(Sys.Date(), pattern = "-"))[c(3,2,1)], collapse = "")
REF ="EPI_ISL_402125"
data <- read.table(paste("summaryfile_",date,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)
snps <- read.table(paste("Stats/snps_summaryfile_",date,"_",REF,".txt",sep=""),sep="\t",header = T,stringsAsFactors = F)

tabmut <- as.data.frame(table(snps$mut))
#tabmut$effect <- unlist(lapply(as.character(tabmut$Var1),FUN = get_mutation))

ind = match(unique(snps$mut), snps$mut)
n=length(ind)
vcffile <- data.frame(CHROM=rep("NC_045512.2",n),	POS = snps$Ref.pos[ind],	ID	= rep(".",n), REF=snps$Ref.sub[ind],	ALT=snps$Qry.sub[ind]	,QUAL= rep(".",n)	,FILTER	= rep(".",n),INFO= rep(".",n))

vcffile$REF[which(vcffile$REF=="TRUE")]="T"
vcffile$ALT[which(vcffile$ALT=="TRUE")]="T"
x_indeces <- which(vcffile$ALT=="X")
if (length(x_indeces)>0){
  vcffile <- vcffile[-which(vcffile$ALT=="X"),]
}

write.table(vcffile,paste("VCFs/mutations_",date,".vcf", sep =  ""),quote = F,row.names = F,col.names = T,sep="\t")
