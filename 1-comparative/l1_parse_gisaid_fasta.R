library(stringr)
library(readxl)

dir.create("SubFiles/")
dir.create("FileLists/")

fdate = paste(unlist(str_split(Sys.Date(), pattern = "-"))[c(3,2,1)], collapse = "")

# update the merged acknowlegement table
if (file.exists("AcknowledgementTables/COMPLETE_acknowledgement_table.xlsx")){
  oldfile <- read_xlsx("AcknowledgementTables/COMPLETE_acknowledgement_table.xlsx")
  currentfile.name <- list.files(path = "AcknowledgementTables/", pattern = paste("gisaid_hcov-19_acknowledgement_table_", gsub(Sys.Date(), pattern = "-", replacement = "_") ,sep = "") )
  currentfile <- read_xls(paste("AcknowledgementTables/", currentfile.name, sep = ""), skip =2)[-1,]
  writexl::write_xlsx(unique(rbind(oldfile, currentfile)), "AcknowledgementTables/COMPLETE_acknowledgement_table.xlsx")
  rm(oldfile, currentfile, currentfile.name)
}
# update the merged patient status metadata

if(file.exists("MetaData/COMPLETE_metadata.txt")){
  oldfile <- read.delim("MetaData/COMPLETE_metadata.txt", sep= "\t")
  currentfile.name <- list.files(path = "MetaData/", pattern = paste("gisaid_hcov-19_", gsub(Sys.Date(), pattern = "-", replacement = "_") ,sep = "") )
  currentfile <- read.delim(paste("MetaData/", currentfile.name, sep = ""), skip =2)[-1,]
  write.table(unique(rbind(oldfile, currentfile)), "MetaData/COMPLETE_metadata.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  rm(oldfile, currentfile, currentfile.name)
}

filename <- list.files(path = "MainFile/", pattern = paste("gisaid_hcov-19_", gsub(Sys.Date(), pattern = "-", replacement = "_") , ".*fasta$",sep = ""))
lines <- readLines(paste("MainFile/", filename, sep= ""), skip =0)

ind <- grep("^>",lines)
length(ind)

virus <- unlist(lapply(lines[ind],function(x){unlist(strsplit(x,"/"))[1]}))
unique(virus)
country <- unlist(lapply(lines[ind],function(x){unlist(strsplit(x,"/"))[2]}))
table(country)
identifier <- unlist(lapply(lines[ind],function(x){unlist(strsplit(x,"/"))[3]}))
names4files <- unlist(lapply(lines[ind],function(x){unlist(strsplit(unlist(strsplit(x,"/"))[4],"\\|"))[2]}))                    
date <- unlist(lapply(lines[ind],function(x){unlist(strsplit(unlist(strsplit(x,"/"))[4],"\\|"))[3]}))                    
#### Manual curation 

if(length(which(is.na(date)))>0 || length(which(is.na(names4files)))>0){
  print(" missing info ")
  quit(1)
}

# lines[ind[which(is.na(date))]]
# i = 1
# lines[ind[which(is.na(date))[i]]]
# date[which(is.na(date))[i]] <- "2020-03-10"
# lines[ind[which(is.na(date))]]
# i=10
# lines[ind[which(is.na(date))[i]]]
# date[which(is.na(date))[i]] <- "2019-12-31"
# lines[ind[which(is.na(date))]]
# 
# lines[ind[which(is.na(names4files))]]
# i=1
# names4files[which(is.na(names4files))[i]] <- "EPI_ISL_418800"  #where do i get this from??
# lines[ind[which(is.na(names4files))]]
# i=10
# names4files[which(is.na(names4files))[i]] <- "EPI_ISL_402125"
# lines[ind[which(is.na(names4files))]]


mini <- ind+1
maxi <- ind-1
maxi <- c(maxi[2:length(maxi)],length(lines))
seqob <- c(1:length(ind))

dir.create(paste("SubFiles/",fdate, sep = ""))
for (i in 1:length(ind)) {
  if(is.na(names4files[i])) {
    cat(names4files[i],"\n")
    next
  }
  i1 = mini[i]-1
  i2 = maxi[i]
  file2create = paste(names4files[i],"_",fdate,".fasta",sep="")
  print(file2create)
  # i1-1 is to include the fasta header
  writeLines(text = lines[i1:i2],paste("SubFiles/",fdate,"/",file2create,sep=""),sep = "\n")

}
# seperated every fatsa file into a new file???????



#ind <- ind[which(!is.na(names4files))]

names4files <- names4files[which(!is.na(names4files))]
write.csv(names4files,paste("FileLists/gisaid_files_list_",fdate,".csv",sep=""),quote = F,row.names = F)

