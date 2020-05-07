require(RCurl)
require(dplyr)

########LIST OF FUNCTIONS###########


pNCBISummary=function(input.file){ #this partitions the NCBI summary file into organisms and protein accesion. It compares the organism lists to find LapGD organisms.
  data.input<-readLines(input.file)
  input.lines<-which(grepl("[", data.input, fixed=TRUE)==TRUE)
  access.lines<-input.lines+2
  access.lines<-strsplit(data.input[access.lines]," ", fixed=TRUE)
  seperate<-strsplit(data.input[input.lines], "[", fixed=TRUE)
  
  ORG=NULL
  ACCESS=NULL
  
  for (i in 1:length(access.lines)){
    org<-substr(seperate[[i]][2],0,nchar(seperate[[i]][2])-1)
    ORG<-c(ORG,org)
    access<-access.lines[[i]][1]
    ACCESS<-c(ACCESS,access)
  }
  out<-data.frame(ORG,ACCESS,stringsAsFactors = FALSE)
  
  return(out)
}


searchable.faa=function(faa.file){ #converts genome faa file into a table with the accession, AA squence, putative function, and AA size. can pull info from this table about RTX status, cleavage site, etc
  test1<-readLines(faa.file)
  inputs<-which(substr(test1,0,1)==">") 
  
  ACC<-NULL
  SEQQ<-NULL
  pFUN<-NULL
  SIZE<-NULL
  for (i in 1:length(inputs)){
    slp.inputs<-gregexpr(" ", test1[inputs[i]], fixed=TRUE)
    species.split<-gregexpr("[", test1[inputs[i]], fixed=TRUE)
    acc<-substr(test1[inputs[i]], 2, slp.inputs[[1]][1]-1)
    p.fun<-substr(test1[inputs[i]],slp.inputs[[1]][1]+1, species.split[[1]]-1)
    first<-inputs[i]+1
    
    if (i+1>length(inputs)){
      last<-length(test1)
    } else { last<-inputs[i+1]-1 }
    
    seq=NULL
    
    for(a in first:last){
      seq<-paste(seq, test1[a], sep="")
      size<-nchar(seq)
    }
    
    ACC<-c(ACC,acc)
    pFUN<-c(pFUN,p.fun)
    SEQQ<-c(SEQQ,seq)
    SIZE<-c(SIZE,size)
  }
  
  tbl<-data.frame(ACC,pFUN, SEQQ, SIZE, stringsAsFactors = FALSE)
  names(tbl)<-c("accession", "putative.function","protein.sequence", "protein.length")
  return(tbl)
}

verify.files=function(input.url){
  require("RCurl")
  temp.url<-paste(input.url, "/", sep="")
  getfiles<-getURL(temp.url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  getfiles<-strsplit(getfiles,"\n")
  sp<-strsplit(input.url,"/", fixed=TRUE)
  get.l<-length(sp[[1]])
  temp.protein<-paste(sp[[1]][get.l], "_protein.faa.gz", sep="")
  temp.table<-paste(sp[[1]][get.l], "_feature_table.txt.gz", sep="")
  find.table<-any(getfiles[[1]]==temp.table)
  find.protein<-any(getfiles[[1]]==temp.protein)
  if (find.table==TRUE & find.protein==TRUE){
    return(TRUE)
  } else { return(FALSE)}
}

s.url=function(url){
  #  verify.files(url)
  sp<-strsplit(url,"/", fixed=TRUE)
  get.l<-length(sp[[1]])
  new.url<-paste(url, "/", sp[[1]][get.l], sep="")
  return(new.url)
}

any.RTX=function(y){ #Returns proteins with RTX motifs
  #D-x-[LI]-x(4)-G-x-D-x-[LI]-x-G-G-x(3)-D
  RTX<-c("D[[:alpha:]]L[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]L[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]I[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]L[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]L[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]I[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D",
         "D[[:alpha:]]I[[:alpha:]][[:alpha:]][[:alpha:]][[:alpha:]]G[[:alpha:]]D[[:alpha:]]I[[:alpha:]]GG[[:alpha:]][[:alpha:]][[:alpha:]]D")
  
  checker<-lapply(RTX, function(x) grepl(x,y))
  results<-any(checker==TRUE)
  return(results)
}

find.cleavage=function(anchor){ 
  ## Searches for the following previously characterized LapG cleavage motifs and returns T or F for whatever AA string you feed it. 
  ## In the script below it's given residues 70-150 of proteins that have RTX motifs and are at least 1000 AA long
  
  TAAG<-grepl("TAAG", anchor)
  PAAG<-grepl("PAAG", anchor)
  AAAG<-grepl("AAAG", anchor)
  TAAV<-grepl("TAAV", anchor)
  
  
  if (TAAG==TRUE|PAAG==TRUE|AAAG==TRUE|TAAV==TRUE){
    return(TRUE)} else { return(FALSE) }
}


##################################
##################################
##################################
##################################

#######STEP 1 PROCESS NCBI SUMMARY FILES####### 
#Download LapD_MoxY_N domain containing strain summary list as "LapD.txt": https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=318615

#Download LapG comtaining strain summary list as "LapG.txt": https://www.ncbi.nlm.nih.gov/sites/entrez?LinkName=cdd_protein&db=cdd&cmd=Link&from_uid=310549

###if you use another format the pNCBI summary function won't work###

####CONVERT SUMMARY FILES#####
pSUM.lapD<-pNCBISummary("Data/LapD.txt")
pSUM.lapG<-pNCBISummary("Data/LapG.txt")


# 
##############################

#######STEP 2 - Identify LapGD-encoding bugs, do a little analysis########
lapD.organism.tally<- length(unique(pSUM.lapD$ORG)) #how many unique genomes have LapD-like protein? 4484
lapG.organism.tally<- length(unique(pSUM.lapG$ORG)) #how many unique genomes have LapG-like protein? 6455
lapGD.organisms<-unique(intersect(pSUM.lapG$ORG,pSUM.lapD$ORG)) #select organisms found on both lists. 3993

#############################################################

#######STEP 3 - Read Downloaded Genome list file#######
all.genomes<-read.csv("Data/ALL_genomes_proks.csv", stringsAsFactors = FALSE, header=TRUE) #file contains all bacterial genomes as of 4/13/2020 downloaded from https://www.ncbi.nlm.nih.gov/assembly?term=%28%22Bacteria%22%5BOrganism%5D%29%20AND%20%28bacteria%5Bfilter%5D%20AND%20%28latest%5Bfilter%5D%20OR%20%22latest%20refseq%22%5Bfilter%5D%29%20AND%20%28all%5Bfilter%5D%20NOT%20%22derived%20from%20surveillance%20project%22%5Bfilter%5D%20AND%20all%5Bfilter%5D%20NOT%20anomalous%5Bfilter%5D%29%29&cmd=DetailsSearch

#############################################################

#######STEP 4 #######

###################GET FTP INFO FOR LAPGD BUGS###################

#############create dataframe containing all the lapD/G encoding species with subgroup and sequence info################


target.rows<-lapply(lapGD.organisms, function(x) which(all.genomes$Organism.Name==x))
get.rows<-sapply(target.rows, function(x) x[1]) #returns only one value per 
get.rows<-na.omit(get.rows)


lapd.g.species.info<-data.frame(all.genomes[get.rows,c(1,2,6,18,21)]) #create dataset of only organisms and info of interest
names(lapd.g.species.info)<-c("organism", "TaxID", "subgroup", "level", "RefSeq.FTP") #name cols

#############################################################

#######STEP 5#######

#############create data frame containing species names and protein and table source files######################

tbl<-sapply(lapd.g.species.info$RefSeq.FTP, function(x) paste(s.url(x), "_feature_table.txt.gz", sep=""))
prt<-sapply(lapd.g.species.info$RefSeq.FTP, function(x) paste(s.url(x), "_protein.faa.gz", sep=""))
new.d<-data.frame(lapd.g.species.info,tbl,prt)

#####download FAA files######

#sometimes the download is interupted or individual records aren't downloaded. errors object stores any that fail to download so they can be manually checked.

errors <- list()

for (i in 1:length(new.d$organism)){ #
  filename = paste("./Output/Genomes/", new.d$organism[i], "_protein.faa.gz", sep="")
  if(file.exists(filename)){
    next
  }
  tryCatch({
    message("Downloading: ", new.d$prt[i],"\n")
    download.file(as.character(new.d$prt[i]), destfile = filename)   #download genome faa
    message("Download Complete.\n\n")
  },error = function(e){ 
    errors <<- c(errors, new.d$organism[i])
  })
}


#######STEP 6 #######

######ONCE YOU HAVE DOWNLOADED THE GENOME FAA FILES YOU CAN LOOK THROUGH THEM. THIS WILL SELECT LARGE PROTEINS WITH RTX MOTIFS AND A LAPG CLEAVAGE SITE###########
#######GENERATE LIST OF LARGE RTX PROTEINS. THEN LOOK FOR CLEAVAGE. SEPERATE BY CLEAVAGE STATUS###############

setwd("./Output/Genomes/") #Set working directory to location of protein records

faa.files<-list.files() 
my.table<-data.frame(ACC=character(), pFUN=character(), SEQQ=character(), SIZE=numeric(), stringsAsFactors=FALSE) 

for (i in 1:length(faa.files)){
  limit.size <- NULL
  closeAllConnections()
  species.name<-strsplit(faa.files[i], "_protein.faa.gz") #get species name
  species.name<-species.name[[1]][1]
  
  message("processing....", species.name, " faa file\t", i, " of ", length(faa.files))
  
  f.proteins<-searchable.faa(gzfile(faa.files[i])) #convert faa file to searchable format.
  
  RTX.status<-lapply(f.proteins$protein.sequence, function(x) any.RTX(x)) #look through protein sequence column for proteins with RTX
  
  
  if (any(RTX.status==TRUE)==TRUE) { ##only interogate genomes with large RTX adhesins.
    putative.RTX<-which(RTX.status==TRUE)
    limit.size<-subset(f.proteins[putative.RTX,], protein.length>=1000)
    org.rep<-rep(species.name, length(limit.size$protein.length))
    limit.size<-data.frame(organism=org.rep,limit.size)
  }
  
  my.table<-rbind(my.table, limit.size) 
  message("complete\n\n")
  
}
motif_status<-lapply(my.table$protein.sequence, function(x) find.cleavage(substr(x, 70,150))) ##search list of large RTX proteins for LapG cleavage site. Return T or F for putative cleavage site.
new_my.table<-mutate(my.table, lapG=motif_status)
new_my.table$lapG <- as.character(new_my.table$lapG)

setwd("../..") #Set working directory back to parent directory

new_my.table.unique <- new_my.table[!duplicated(new_my.table$accession),]

closeAllConnections()


substrates2<-subset(new_my.table[which(new_my.table$lapG==TRUE),1:5], stringsAsFactors=FALSE)
no.substrates2<-subset(new_my.table[which(new_my.table$lapG==FALSE),])

substrates.unique <- substrates2[!duplicated(substrates2$accession),]
non.substrates.unique <- no.substrates2[!duplicated(no.substrates2$accession),]
non.substrates.unique$organism <- as.character(non.substrates.unique$organism)
row.names(non.substrates.unique) <- seq(length(row.names(non.substrates.unique)))
substrates.unique$organism <- as.character(substrates.unique$organism)
row.names(substrates.unique) <- seq(length(row.names(substrates.unique)))
protein.count <- data.frame(table(unlist(substrates.unique$organism)))
non.lapA.but.lapDG.orgs <- data.frame(table(unlist(non.substrates.unique$organism)))

write.csv(substrates.unique, "Output/All_lapg_targets.csv", row.names = FALSE)
write.csv(no.substrates2, "Output/all_non_lapg_targets.csv", row.names = FALSE)
