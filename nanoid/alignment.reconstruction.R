### environment set-up ###
source("environmentSetUp.R")

### alignment.reconstruction.function ###
alignment.reconstruction.function = function(read.name){
  fast5.name = read.to.fast5.name[read.name]
  read = read.sequence.list[[read.name]]
  which.strand = as.character(strand(bam[read.name]))
  if (which.strand == "-"){read = as.character(reverseComplement(as(read,'DNAString')))}
  reference = aligned.sequence.list[[read.name]]
  cigar = cigar.list[[read.name]]
  exploded.cigar = as.vector(Rle(names(cigar),cigar))
  alignment.length = length(exploded.cigar)
  aligned.read = rep("-",alignment.length)
  aligned.read[which(exploded.cigar %in% c("M","I","S","H"))] = strsplit(read,split = "")[[1]]
  aligned.reference = rep("-",alignment.length)
  aligned.reference[which(exploded.cigar %in% c("M","D","N"))] = strsplit(reference,split = "")[[1]]
  alignment = rbind(aligned.read,aligned.reference)
  colnames(alignment) = exploded.cigar
  alignment.uncut = alignment
  alignment = alignment[,which(colnames(alignment) != "H"),drop = FALSE] # cigar H is for hard-clipping
  alignment = alignment[,which(colnames(alignment) != "S"),drop = FALSE] # cigar S is for soft-clipping
  alignment = alignment[,which(colnames(alignment) != "N"),drop = FALSE] # cigar N is for splicing
  alignment = alignment[,which(alignment["aligned.reference",] != "-"),drop = FALSE] # "-" is for insertions in the reference
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  return(alignment)
}

### alignment.reconstruction ###
bam = get(load(file.path("data","bam.RData")))
cigar.list = get(load(file.path("data","cigar.list.RData")))
read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))

# aligned.sequence.list = get(load(file.path("data","aligned.sequence.list.RData")))
# read.sequence.list = get(load(file.path("data","read.sequence.list.RData")))

bamNames <- names(bam)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

if(file.exists(where="data/alignment.reconstruction.list.llo"))
{
   start = getListLength("data/alignment.reconstruction.list.llo")/100 + 1
   print(paste0("First index: ",start))
}else{
  alignment.reconstruction = list()
  saveList(object = alignment.reconstruction, file = "data/alignment.reconstruction.list.llo", append = FALSE, compress = FALSE)
  start = 1
}

registerDoParallel(cores = mc.cores)

if(!file.exists("aligned.sequence.list.llo"))
{
  load("data/aligned.sequence.list.RData")
  saveList(aligned.sequence.list,"data/aligned.sequence.list.llo")
  rm(aligned.sequence.list)
  gc()
}

for (j in start:length(index.list)){
  if(j%%50==0){print(j);system("grep MemFree /proc/meminfo")}

  aligned.sequence.list <- readList("data/aligned.sequence.list.llo",index.list[[j]])
  names(aligned.sequence.list) <- bamNames[index.list[[j]]]

  read.sequence.list <- readList("data/read.sequence.list.llo",index.list[[j]])
  names(read.sequence.list) <- bamNames[index.list[[j]]]

  alignment.reconstruction.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% alignment.reconstruction.function(n)
  saveList(object = alignment.reconstruction.tmp, file = "data/alignment.reconstruction.list.llo", append = TRUE, compress = FALSE)

  rm(alignment.reconstruction.tmp)
  rm(aligned.sequence.list)
  rm(read.sequence.list)
  gc()
}
alignment.reconstruction <- readList(file = "data/alignment.reconstruction.list.llo")
names(alignment.reconstruction) = bamNames

save(alignment.reconstruction,file=file.path("data","alignment.reconstruction.RData"))
saveList(alignment.reconstruction,"data/alignment.reconstruction.list.llo")
