### environment set-up ###
source("environmentSetUp.R")

### load objects ###
load(file.path("data","fast5.files.RData"))

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
aligned.sequence.list = get(load(file.path("data","aligned.sequence.list.RData")))
cigar.list = get(load(file.path("data","cigar.list.RData")))
read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("data","read.sequence.list.RData")))

bamNames <- names(bam)

registerDoParallel(cores = mc.cores)
alignment.reconstruction = foreach(n = bamNames,.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.sequence.list","read.to.fast5.name"))) %dopar% alignment.reconstruction.function(n)
names(alignment.reconstruction) = bamNames

save(alignment.reconstruction,file=file.path("data","alignment.reconstruction.RData"))
saveList(alignment.reconstruction,"data/alignment.reconstruction.list.llo")