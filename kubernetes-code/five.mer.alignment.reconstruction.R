### environment set-up ###
source("environmentSetUp.R")

# MF: all these objects are required by the following function but never loaded
load(file.path("Supplemental_Code_S7/BasicObjects/","five.mers.minus.RData"))
load(file.path("Supplemental_Code_S7/BasicObjects/","cigar.five.mers.RData"))
load(file.path("Supplemental_Code_S7/BasicObjects/","cigar.five.mers.S.containing.RData"))
load(file.path("Supplemental_Code_S7/BasicObjects/","cigar.five.mers.H.S.N.I.centered.RData"))

### five.mer.alignment.reconstruction.function ###
five.mer.alignment.reconstruction.function = function(read.name){
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
  
  if (which.strand == "-"){
    complement.map = c("A" = "T","C" = "G","G" = "C","T" = "A","N" = "N","-" = "-")
    
    alignment["aligned.read",] = complement.map[alignment["aligned.read",]]
    alignment["aligned.reference",] = complement.map[alignment["aligned.reference",]]
  }
  
  aligned.read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.read",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  aligned.reference.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[alignment["aligned.reference",which(!colnames(alignment) %in% c("H","S"))]]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
  
  exploded.cigar.sequence.five.mers = as.character(cigar.five.mers[as.numeric(stats::filter(as.numeric(c("M" = 1,"I" = 2,"D" = 3,"H" = 4,"S" = 5,"N" = 6)[exploded.cigar]),6^(0:8)[1:5],sides = 1)-(sum(6^(0:8)[1:5])-1))[-c(1:4)]])
  
  five.mers.alignment = rbind(exploded.cigar.sequence.five.mers,
                              "aligned.read.sequence.five.mers" = NA,
                              "aligned.reference.sequence.five.mers" = NA)
  five.mers.alignment["aligned.read.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.read.sequence.five.mers
  five.mers.alignment["aligned.reference.sequence.five.mers",!(exploded.cigar.sequence.five.mers %in% c(cigar.five.mers.S.containing))] = aligned.reference.sequence.five.mers
  five.mers.alignment = five.mers.alignment[,!(exploded.cigar.sequence.five.mers %in% cigar.five.mers.H.S.N.I.centered)]
  
  return(five.mers.alignment)
}

### five.mer.alignment.reconstruction ###
bam = get(load(file.path("data","bam.RData")))
aligned.sequence.list = get(load(file.path("data","aligned.sequence.list.RData")))
cigar.list = get(load(file.path("data","cigar.list.RData")))
read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("data","read.sequence.list.RData")))

five.mer.alignment.reconstruction.list = list()
saveList(object = five.mer.alignment.reconstruction.list, file = "data/five.mer.alignment.reconstruction.list.llo", append = FALSE, compress = FALSE)

bamNames <- names(bam)

index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

registerDoParallel(cores = mc.cores)
for (j in 1:length(index.list)){
  	five.mer.alignment.reconstruction.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","aligned.sequence.list","cigar.list","read.sequence.list","read.to.fast5.name"))) %dopar% five.mer.alignment.reconstruction.function(n)
	saveList(object = five.mer.alignment.reconstruction.list.tmp, file = "data/five.mer.alignment.reconstruction.list.llo", append = TRUE, compress = FALSE)
}

five.mer.alignment.reconstruction.list = readList("data/five.mer.alignment.reconstruction.list.llo")
names(five.mer.alignment.reconstruction.list) <- bamNames

save(five.mer.alignment.reconstruction.list,file=file.path("data","five.mer.alignment.reconstruction.RData"))
saveList(object = five.mer.alignment.reconstruction.list, file = "data/five.mer.alignment.reconstruction.list.llo", append = FALSE, compress = FALSE)
# system("rm data/five.mer.alignment.reconstruction.list.llo")
