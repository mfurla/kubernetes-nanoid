### environment set-up ###
source("environmentSetUp.R")

### load objects ###
load(file.path("data","fast5.files.RData"))

### extract aligned sequences ###
reference.genome = readDNAStringSet(filepath = "Mus_musculus.GRCm38.dna.primary_assembly.fa")
## MF: I don't know why but the chromosomes names in the reference genome are strange...
names(reference.genome) <- sapply(strsplit(names(reference.genome)," "),"[[",1)

bam = readGAlignments("data/minimap.bam",param = NULL,use.names = TRUE)
bam = bam[isUnique(names(bam))]

aligned.sequence.list = as.list(reference.genome[as(bam,"GRanges")])
aligned.sequence.list = lapply(aligned.sequence.list,as.character)
names(aligned.sequence.list) = names(bam)

cigar.strings = cigar(bam)
cigar.ops = explodeCigarOps(cigar.strings)
cigar.list = explodeCigarOpLengths(cigar.strings)
cigar.list = lapply(1:length(cigar.strings),function(x){vec = cigar.list[[x]];names(vec) = cigar.ops[[x]];vec})
names(cigar.list) = names(bam)

save(bam,file=file.path("data","bam.RData"))
save(aligned.sequence.list,file=file.path("data","aligned.sequence.list.RData"))
save(cigar.list,file=file.path("data","cigar.list.RData"))

saveList(aligned.sequence.list,"data/aligned.sequence.list.llo")
