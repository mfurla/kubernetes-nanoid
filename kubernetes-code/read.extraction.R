### environment set-up ###
source("environmentSetUp.R")

### load objects ###
load(file.path("data","fast5.files.RData"))

load(file=file.path("data","bam.RData"))
bamNames <- names(bam)
rm(bam)
gc()

### for any read we associate the fast5 file and the nucleotide sequence predicted by guppy ###

read.to.fast5.name = c()
read.sequence.list = list()
saveList(object = read.sequence.list, file = "data/read.sequence.list.llo", append = FALSE, compress = FALSE)

for(fast5.file in fast5.files){
  try({
	## extraction of fast5 reads names
	fast5ls <- h5ls(fast5.file)
	fast5reads <- unique(fast5ls$name[grep("read_",fast5ls$name)])

	print(fast5.file)
	## instead of a for we use a mclapply to parallelize
	tmp <- unlist(mclapply(fast5reads, function(fast5read)
	{
		h5 <- h5read(fast5.file,paste0(fast5read,"/Analyses/Basecall_1D_001/BaseCalled_template/Fastq"))

		fast5read <- strsplit(fast5read,"_")[[1]][[2]]

		out <- gsub("U","T",strsplit(h5,split = "\n")[[1]][2])
		names(out) <- fast5.file

		out
	},mc.cores=mc.cores))

 	namesTmp <- names(tmp)
 	names(namesTmp) <- sapply(strsplit(fast5reads,"_"),function(i)i[[2]])
 	names(tmp) <- sapply(strsplit(fast5reads,"_"),function(i)i[[2]])

  	read.to.fast5.name = c(read.to.fast5.name, namesTmp)
  	# read.sequence.list = c(read.sequence.list, tmp)

  	saveList(object = as.list(tmp), file = "data/read.sequence.list.llo", append = TRUE, compress = FALSE)

  },silent = TRUE)
}

### results save
save(read.to.fast5.name,file=file.path("data","read.to.fast5.name.RData"))

read.sequence.list <- readList(file = "data/read.sequence.list.llo")
names(read.sequence.list) <- names(read.to.fast5.name)

# MF I decide to keep only the aligned reads
#
read.to.fast5.name <- read.to.fast5.name[bamNames]
read.sequence.list <- read.sequence.list[bamNames]

save(read.to.fast5.name,file=file.path("data","read.to.fast5.name.RData"))
save(read.sequence.list,file=file.path("data","read.sequence.list.RData"))
saveList(read.sequence.list,"data/read.sequence.list.llo")

# system("rm data/read.sequence.list.llo")
