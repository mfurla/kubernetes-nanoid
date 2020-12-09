### environment set-up ###
source("environmentSetUp.R")

### custom function ###
rna.to.dna = function(h5,suffix = ""){
	  fastq = strsplit(h5,split = "\n")
  read.name = strsplit(fastq[[1]][1],split = " ")
    read.name[[1]][1] = paste0(read.name[[1]][1],suffix)
    fastq[[1]][1] = paste(c(unlist(read.name),""),collapse = " ")
      fastq[[1]][2] = gsub("U","T",fastq[[1]][2])
      paste(c(unlist(fastq),""),collapse = "\n")
}

### create fast5.files vector ###
fast5.files = list.files(FAST5paths,pattern = ".fast5$",full.names = TRUE,recursive = TRUE)
save(fast5.files,file=file.path("data","fast5.files.RData"))

## for any multi-reads fast5 we create a fastq

registerDoParallel(cores = mc.cores)
foreach(fast5.file=fast5.files) %dopar% {
	print(fast5.file)
	## extraction of fast5 reads names
	fast5ls <- h5ls(fast5.file)
	fast5reads <- unique(fast5ls$name[grep("read_",fast5ls$name)])
	## extraction of the Fastq slot in the fastq file
	# sink(paste0(strsplit(fast5.file,split='[.]')[[1]][[1]],".fastq"))
		
	nameTmp <- strsplit(fast5.file,split='FAST5')
	sink(paste0(nameTmp[[1]][[1]],"FASTQ_DNA",strsplit(nameTmp[[1]][[2]],".fast5")[[1]],".fastq"))
			
	for(fast5read in fast5reads)
	{
		cat(rna.to.dna(h5read(fast5.file,paste0(fast5read,"/Analyses/Basecall_1D_001/BaseCalled_template/Fastq")),suffix = ""))
	}
	sink()
}

