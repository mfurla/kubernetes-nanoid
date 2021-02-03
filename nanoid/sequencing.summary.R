### environment set-up ###
source("environmentSetUp.R")

# Bam names
#
load(file=file.path("data","bam.RData"))
bamNames <- names(bam)
rm(bam)
gc()

### create sequencing summary matrices ###
sequencing.files = list.files(strsplit(FAST5paths,"/FAST5")[[1]],pattern = "sequencing_summary",full.names = TRUE,recursive = TRUE)

sequencing.summary <- read.delim(sequencing.files[[1]],sep="\t",row.names = NULL,header = TRUE,stringsAsFactors = FALSE)

if(length(sequencing.files>1))
{
	for(sequencing.file in sequencing.files[-1])
	{
		sequencing.summary <- rbind(sequencing.summary,read.delim(sequencing.file,sep="\t",row.names = NULL,header = TRUE,stringsAsFactors = FALSE))
	}
}

rownames(sequencing.summary) = sequencing.summary[,"read_id"]

sequencing.summary <- sequencing.summary[bamNames,]

save(sequencing.summary,file=file.path("data","sequencing.summary.RData"))