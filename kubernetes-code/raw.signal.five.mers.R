### environment set-up ###
source("environmentSetUp.R")

### load objects ###
objects = list.files(file.path("Supplemental_Code_S7","BasicObjects"))

# MF: I changed source with load because of the following error occurred: "Error in eval(expr, envir, enclos) : object 'RDX2' not found"
for(object in objects){load(file.path("Supplemental_Code_S7/BasicObjects",object))}

### build.raw.signal.five.mers.list ###
build.raw.signal.five.mers.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    
    # read.number = strsplit(strsplit(fast5.file,split = "read_")[[1]][2],split = "_")[[1]][1]
    # if (is.na(read.number)){read.number = strsplit(strsplit(fast5.file,split = "read")[[1]][3],split = "_")[[1]][1]}
    
    # if (is.na(read.number)){raw.signal = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Raw/Signal"))} else {raw.signal = h5read(file.path(prewd,fast5.file),paste0("/Raw/Reads/Read_",read.number,"/Signal"))}
    raw.signal = h5read(fast5.file,paste0("read_",read.name,"/Raw/Signal"))

    # if (is.na(read.number)){move = h5read(file.path(prewd,fast5.file),paste0("read_",strsplit(read.name,split = "\\.")[[1]][1],"/Analyses/Basecall_1D_000/BaseCalled_template/Move"))} else {move = h5read(file.path(prewd,fast5.file),"/Analyses/Basecall_1D_000/BaseCalled_template/Move")}
    move = h5read(fast5.file,paste0("read_",read.name,"/Analyses/Basecall_1D_001/BaseCalled_template/Move"))
    move.rle = Rle(move)

    read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
    
    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    
    stride = 10
    
    num_events_template = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events_template"]
    num_events = sequencing.summary[strsplit(read.name,split = "\\.")[[1]][1],"num_events"]
    
    anchor = (num_events - num_events_template)*stride
    anchored.raw.signal = raw.signal[anchor:length(raw.signal)]
    
    read.sequence.five.mers = as.character(five.mers.minus[as.numeric(stats::filter(as.numeric(c("A" = 1,"C" = 2,"G" = 3,"T" = 4,"-" = 5)[read.sequence]),5^(0:8)[1:5],sides = 1)-(sum(5^(0:8)[1:5])-1))[-c(1:4)]])
    read.sequence.five.mers = c(NA,NA,read.sequence.five.mers,NA,NA)
    
    events.to.raw = cut(1:length(anchored.raw.signal),(0:num_events_template)*stride,right = FALSE,labels = FALSE)
    event.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = events.to.raw,identity),function(x){anchored.raw.signal[x]})
    raw.means = sapply(event.raw.signal,mean)
    
    cumsum.event.repeats = cumsum(event.repeats)
    
    reevent = cbind(c(1,cumsum.event.repeats[-length(cumsum.event.repeats)] + 1),cumsum.event.repeats)
    colnames(reevent) = c("start","end")

    reevents.to.raw = cut(1:length(anchored.raw.signal),c(((0:num_events_template)*stride)[reevent[,"start"]],num_events_template*stride+1),right = FALSE,labels = FALSE)
    reevent.raw.signal = lapply(tapply(1:length(anchored.raw.signal),INDEX = reevents.to.raw,identity),function(x){anchored.raw.signal[x]})
    reevent.raw.means = sapply(reevent.raw.signal,mean)

    events = cbind("move" = move,
                   "read.sequence" = rep(read.sequence,times = event.repeats),
                   "read.sequence.five.mers" = rep(read.sequence.five.mers,times = event.repeats))
    
    means = cbind("raw.means" = raw.means,
                   "reevent.raw.means" = rep(reevent.raw.means,times = event.repeats))
    
    events.move = events[move == 1,]
    means.move = means[move == 1,]
    
    events.no.move = events[move == 0,]
    means.no.move = means[move == 0,]
    
    results = c(read.name,
                sapply(five.mers,function(x){apply(means.move[events.move[,"read.sequence.five.mers"] %in% x,"reevent.raw.means",drop = FALSE],2,mean)})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              mkAllStrings(c("A","C","G","T"),5))

z.score.rescaling = function(x,y,z){((z - mean(y,na.rm = TRUE))*(sd(x,na.rm = TRUE)/sd(y,na.rm = TRUE))) + mean(x,na.rm = TRUE)}

non.T.containing = rownames(five.mer.RNA.pore.model) %in% five.mers.non.T.containing

### raw signal of kmers ###
bam = get(load(file.path("data","bam.RData")))
read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))

sequencing.summary = get(load(file.path("data","sequencing.summary.RData")))

bamNames <- names(bam)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

raw.signal.five.mers.list = list()
saveList(object = raw.signal.five.mers.list, file = "data/raw.signal.five.mers.list.llo", append = FALSE, compress = FALSE)

registerDoParallel(cores = mc.cores)

for (j in 1:length(index.list)){
  if(j%%25==0){print(j);system("grep MemFree /proc/meminfo")}

  read.sequence.list = readList("data/read.sequence.list.llo",index=index.list[[j]])
  names(read.sequence.list) <- bamNames[index.list[[j]]]

  raw.signal.five.mers.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list","sequencing.summary"))) %dopar% build.raw.signal.five.mers.list(n)
  saveList(object = raw.signal.five.mers.tmp, file = "data/raw.signal.five.mers.list.llo", append = TRUE, compress = FALSE)

  rm(raw.signal.five.mers.tmp)
  rm(read.sequence.list)
  gc()
}
# MF: Warning message:
# In event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) ==  : number of items to replace is not a multiple of replacement length

raw.signal.five.mers.list <- readList("data/raw.signal.five.mers.list.llo")
names(raw.signal.five.mers.list) <- bamNames

raw.signal.five.mers = t(sapply(raw.signal.five.mers.list,c))
colnames(raw.signal.five.mers) = col.names
rownames(raw.signal.five.mers) = raw.signal.five.mers[,"read.name"]

raw.signal.five.mers <- raw.signal.five.mers[,setdiff(col.names,"read.name")]
storage.mode(raw.signal.five.mers) <- "numeric"

save(raw.signal.five.mers,file=file.path("data","raw.signal.five.mers.RData"))

saveList(object = raw.signal.five.mers.list, file = "data/raw.signal.five.mers.list.llo", append = FALSE, compress = FALSE)
# system("rm data/raw.signal.five.mers.list.llo")