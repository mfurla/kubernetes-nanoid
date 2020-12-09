### environment set-up ###
source("environmentSetUp.R")

### load objects ###
objects = list.files(file.path("Supplemental_Code_S7","BasicObjects"))

# MF: I changed source with load because of the following error occurred: "Error in eval(expr, envir, enclos) : object 'RDX2' not found"
for(object in objects){load(file.path("Supplemental_Code_S7/BasicObjects",object))}

### build.trace.add.on.list ###
build.trace.add.on.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
    trace = h5read(fast5.file,paste0("read_",read.name,"/Analyses/Basecall_1D_001/BaseCalled_template/Trace"))
    rownames(trace) = c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T")
    move = h5read(fast5.file,paste0("read_",read.name,"/Analyses/Basecall_1D_001/BaseCalled_template/Move"))
    move.rle = Rle(move)
    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    colnames(trace) = rep(read.sequence,times = event.repeats)
    
    move.rle.A = Rle(move[colnames(trace) == "A"])
    move.rle.C = Rle(move[colnames(trace) == "C"])
    move.rle.G = Rle(move[colnames(trace) == "G"])
    move.rle.T = Rle(move[colnames(trace) == "T"])
    
    results = c(read.name,
                
                mean(runLength(move.rle)[runValue(move.rle) == 1]),
                quantile(runLength(move.rle)[runValue(move.rle) == 1],seq(0,1,length.out = 100)),
                
                mean(runLength(move.rle)[runValue(move.rle) == 0]),
                quantile(runLength(move.rle)[runValue(move.rle) == 0],seq(0,1,length.out = 100)),
                
                tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.A)[runValue(move.rle.A) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.A)[runValue(move.rle.A) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.C)[runValue(move.rle.C) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.C)[runValue(move.rle.C) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.G)[runValue(move.rle.G) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.G)[runValue(move.rle.G) == 0],seq(0,1,length.out = 100))},error = function(x){0}),
                
                tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 1])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 1],seq(0,1,length.out = 100))},error = function(x){0}),
                tryCatch({mean(runLength(move.rle.T)[runValue(move.rle.T) == 0])},error = function(x){0}),
                tryCatch({quantile(runLength(move.rle.T)[runValue(move.rle.T) == 0],seq(0,1,length.out = 100))},error = function(x){0})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              paste0("1s","_mean"),
              paste0("1s_",1:100),
              paste0("0s","_mean"),
              paste0("0s_",1:100),
              paste0("1s","_A_mean"),
              paste0("1s_A_",1:100),
              paste0("0s","_A_mean"),
              paste0("0s_A_",1:100),
              paste0("1s","_C_mean"),
              paste0("1s_C_",1:100),
              paste0("0s","_C_mean"),
              paste0("0s_C_",1:100),
              paste0("1s","_G_mean"),
              paste0("1s_G_",1:100),
              paste0("0s","_G_mean"),
              paste0("0s_G_",1:100),
              paste0("1s","_T_mean"),
              paste0("1s_T_",1:100),
              paste0("0s","_T_mean"),
              paste0("0s_T_",1:100))

### probability of traces ###
bam = get(load(file.path("data","bam.RData")))
bamNames <- names(bam)
rm(bam)
gc()

read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))
read.sequence.list = get(load(file.path("data","read.sequence.list.RData")))

# MF: the splitvector function is missing
# index.list = splitvector(1:length(bamNames),100)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

trace.add.on.list = list()
saveList(object = trace.add.on.list, file = "data/trace.add.on.list.llo", append = FALSE, compress = FALSE)

registerDoParallel(cores = mc.cores)
for (j in 1:length(index.list)){

  if(j%%50==0){print(j);gc();system("grep MemFree /proc/meminfo")} 
  trace.add.on.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.add.on.list(n)
  saveList(object = trace.add.on.list.tmp, file = "data/trace.add.on.list.llo", append = TRUE, compress = FALSE)
  rm(trace.add.on.list.tmp)
  gc()
}

trace.add.on.list <- readList("data/trace.add.on.list.llo")

traces.add.on = t(sapply(trace.add.on.list,c))
colnames(traces.add.on) = col.names
rownames(traces.add.on) = traces.add.on[,"read.name"]

traces.add.on <- traces.add.on[,setdiff(col.names,"read.name")]
storage.mode(traces.add.on) <- "numeric"

save(traces.add.on,file=file.path("data","traces.add.on.RData"))

names(trace.add.on.list) <- bamNames
saveList(object = trace.add.on.list, file = "data/trace.add.on.list.llo", append = FALSE, compress = FALSE)

# system("rm data/trace.add.on.list.llo")
