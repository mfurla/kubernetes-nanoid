### environment set-up ###
source("environmentSetUp.R")

### load objects ###
objects = list.files(file.path("Supplemental_Code_S7","BasicObjects"))

# MF: I changed source with load because of the following error occurred: "Error in eval(expr, envir, enclos) : object 'RDX2' not found"
for(object in objects){load(file.path("Supplemental_Code_S7/BasicObjects",object))}

### build.trace.list ###
build.trace.list = function(read.name){
  try({
    fast5.file = read.to.fast5.name[read.name]
    read.sequence = rev(strsplit(read.sequence.list[[read.name]],split = "")[[1]])
    trace = h5read(fast5.file,paste0("read_",read.name,"/Analyses/",slotName,"/BaseCalled_template/Trace"))
    rownames(trace) = c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T")
    move = h5read(fast5.file,paste0("read_",read.name,"/Analyses/",slotName,"/BaseCalled_template/Move"))
    move.rle = Rle(move)
    event.repeats = move[move == 1]
    event.repeats[cumsum(runLength(move.rle)[runValue(move.rle) == 1])] = runLength(move.rle)[runValue(move.rle) == 0] + 1
    if (move[length(move)] == 1){event.repeats[length(event.repeats)] = 1}
    colnames(trace) = rep(read.sequence,times = event.repeats)
    trace.move = trace[,move > 0]
    trace.no.move = trace[,move == 0]
    
    results = c(read.name,
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "A",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "C",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "G",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) == "T",drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                
                tryCatch({apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                tryCatch({apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,mean)},error = function(x){0}),
                
                tryCatch({as.vector(t(apply(trace[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)}),
                tryCatch({as.vector(t(apply(trace.no.move[c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),colnames(trace.no.move) %in% c("A","C","G"),drop = FALSE],1,function(x){quantile(x,seq(0,1,length.out = 10))})))},error = function(x){rep(0,10)})
    )
    names(results) = NULL
    return(results)
  },silent = TRUE)
}

col.names = c("read.name",
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_A_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_C_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_G_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_T_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_move_ACG_"),seq(0,1,length.out = 10),paste0)),
              paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_mean"),
              as.vector(outer(paste0(c("flip_A","flip_C","flip_G","flip_T","flop_A","flop_C","flop_G","flop_T"),"_no.move_ACG_"),seq(0,1,length.out = 10),paste0)))

### probability of traces ###
bam = get(load(file.path("data","bam.RData")))
bamNames <- names(bam)
rm(bam)
gc()

read.to.fast5.name = get(load(file.path("data","read.to.fast5.name.RData")))

# MF: the splitvector function is missing
# index.list = splitvector(1:length(bamNames),100)
index.list = split(1:length(bamNames), ceiling(1:length(bamNames)/100))

if(file.exists(where="data/trace.list.llo"))
{
   start = getListLength("data/trace.list.llo")/100 + 1
   print(paste0("First index: ",start))
}else{
  trace.list = list()
  saveList(object = trace.list, file = "data/trace.list.llo", append = FALSE, compress = FALSE)

  start = 1
}

registerDoParallel(cores = mc.cores)
for (j in start:length(index.list)){
  if(j%%50==0){print(j);system("grep MemFree /proc/meminfo")}

  read.sequence.list <- readList("data/read.sequence.list.llo",index.list[[j]])
  names(read.sequence.list) <- bamNames[index.list[[j]]]

  trace.list.tmp = foreach(n = bamNames[index.list[[j]]],.noexport = setdiff(ls(),c("bam","read.to.fast5.name","read.sequence.list"))) %dopar% build.trace.list(n)
  saveList(object = trace.list.tmp, file = "data/trace.list.llo", append = TRUE, compress = FALSE)
  
  rm(trace.list.tmp)
  rm(read.sequence.list)
  gc()
}

trace.list <- readList("data/trace.list.llo")

traces = t(sapply(trace.list,c))
colnames(traces) = col.names
rownames(traces) = traces[,"read.name"]

traces <- traces[,setdiff(col.names,"read.name")]
storage.mode(traces) <- "numeric"

save(traces,file=file.path("data","traces.RData"))

names(trace.list) <- bamNames
saveList(object = trace.list, file = "data/trace.list.llo", append = FALSE, compress = FALSE)
# system("rm data/trace.list.llo")

