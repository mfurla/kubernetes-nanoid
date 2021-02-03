### environment set-up ###
source("environmentSetUp.R")

### data loading ###
load(file.path("data","read.based.mismatch.identification.RData"))
load(file.path("data","traces.RData"))
load(file.path("data","traces.add.on.RData"))
traces = cbind(traces,traces.add.on[rownames(traces),])
load(file.path("data","raw.signal.RData"))
load(file.path("data","raw.signal.five.mers.RData"))
load(file.path("data","raw.signal.five.mers.add.on.RData"))
raw.signal.five.mers = cbind(raw.signal.five.mers,raw.signal.five.mers.add.on[rownames(raw.signal.five.mers),])
read.based.parameters = cbind(read.based.mismatch.identification,traces[rownames(read.based.mismatch.identification),],raw.signal[rownames(read.based.mismatch.identification),],raw.signal.five.mers[rownames(read.based.mismatch.identification),])
colnames(read.based.parameters) = paste0("P",substr(rep("0000",ncol(read.based.parameters)),1,4-nchar(1:ncol(read.based.parameters))),1:ncol(read.based.parameters),"-",colnames(read.based.parameters))

save(read.based.parameters,file=file.path("data","read.based.parameters.RData"))