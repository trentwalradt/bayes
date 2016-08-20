####### This is doing 1 segment at a time!!!########

library(data.table)
library(pscl)
#Notes on how data file was produced are located at: /data/research/mski_lab/projects/Chantal/blog.R

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

load("/gpfs/commons/groups/imielinski_lab/home/twalradt/cov_seg_input.RData", verbose = TRUE)
data <- data.table(output3)
names(data) <- c("data", "segment")
data <- na.omit(data)
data <- data[is.finite(data),]
data = data[1:1000,]

###SLOIDY: Global sloidy being used.  Is that right???
sloidy = mean(data$data)

#Choose values for purity, ploidy, and variance you want to search. Make sure doing the CN thing right.  Should this be set here and then only modified to offset and spacing at the end when finally graphing????
combo = data.table(expand.grid(ploidy = seq(0.5, 5, 0.1), purity = seq(0.0, 1.0, 0.01), variance = c(0.1,0.2), CN = 0:10))

pp_model <- function(num.segs = 1, data){

    segs = ceiling(length(unique(data$segment)) / num.segs)
    my.list<-vector("list",segs)
    for(i in 1:segs){
        x = data[segment %in% ((i * num.segs) - num.segs + 1):(i * num.segs)]
        input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (sqrt(variance), nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
        input[, CN_prior:= dnorm(CN, ploidy, 1)]
        input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
        input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD) * CN_prior]
        input<-input[,.(logsum = sum(log(prob))),by=.(purity, ploidy, SD, CN, segment)]
        input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
        my.list[[i]]<-input
        print(i)
    }
    return(my.list)
}

input<-rbindlist(pp_model(num.segs=1,data=data))
input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]

##Maximum a posteriori for variance
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
    

#Get purity, ploidy and prob for largest value of SD in MAP
hm <- out[SD == max(MAP$MAP),.(purity,ploidy,prob,lse)]
hm <- acast(hm,ploidy~purity,value.var="prob")


###############################################################################################################################
###### Holding purity, ploidy and SD constant, cycle through all CN and segment ###############

library(data.table)
library(pscl)

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

#Notes on how data file was produced are located at: /data/research/mski_lab/projects/Chantal/blog.R
# This is the first time using data from whole sample, not just chr 1
load("/gpfs/commons/groups/imielinski_lab/home/twalradt/whole_sample_cov_seg.RData", verbose = TRUE)
data <- final_data
data <- na.omit(data)
data <- data[is.finite(data),]
x = data[1:1000,]

###SLOIDY: Global sloidy being used.  Is that right???
sloidy = mean(data$data)

#Choose values for purity, ploidy, and variance you want to search.
combo_whole = data.table(expand.grid(ploidy = seq(0.5, 5, 0.1), purity = seq(0.0, 1.0, 0.01), variance = c(0.1,0.2), CN = 0:10))
setkeyv(combo_whole,c("ploidy", "purity", "variance"))

num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,variance)]))
ppv <- unique(combo_whole[,list(ploidy, purity,variance)])


my.list<-vector("list",num.ppv)
for(i in 1:num.ppv){
    combo = merge(combo_whole, ppv[i,])
    input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (sqrt(variance), nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
    input[, CN_prior:= dnorm(CN, ploidy, 1)]
    input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
    input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD)* CN_prior]
    input<-input[,.(logsum = sum(log(prob))),by=.(purity, ploidy, SD, CN, segment)]
    input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
    input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
    my.list[[i]]<-input
    print(i)
}


input <- rbindlist(my.list)

        
setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]
out.v2 <- out
##Maximum a posteriori for variance
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
MAP.v2<-MAP    

#Get purity, ploidy and prob for largest value of SD in MAP
hm <- out[SD == max(MAP$MAP),.(purity,ploidy,prob,lse)]
hm <- acast(hm,ploidy~purity,value.var="prob")


###############################################################################################################################
###### Holding purity, ploidy and SD constant, cycle through all CN and segment ###############
###### Now with LOOPING! ##################

library(data.table)
library(pscl)

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

#Notes on how data file was produced are located at: /data/research/mski_lab/projects/Chantal/blog.R
# This is the first time using data from whole sample, not just chr 1
load("/gpfs/commons/groups/imielinski_lab/home/twalradt/whole_sample_cov_seg.RData", verbose = TRUE)
data <- final_data
data <- na.omit(data)
data <- data[is.finite(data),]
#x = data[1:1000,]
#Take every 100th element
x = data[1:100==1]

###SLOIDY: Global sloidy being used.  Is that right???
sloidy = mean(data$data)

#Choose values for purity, ploidy, and variance you want to search.
combo_whole = data.table(expand.grid(ploidy = seq(0.5, 5, 0.1), purity = seq(0.0, 1.0, 0.01), variance = c(0.1,0.2), CN = 0:10))
setkeyv(combo_whole,c("ploidy", "purity", "variance"))

num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,variance)]))
ppv <- unique(combo_whole[,list(ploidy, purity,variance)])


bin.size = 10
num.ppv = ceiling(num.ppv/bin.size)


my.list<-vector("list",num.ppv)
system.time(
for(i in 1:num.ppv){
    combo = merge(combo_whole, ppv[((i * bin.size) - bin.size + 1) : (i * bin.size),])
    input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (sqrt(variance), nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
    input[, CN_prior:= dnorm(CN, ploidy, 1)]
    input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
    input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD)* CN_prior]
    input<-input[,.(logsum = sum(log(prob))),by=.(purity, ploidy, SD, CN, segment)]
    input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
    input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
    my.list[[i]]<-input
    print(i)
}
)

input <- rbindlist(my.list)

        
setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]
out.v2 <- out
##Maximum a posteriori for variance. Note, MAP value is in log form
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
MAP.v2<-MAP    

#Get purity, ploidy and prob for largest value of SD in MAP
hm <- out[SD == max(MAP$MAP),.(purity,ploidy,prob,lse)]
hm <- acast(hm,ploidy~purity,value.var="prob")



#################################################### Generate input data #################################################

library(skitools)
library(DNAcopy)

#NOTE, this is some sort of stocastic process, use set.seed() if you want consistent results
.segment = function(tcov){
    ix = which(!is.na(tcov$ratio))
    cat('sending ', length(ix), ' segments\n')
    cna = CNA(log(tcov$ratio[ix]), as.character(seqnames(tcov))[ix], start(tcov)[ix], data.type = 'logratio')
    gc()
    cat('finished making cna\n')
    Seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = T)
    cat('finished segmenting\n')
    ## out = seg2gr(print(seg), new.sl) ## remove seqlengths that have not been segmented
    ## out = gr.fix(out, new.sl, drop = T)
    return(Seg)
}

#"/data/research/mski_lab/projects/Chantal/Flow/CBS/29331-B_29331-T/cov.rds"

set.bin.size <- function(n = 500, data.path = "/data/research/mski_lab/projects/Chantal/Flow/CBS/29407_B_29407_T1/cov.rds"){
    CBS_cov<-readRDS(data.path)
    raw <- cbind(CBS_cov$ratio,as.data.table(ranges(CBS_cov))$start,as.data.table(seqnames(CBS_cov)))
    raw[V1 ==Inf,V1:=NA]
    new <- raw[,.(ratio = mean(V1,na.rm=T),start = min(V2)),by=.(value,(seq(nrow(raw)) - 1) %/% n)]
    new <- new[!is.na(ratio)]
    new <- new[ratio > 0]
    seq.names <- as.data.table(table(new$value))
    #Create new grange table which is required for .segment function
    grange_new <- GRanges(seqnames = Rle(seq.names$V1, seq.names$N), ranges = IRanges(new$start,end=new$start+n*200-1,width=NULL,names=NULL), ratio = new$ratio)
    new.seg<-.segment(grange_new)
    final <- cbind(new, rep(seq(1, length(new.seg$output$num.mark), 1), new.seg$output$num.mark))
    names(final) <- c("chr","seq","data","start","segment")
    return(final)
}


###############################################################################################################################
###### Run this version do determine the SD ###################

library(data.table)
library(pscl)

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

set.seed(1)
y = set.bin.size()
# Now y has parts 1 and 2 y[[1]] and y[[2]] which are the data.table and grange object
x <- y
x <- x[sample(nrow(x),nrow(x)/10)]
x <- x[sample(nrow(x),nrow(x)/100)]

###SLOIDY: Global sloidy being used.  Is that right???
sloidy = mean(x$data)

#Choose values for purity, ploidy, and variance you want to search.
combo_whole = data.table(expand.grid(ploidy = seq(0.5, 5, 0.5), purity = seq(0.1, 1.0, 0.1), SD = seq(sd(x$data/100),sd(x$data),sd(x$data/100)), CN = 0:10))

setkeyv(combo_whole,c("ploidy", "purity", "SD"))

num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
ppv <- unique(combo_whole[,list(ploidy, purity, SD)])


bin.size = 50
num.ppv = ceiling(num.ppv/bin.size)

#Maybe try to only look at copy states <=3 the ploidy.   E.g. go up to CN=6 for ploidy=2

pp.list<-vector("list",num.ppv)
system.time(
for(i in 1:num.ppv){
    combo = merge(combo_whole, ppv[((i * bin.size) - bin.size + 1) : (i * bin.size),])
    input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
    input[, CN_prior:= dnorm(CN, ploidy, 1)]
    input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
    #Added log=T to dnorm. Also now adding low of CN_prior
    input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
    #Removed the log from inside the sume here
    input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
    input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
    input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
    pp.list[[i]]<-input
    print(i)
}
)

input <- rbindlist(pp.list)


setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]
##Maximum a posteriori for variance. Note, MAP value is in log form
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
SD.MAP = MAP[MAP == max(MAP)]$SD
#######Sample 1
## SD.MAP is  0.0452515976817749 for whole sample
## SD.MAP is 0.04275849 for sample/10
## SD.MAP is 0.05260518 for sample/100
#######Sample /data/research/mski_lab/projects/Chantal/Flow/CBS/29407_B_29407_T1/cov.rds
## SD.MAP is 0.06243132 for whole sample
## SD.MAP is 0.06421183 for sample/10
## SD.MAP is 0.06128031 for sample/100
## SD.MAP is 0.05838234 for sample w/ n = 100 sample/100
## SD.MAP is 0.06633568 for sample w/ n = 50 sample/100
## SD.MAP is 0.07052021 for sample w/ n = 50 sample/10
## SD.MAP is 0.06849585 for sample w/ n = 50 whole
## SD.MAP is 0.06276308 for sample w/ n = 100 sample/10
## SD.MAP is 0.06457467 for sample w/ n = 100 whole

#############################################################
############ Run this w/ max SD to get purity and ploidy#####


#Choose values for purity, ploidy, and variance you want to search.
combo_whole = data.table(expand.grid(ploidy = seq(0.55, 5.5, 0.05), purity = seq(0.01, 1.0, 0.01), SD = SD.MAP, CN = 0:10))
setkeyv(combo_whole,c("ploidy", "purity", "SD"))

num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
ppv <- unique(combo_whole[,list(ploidy, purity, SD)])


bin.size = 50
num.ppv = ceiling(num.ppv/bin.size)

#Maybe try to only look at copy states <=3 the ploidy.   E.g. go up to CN=6 for ploidy=2

CN.list<-vector("list",num.ppv)
pp.list<-vector("list",num.ppv)
system.time(
for(i in 1:num.ppv){
    combo = merge(combo_whole, ppv[((i * bin.size) - bin.size + 1) : (i * bin.size),])
    input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
    input[, CN_prior:= dnorm(CN, ploidy, 1)]
    input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
    #Added log=T to dnorm. Also now adding low of CN_prior
    input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
    #removed the log from inside the sume here
    input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
    #Maybe store the data at this stage to get the prob for every CN state for every segment
    CN.list[[i]] <- dcast(input,purity + ploidy + segment ~ CN, value.var = 'logsum')
    input <- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
    input <- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
    pp.list[[i]] <- input
    print(i)
}
)

input <- rbindlist(pp.list)
CN.matrix <- rbindlist(CN.list)


input.pp.trial <- input
save(input.pp.trial, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial_whole.RData")
input.pp.trial.10 <- input
save(input.pp.trial.10, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial_10.RData")
input.pp.trial.100 <- input
save(input.pp.trial.100, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial_whole.RData")
#### /data/research/mski_lab/projects/Chantal/Flow/CBS/29407_B_29407_T1/cov.rds
input.pp.trial2 <- input
save(input.pp.trial2, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_whole.RData")
input.pp.trial2.10 <- input
save(input.pp.trial2.10, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_10.RData")
input.pp.trial2.100 <- input
save(input.pp.trial2.100, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_100.RData")
input.pp.trial2.n100.whole <- input
save(input.pp.trial2.n100.whole, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n100_whole.RData")
input.pp.trial2.n100.100 <- input
save(input.pp.trial2.n100.100, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n100_100.RData")
input.pp.trial2.n50.100 <- input
save(input.pp.trial2.n50.100, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n50_100.RData")
input.pp.trial2.n50.10 <- input
save(input.pp.trial2.n50.10, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n50_10.RData")
input.pp.trial2.n50.whole <- input
save(input.pp.trial2.n50.whole, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n50_whole.RData")

input.pp.trial2.n100.10 <- input
save(input.pp.trial2.n100.10, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial2_n100_10.RData")




input<-input.pp.trial2.100

setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]

#Make heat map
hm <- acast(out,ploidy~purity,value.var="prob")
library(d3heatmap)
wij(d3heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE),filename="sample2_n50_whole.html")

#Try something like this:
wij(d3heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE,main="Sample 2 100", xlab = "Purity", ylab = "Ploidy"),filename="sample2_100_v2.html")

### Select purity ploidy solutions w/ probability greater than some cutoff value and subset CN.matrix
prob.cutoff = 0.01
out<-out[prob > prob.cutoff]
setkeyv(out,c("purity", "ploidy"))
setkeyv(CN.matrix,c("purity", "ploidy"))
high.prob.CN.matrix <- merge(out, CN.matrix)

#SUPER hacky way of doing this, try again later
CN<-high.prob.CN.matrix[,(exp(.SD)/rowSums(exp(.SD))), .SDcols = 8:18]
CN.output <- cbind(high.prob.CN.matrix[,.(purity, ploidy, segment, prob)], CN)


library(MASS)

# Create list before running loop
# Not vectorized, definitely a faster way to do this
segments <- unique(CN.output$segment)
hist.list <- vector("list", length(segments))
alpha.list <- vector("list", length(segments))
beta.list <- vector("list", length(segments))
for(z in 1:length(segments)){

    combo_spec <- CN.output[segment==segments[z]]

## Create 2 column matrix now: P(alpha)P(k) and alpha(K*/K)

    num.CN = 10
    my.list<-vector("list",num.CN)
    for(i in 1:10){
        prob <- rep(combo_spec$prob * combo_spec[[i + 5]], each = i)
        loc <- c(t(do.call(base:::'%o%',list(combo_spec$purity, ((1:i) / i)))))
        out <- cbind(prob, loc)
        my.list[[i]] <- out
    }

    hist <- data.table(do.call(rbind, my.list))
    hist[,prob2 := prob/sum(prob)]
    hist.matrix = hist[, sample(loc, 10000, prob = prob2, replace = TRUE)]
    beta = fitdistr(hist.matrix, dbeta, start = list(shape1 = 1, shape2 = 20))
    hist.list[[z]] <- hist.matrix
    alpha.list[[z]] <- beta$estimate[1]
    beta.list[[z]] <- beta$estimate[2]
    print(z)
}

hist.out <- data.table(do.call(cbind, hist.list))
colnames(hist.out) <- as.character(segments)
alpha.out <- data.table(cbind(alpha.list))
beta.out <- data.table(cbind(beta.list))
alpha.beta <- cbind( segments, alpha.out, beta.out)
setkey(alpha.beta, segments)

setkey(x, segment)

# Do this if you want ALL segments annotated with alpha and beta -> this produces NULL for segments that weren't subsampled if subsapmling was performed
x <- alpha.beta[x]

#Do this if you want only segments with a non NULL alpha and beta
x <- x[alpha.beta]

#Make sure when calling n here you referene the relevatn opt$ in run.R
seq.names <- as.data.table(table(x$chr))
    #Create new grange table which is required for .segment function
grange.out <- GRanges(seqnames = Rle(seq.names$V1, seq.names$N), ranges = IRanges(x$start,end=x$start+n*200-1,width=NULL,names=NULL), ratio = x$data, alpha = as.numeric(x$alpha.list), beta = as.numeric(x$beta.list), segment = as.numeric(x$segment))


##### To do
# Create granges object with segments and local alpha and beta
# Make matrix with a column for histogram values for each segment
# Comment code so that is is all easy to read and understand







##########################################################################
##########################################################################
## Wrapped up in a package

library(skitools)
library(DNAcopy)
library(data.table)
library(pscl)

#################################################### Generate input data #################################################
p.p.model <- function(bin.red = 500, data.path = "/data/research/mski_lab/projects/Chantal/Flow/CBS/29407_B_29407_T1/cov.rds", sub.sample.fold = 10, bin.num = 50){
#output.name as another argument?

    .segment = function(tcov){
        ix = which(!is.na(tcov$ratio))
        cat('sending ', length(ix), ' segments\n')
        cna = CNA(log(tcov$ratio[ix]), as.character(seqnames(tcov))[ix], start(tcov)[ix], data.type = 'logratio')
        gc()
        cat('finished making cna\n')
        Seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = T)
        cat('finished segmenting\n')
        ## out = seg2gr(print(seg), new.sl) ## remove seqlengths that have not been segmented
        ## out = gr.fix(out, new.sl, drop = T)
        return(Seg)
    }

    #"/data/research/mski_lab/projects/Chantal/Flow/CBS/29331-B_29331-T/cov.rds"

    set.bin.size <- function(n = bin.red, data_path = data.path){
        CBS_cov<-readRDS(data_path)
        raw <- cbind(CBS_cov$ratio,as.data.table(ranges(CBS_cov))$start,as.data.table(seqnames(CBS_cov)))
        raw[V1 ==Inf,V1:=NA]
        new <- raw[,.(ratio = mean(V1,na.rm=T),start = min(V2)),by=.(value,(seq(nrow(raw)) - 1) %/% n)]
        new <- new[!is.na(ratio)]
        new <- new[ratio > 0]
        seq.names <- as.data.table(table(new$value))
                                        #Create new grange table which is required for .segment function
        grange_new <- GRanges(seqnames = Rle(seq.names$V1, seq.names$N), ranges = IRanges(new$start,end=new$start+n*200-1,width=NULL,names=NULL), ratio = new$ratio)
        new.seg<-.segment(grange_new)
        final <- cbind(new, rep(seq(1, length(new.seg$output$num.mark), 1), new.seg$output$num.mark))
        names(final) <- c("chr","seq","data","start","segment")
        return(final)
    }

    set.seed(1)
    y = set.bin.size()

    sub.sample <- function(num = sub.sample.fold, z = y){
        x <- z[sample(nrow(z),nrow(z)/num)]
        return(x)
    }

    log.sum.exp<- function(x){
        offset <- max(x)
        log(sum(exp(x - offset))) + offset
    }

    sub.sample()


###############################################################################################################################
###### Run this version do determine the SD ###################


###SLOIDY: Global sloidy being used.  Is that right???
    sloidy = mean(x$data)

    #Choose values for purity, ploidy, and variance you want to search.
    combo_whole = data.table(expand.grid(ploidy = seq(0.5, 5, 0.5), purity = seq(0.1, 1.0, 0.1), SD = seq(sd(x$data/100),sd(x$data),sd(x$data/100)), CN = 0:10))

    setkeyv(combo_whole,c("ploidy", "purity", "SD"))

    num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
    ppv <- unique(combo_whole[,list(ploidy, purity, SD)])



    num.ppv = ceiling(num.ppv/bin.num)

    #Maybe try to only look at copy states <=3 the ploidy.   E.g. go up to CN=6 for ploidy=2

    my.list<-vector("list",num.ppv)
    system.time(
        for(i in 1:num.ppv){
            combo = merge(combo_whole, ppv[((i * bin.num) - bin.num + 1) : (i * bin.num),])
            input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
            input[, CN_prior:= dnorm(CN, ploidy, 1)]
            input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
                                        #Added log=T to dnorm. Also now adding low of CN_prior
            input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
                                        #Removed the log from inside the sume here
            input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
            input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
            input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
            my.list[[i]]<-input
            print(i)
        }
        )

    input <- rbindlist(my.list)

    setkey(input,SD)
    int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
    out = merge(input, int)[, prob := exp(final - lse)]
    ##Maximum a posteriori for variance. Note, MAP value is in log form
    MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
    SD.MAP = MAP[MAP == max(MAP)]$SD

#############################################################
############ Run this w/ max SD to get purity and ploidy#####

    #Choose values for purity, ploidy, and variance you want to search.
    combo_whole = data.table(expand.grid(ploidy = seq(0.55, 5.5, 0.05), purity = seq(0.01, 1.0, 0.01), SD = SD.MAP, CN = 0:10))
    setkeyv(combo_whole,c("ploidy", "purity", "SD"))

    num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
    ppv <- unique(combo_whole[,list(ploidy, purity, SD)])

    num.ppv = ceiling(num.ppv/bin.num)

    #Maybe try to only look at copy states <=3 the ploidy.   E.g. go up to CN=6 for ploidy=2

    my.list<-vector("list",num.ppv)
    system.time(
        for(i in 1:num.ppv){
            combo = merge(combo_whole, ppv[((i * bin.num) - bin.num + 1) : (i * bin.num),])
            input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
            input[, CN_prior:= dnorm(CN, ploidy, 1)]
            input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
            input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
            input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
            input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
            input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
            my.list[[i]]<-input
            print(i)
        }
        )

    input <- rbindlist(my.list)

    #input.pp.trial <- input
    #save(input.pp.trial, file = "/gpfs/commons/groups/imielinski_lab/home/twalradt/input_pp_trial_whole.RData")

    setkey(input,SD)
    int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
    out = merge(input, int)[, prob := exp(final - lse)]

    return(list(SD.MAP,out))

}
#Make heat map
hm <- acast(out,ploidy~purity,value.var="prob")

wij(d3heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE),filename="sample2_100.html")

#Try something like this:
wij(d3heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE,main="Sample 2 100", xlab = "Purity", ylab = "Ploidy"),filename="sample2_100_v2.html")



