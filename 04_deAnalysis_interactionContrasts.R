# File: 04_deAnalysis_interactionContrasts.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: selecting DE genes on interaction contrasts
# Date: 01/5/2020

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# load the count matrix
dfData = read.csv('dataExternal/Partek_rna2_Quantify_to_annotation_model_(Partek_E_M)_Gene_counts.csv',
                  sep='\t', header=T, stringsAsFactors = F)
mCounts = as.matrix(dfData[,-c(1:5)])
rownames(mCounts) = dfData$Gene.Symbol
dim(mCounts)
## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 45')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

# sanity check
table(dfSample$title %in% colnames(mCounts))
mCounts = mCounts[,dfSample$title]
identical(dfSample$title, (colnames(mCounts)))

mData = mCounts
dim(mData)
colnames(mData) = dfSample$id

## merge the technical replicates i.e. lanes 
fReplicates = strsplit(dfSample$description, ';')
fReplicates = factor(sapply(fReplicates, function(x) return(x[7])))
nlevels(fReplicates)
fLanes = strsplit(dfSample$description, ';')
fLanes = factor(sapply(fLanes, function(x) return(x[2])))
table(fReplicates, fLanes)

dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))
dim(dfSample.2)
dfSample.2 = droplevels.data.frame(dfSample.2)

## normalise the data
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 15480  1739 
mData = mData[!(i< 3),]
dim(mData)
# [1] 15480    16

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))

## delete sample section after testing
mData.norm = round(mData.norm, 0)

# set.seed(123)
# i = sample(1:nrow(mData.norm), 10, replace = F)
# dfData = data.frame(t(mData.norm[i,]))

dfData = data.frame(t(mData.norm))
dim(dfData)
dfData = stack(dfData)

## create covariates for modelling
str(dfSample.2)
f = paste(dfSample.2$group1, dfSample.2$group2, sep='_')
dfData$fTreatment = factor(f)
fDesc = strsplit(dfSample.2$description, ';')
fDesc = (sapply(fDesc, function(x) return(x[11])))
fNewDay = rep('Low', times=length(fDesc))
fNewDay[fDesc %in% c('Day_5', 'Day_6', 'Day_7', 'Day_8')] = 'High'
fNewDay = factor(fNewDay, levels = c('Low', 'High'))

dfData$fDay = fNewDay
dfData = droplevels.data.frame(dfData)

dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fDay:dfData$ind)
str(dfData)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load('results/fit.stan.nb_22Apr.rds')
print(fit.stan, c('sigmaRan1'), digits=3)

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

getDifferenceVector = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  return(d)
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef.1))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base1='Sham_WT', deflection1='MI_WT', base2='Sham_Tg', deflection2='MI_Tg') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif1 = getDifferenceVector(ivData = mCoef[,c[deflection1]], ivBaseline = mCoef[,c[base1]])
  dif2 = getDifferenceVector(ivData = mCoef[,c[deflection2]], ivBaseline = mCoef[,c[base2]])
  dif = getDifference(ivData = dif1, ivBaseline = dif2)
  r = data.frame(ind= as.character(d$ind[c[base1]]), coef.base1=mean(mCoef[,c[base1]]), 
                 coef.deflection1=mean(mCoef[,c[deflection1]]),
                 coef.base2=mean(mCoef[,c[base2]]), coef.deflection2=mean(mCoef[,c[deflection2]]),
                 zscore=dif$z, pvalue=dif$p)
  r$difference = (r$coef.deflection1 - r$coef.base1) - (r$coef.deflection2 - r$coef.base2)
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
# library(org.Hs.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))
# df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENSEMBL')
# i = match(dfResults$ind, df$ENTREZID)
# df = df[i,]
dfResults$SYMBOL = as.character(dfResults$ind)
# identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'MI_WT vs Sham_WT')#, fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'MI_WT vs Sham_WT', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'MI_Tg vs Sham_Tg')
table(dfResults$adj.P.Val < 0.1)
table(dfResults$pvalue < 0.01)
## save the results 
write.csv(dfResults, file='results/DEAnalysisMIWtVsShamWt.xls')

######### do a comparison with deseq2
str(dfSample.2)
f = paste(dfSample.2$group1, dfSample.2$group2, sep='_')
fTreatment = factor(f)
fDesc = strsplit(dfSample.2$description, ';')
fDesc = (sapply(fDesc, function(x) return(x[11])))
fNewDay = rep('Low', times=length(fDesc))
fNewDay[fDesc %in% c('Day_5', 'Day_6', 'Day_7', 'Day_8')] = 'High'
fNewDay = factor(fNewDay, levels = c('Low', 'High'))

dfDesign = data.frame(Treatment = fTreatment, Patient=fNewDay,
                      row.names=colnames(mData))
mData2 = apply(mData, 2, as.integer)
rownames(mData2) = rownames(mData)

oDseq = DESeqDataSetFromMatrix(mData2, dfDesign, design = ~ Treatment)# + Patient)
oDseq = DESeq(oDseq)

plotDispEsts(oDseq)
oRes = results(oDseq, contrast = c('Treatment', 'MI_WT', 'Sham_WT'))
plotMA(oRes)
temp = as.data.frame(oRes)
table(rownames(temp) %in% dfResults$ind)
dfResults = dfResults[dfResults$ind %in% rownames(temp),]
temp = temp[dfResults$ind,]
identical((dfResults$ind), rownames(temp))
plot(dfResults$logFC, log(2^temp$log2FoldChange), pch=20)
table(oRes$padj < 0.01)
write.csv(oRes, file='results/DESeq.xls')

r1 = dfResults
r2 = oRes

r1 = r1[order(abs(r1$logFC), decreasing = T),]
head(r1)
r2$logFC = log(2^r2$log2FoldChange)
r2 = r2[order(abs(r2$logFC), decreasing = T),]

r1.top = r1$ind[1:100]
r2.top = rownames(r2)[1:100]
table(r1.top %in% r2.top)