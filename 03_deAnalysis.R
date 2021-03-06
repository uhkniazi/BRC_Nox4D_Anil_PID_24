# File: 03_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 17/4/2020

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

# # setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef.1) + (1 | Coef.2), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + (1 | Coef.1), data=dfData)
# summary(fit.lme2)
# anova(fit.lme1, fit.lme2)
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
# l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 1, sigmaRan2=1)
# }

## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef.1), ]
d2 = dfData[!duplicated(dfData$Coef.2), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef.1),
                 Nclusters2=nlevels(dfData$Coef.2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches2 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef.1),
                 NgroupMap2=as.numeric(dfData$Coef.2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 NBatchMap2=as.numeric(d2$ind), # this is where we use the second level mapping
                 Nphi=nlevels(dfData$ind),
                 NphiMap=as.numeric(dfData$ind),
                 y=dfData$values, 
                 #gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

ptm = proc.time()

fit.stan = sampling(stanDso, data=lStanData, iter=1500, chains=2,
                    pars=c('sigmaRan1',
                           'sigmaRan2',
                           'phi',
                           #'mu',
                           'rGroupsJitter1',
                           'rGroupsJitter2',
                           'betas'
                           #'phi_scaled'
                    ),
                    cores=2, control=list(adapt_delta=0.99, max_treedepth = 11))#, init=initf)
save(fit.stan, file='results/fit.stan.nb_17Apr.rds')
ptm.end = proc.time()
print(fit.stan, c('sigmaRan1'), digits=3)
print(fit.stan, c('sigmaRan2'), digits=3)
print(fit.stan, c('phi'), digits=3)
print(fit.stan, c('rGroupsJitter1'))
traceplot(fit.stan, c('sigmaRan1[1]'))
traceplot(fit.stan, c('sigmaRan1[2]'))
traceplot(fit.stan, c('rGroupsJitter1[1]', 'sigmaRan1[1]'))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# # ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

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
l = tapply(d$cols, d$split, FUN = function(x, base='Sham_WT', deflection='MI_WT') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
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