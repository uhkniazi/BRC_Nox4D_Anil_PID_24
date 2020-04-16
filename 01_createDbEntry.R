# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 14/4/2020


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
# setwd('dataExternal/')
# setwd('fastq/')
# 
# # list the files
# cvFiles = list.files(pattern = 'fastq.gz', recursive = T)
# 
# # each sample has 2 files 
# fSplit = gsub('_R[1|2]', '', cvFiles)
# 
# lFiles = split(cvFiles, fSplit)
# 
# setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/Copy of Anil_metadata_Nox4D_total_OM_final.csv', header=T, stringsAsFactors = F)
str(dfMeta)

# remove white spaces 
cn = colnames(dfMeta)
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])

# sanity check
# table(as.character(dfMeta$Name) %in% cvFiles)

# ## order the table in the same sequence as file names
# i = match(cvFiles, dfMeta$Name)
# dfMeta = dfMeta[i,]
# identical(as.character(dfMeta$Name), cvFiles)
# dfMeta$fSplit = fSplit

## extract reduced sample table by removing one pair of the fastq file
i = grepl('_1.fastq.gz', x = dfMeta$File.name)
# cvSample = unique(dfMeta$fSplit)
# i = match(cvSample, dfMeta$fSplit)
dfMeta.sam = dfMeta[i,]
str(dfMeta.sam)
dfMeta.sam$File.name = gsub('_1.fastq.gz', replacement = '', dfMeta.sam$File.name)

xtabs( ~ File.name + Animal_unique_ID, data=dfMeta.sam)
xtabs( ~ Treatment + Genotype, data=dfMeta.sam)
xtabs( ~ Treatment + Parents_ID, data=dfMeta.sam)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta.sam$File.name, 
                       description= paste('sequencing_lane', as.character(dfMeta.sam$Lane),
                                          'group1 is Treatment',
                                          'group2 is Genotype',
                                          'group3 is RNA_isolation_batch',
                                          'Animal_ID', as.character(dfMeta.sam$Animal_unique_ID),
                                          'Parents_ID', as.character(dfMeta.sam$Parents_ID),
                                          'MI_procedure_day', as.character(dfMeta.sam$MI_procedure_day),
                                          sep=';'),
                       group1 = dfMeta.sam$Treatment, group2= dfMeta.sam$Genotype, 
                       group3=dfMeta.sam$RNA_isolation_batch)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 45;'))

# # create entries for these files in the database
# dbListTables(db)
# cn = dbListFields(db, 'File')[-1]
# cn
# identical(names(lFiles), dfMeta.sam$fSplit)
# names(lFiles) = dfSamples$id
# 
# # get the names of the samples
# temp = lapply(as.character(dfSamples$id), function(x){
#   # get the file names
#   df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[as.character(dfSamples$id) == x, 'id'])
#   return(df)
# })
# 
# dfFiles = do.call(rbind, temp)
# rownames(dfFiles) = NULL
# 
# # write this table to database
# ## note: do not execute as it is already done
# #dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
