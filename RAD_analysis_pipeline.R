# RAD_analysis_pipeline.R
#
# written by Steven Myers 
# 26.ix.2016
#
# Acknowledgements to Kyle Armstrong, Tessa Bradford, and Terry Bertozzi
#
# description: set missing data limits (filter based on missing data allowance), perform HWE and LD tests, PCA, generate summary stats, and export data for downstream analysis
#
# input: output from 'DArTtoR_v2.R'. 
# 
# R dependencies: adegenet, ggplot2, parallel



#############
### Steps ###

#1 load packages and data if necessary. 

#2. recode the loci as minor alleles and remove monomorphic and chimeric loci, and calculate minor allele frequencies

#3. remove individuals and loci based on user-defined level of missing data

#4. Identify putative clusters by PCA

#5. remove loci with missing data bias across clusters
# --> WRITE DATA TO FILE FOR USE IN NEIGHBOUR JOINING TREE PIPELINE

#6. calculate minor allele freq for reporting -- to accompany PCA and tree data

#7. Identify genetic clusters by PCA





setwd("C:/path/to/data/output/from/DArTtoR_v2.R")
options(max.print=99999999) 


###############
###  STEP 1 ###    <<< load libraries and data >>>

#install.packages("ape")
#install.packages("adegenet")
#install.packages("ggplot2")


library(ape)
library(adegenet)
library(ggplot2)
library(parallel)


# following commented section only necessary if you need to reload the data onto the gloabl environment (i.e. if it's not already loaded from running 'DArTtoR_v2.R')
#data = read.table('data_Filtered1row.csv', sep=',', header=F, row.names=NULL, stringsAsFactors=FALSE) # Load the data
#colnames(data) <- as.character(unlist(data[1,]))
#data <- t(data[2:nrow(data),1:ncol(data)])
#colnames(data) <- as.character(unlist(data[1,]))
#data <- data[2:nrow(data),1:ncol(data)]

x <- as.data.frame(t(data))                # quicker to process with less columns
x <- data.frame(sapply(x, as.integer))  

# DArT Report data is coded: "0" = Reference allele homozygote (major), "1"= SNP allele homozygote (minor), "2"= heterozygote and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# at this stage in the pipeline 0=1, 1=2, and 2=3
# However, this script was written for data coded: "0" = homozygote major, "1"= heterozygote, "2"= homozygote minor, NA = missing data; so we need to convert the data with this in mind

x[x=="1"] <- 0   #dataset will go from 1,2,3 to 0,2,3 
x[x=="3"] <- 1   #dataset will go from 0,2,3 to 0,2,1 

data <- t(data)

rownames(x) <- rownames(data)
colnames(x) <- colnames(data)

rm(data)



##############
### STEP 2 ###    <<< recode the loci as minor alleles and remove monomorphic and chimeric loci >>>

# calculate some important variables
a <- rowSums(x, na.rm=TRUE) # this calculates the sum of allele identities (0=major homo, 1=hetero, 2=min homo) for each locus 
b <- apply(x, 1, FUN=function(a)length(which(!is.na(a)))) # this calculates the number of individuals with successful calls per locus

# a) Remove monomorphic loci (DArT pipeline should've taken care of these; unless you've since split up your data).
y <- nrow(x)  # number of loci
x$a <-a   #adding the list of summed allele identities as a column in x
x$c <-2*b   #adding a column in x equal to 2 * number of individuals with successful calls per locus
x <-x[!(x$a=="0"),]   #removing all loci monomorphic for major allele
x <-x[!(x$a==x$c),]   #removing all loci monomorphic for minor allele (which would actually make it the major allele, but we haven't coded minor allele properly yet)

#how many loci were removed?
y-nrow(x)  

#tidy up (removing the added columns)
x <-subset(x, select=-a)
x <-subset(x, select=-c)

# b) calculate the minor allele frquency
minAlleleFreq <- a/(2*b)

# check how many if any loci are coded the wrong way (i.e. by major allele freq).  If returns 'NA', then you could have monomorphic loci or empty rows (although right now I'm not getting any NAs and I purposely put monomorphic loci in there!!!)
maf1 <-sum(minAlleleFreq > 0.5); maf1

# plot the minor allele frequencies
hist(minAlleleFreq, main=NULL,xlab="Minor allele frequency",ylab="Number of Loci",col=24,breaks=50,xaxt = "n")
axis(1, at = seq(0, 1, 0.1))  
# does it resemble a Rorschach test as expected? if not what are possible reasons for that?

# c) recode the loci as minor alleles
# first, save loci that will be recoded for changing back for conversion to IUPAC code
y <- x[minAlleleFreq > 0.5,]
y <- rownames(y)
# export list of loci for later reference
write.table(y, "data_Filtered1row_recodedLoci.csv", sep=",", na = "NA", row.names = F, col.names = F) 

for(n in 1:nrow(x)){     
  if (minAlleleFreq[n] > 0.5){
    recode <- (gsub(0,9,x[n,]))  #could use suppresswarnings() here because it will give warning 'NAs introduced by coercion' in every case a change was made.
    recode <- (gsub(2,0,recode))
    x[n,]  <- (as.integer(gsub(9, 2, recode)))
  }else{
    next #last loop doesn't seem to work properly unless this is here
  }
}
# this can take some time; up to ~1min per 10000 loci; you will get warnings that 'NAs introduced by coercion'

# d) re-check how many if any loci are coded the wrong way (hint: should be none!!).
a <- rowSums(x, na.rm=TRUE)  #sum of allele identities
b <- apply(x, 1, FUN=function(a)length(which(!is.na(a))))  # number of called individuals for given locus
minAlleleFreqFixed <- a/(2*b)
maf3 <-sum(minAlleleFreqFixed > 0.5); maf3 #maf3 should have no entries!


# e) remove chimeras
# search all loci with freq of 0.5 and remove loci with all individuals coded as a heterozygote (i.e. as 1s)
# this is same as last step in fixing min alleles but put here in case you skipped straight to this step
a <- rowSums(x, na.rm=TRUE)  #sum of allele identities
b <- apply(x, 1, FUN=function(a)length(which(!is.na(a))))  
minAlleleFreqFixed <- a/(2*b)

maf4 <-sum(minAlleleFreqFixed == 0.5); maf4 #calculates number of loci with allele frequencies of 0.5; if maf4>0 then you need to check for chimeras; if not you can skip to (f)

chimeraBool <- x=="1"   # converts table to boolean account of if allele call == 1
d <- rowSums(chimeraBool, na.rm=TRUE)    # sums the number of heterozygous allele calls per locus

x$d <-d   #adds the list of number of allele identities '1' as a column in x
x$b <-b   #adds the list of summed allele identities as a column in x
x <-x[!(x$b==x$d),]   #removing loci with every locus having all major hets or major homos

# tidy up (removing the added columns)
x <-subset(x, select=-b)
x <-subset(x, select=-d)

# how many chimeric loci were there?
a <- rowSums(x, na.rm=TRUE)  #sum of allele identities
b <- apply(x, 1, FUN=function(a)length(which(!is.na(a))))  # number of called individuals for given locus  (I think I found a simpler way to do this: b <- rowSums(!is.na(x))   !!!!)
minAlleleFreqFixed <- a/(2*b)
maf5 <-sum(minAlleleFreqFixed == 0.5)
maf4 - maf5

# plot the minor allele frequencies
hist(minAlleleFreqFixed, main=NULL,xlab="Minor allele frequency",ylab="Number of Loci",col=24,breaks=50,xaxt = "n")
axis(1, at = seq(0, 0.5, 0.1))   

# how many loci have a frequency >= 0.1?
sum(minAlleleFreqFixed >= 0.1)
# what percentage of the loci is that?
round(((sum(minAlleleFreqFixed >= 0.1)/nrow(x))*100),digits=2)

# tidy up the global environment
rm(a,b,maf1,maf3,maf4,maf5,minAlleleFreqFixed,minAlleleFreq,chimeraBool,d,n,recode)

# Export your data here to return to if needed
write.table(x, "data_Filtered1row_RC.csv", sep=",", na = "NA", row.names = TRUE, col.names = TRUE)

# if you need to reload the data, use the following:
#x <- read.table('data_Filtered1row_RC.csv', sep=',', header=T, row.names=1, stringsAsFactors=FALSE) # Load the data
#colnames(x) <- sub("X", "", colnames(x))




##############
### STEP 3 ###    <<< remove individuals and loci based on user-defined level of missing data >>>


# If you want, you can look at your data; how many loci have what proportion of missing data
loci <- nrow(x); loci
sample <- ncol(x); sample

locNumNAs <- cbind(apply(x, 1, function(y) sum(is.na(y)))) # calculate the number of loci with missing data for each sample (i.e. per sample amount of missing data)
x$locPropNAs <- locNumNAs / sample * 100 # calculate the proportion of missing data for each locus (i.e. per locus proportion of missing data) and add to dataframe as column
max <- max(x$locPropNAs, na.rm = FALSE); max   #use this value to set your breaks in the next step (it'll probably match the value you used for missing data in the step above)
# what's the max missing data?

hist(x$locPropNAs, c(seq(from = 0, to = ceiling(max), by = ceiling(max)/ceiling(max))), xlab="Proportion of missing data", ylab="number of loci")  # How many loci have what numNAs out of total loci?  You can use this to approximate how many loci you will get back if you choose a given percentage of missing data (or more importantly what proportion of loci you will be excluding by choosing that). 
sum(x$locPropNAs < 10)   # how many loci with < 10% missing data
sum(x$locPropNAs <= 10)   # how many with <=10% missing data

x <-subset(x, select=-locPropNAs) 
x = as.data.frame(t(x))

sampNumNAs <- cbind(apply(x, 1, function(y) sum(is.na(y)))) # calculate the number of loci with missing data for each sample (i.e. per sample amount of missing data)
x$sampPropNAs <- sampNumNAs / loci * 100 # calculate the proportion of missing data for each locus (i.e. per locus proportion of missing data) and add to dataframe as column
max <- max(x$sampPropNAs, na.rm = FALSE); max

hist(x$sampPropNAs, c(seq(from = 0, to = ceiling(max), by = ceiling(max)/ceiling(max))), xlab="Proportion of missing data", ylab="number of samples")  # How many loci have what numNAs out of total loci?  You can use this to approximate how many loci you will get back if you choose a given percentage of missing data (or more importantly what proportion of loci you will be excluding by choosing that). 
# examine the plot. investigate any samples that don't fit expected patterns or have lots of missing data.

#isolate samples that don't fit expected patterns or have lots of missing data (in this case it's samples with > 30% missing data).
y <-subset(x, !x$sampPropNAs < 30)
y <- as.data.frame(rownames(y))

# examine the samples in y.

# choose to leave or remove the samples (edit the script accordingly -- currently it's written as though the samples are being left in the data)

# prepare the data
x <-subset(x, select=-sampPropNAs) 
x <- as.data.frame(t(x))
y <- nrow(x)
locNumNAs <- cbind(apply(x, 1, function(y) sum(is.na(y)))) # calculate the number of loci with missing data for each sample (i.e. per sample amount of missing data)
x$locPropNAs <- locNumNAs / sample * 100 # calculate the proportion of missing data for each locus (i.e. per locus proportion of missing data) and add to dataframe as column

# set the allowed missing data
allowedMissing <- 15

#remove loci
x <-subset(x, x$locPropNAs <= allowedMissing)  #subsetting to include only LOCI with missing data <= the allowed missing data threshold

# how many loci were removed? (check the console for answer)
y-nrow(x)   
rm(y)
# how many loci are left? (check the console for answer)
loci <- nrow(x); loci     

# tidy up the data
x <-subset(x, select=-locPropNAs) 
x <- as.data.frame(t(x))   #transpose back

# how has removing those loci influenced the % missing data for our samples?

sampNumNAs <- cbind(apply(x, 1, function(y) sum(is.na(y)))) # calculate the number of loci with missing data for each sample (i.e. per sample amount of missing data)
x$sampPropNAs <- sampNumNAs / loci * 100 # calculate the proportion of missing data for each locus (i.e. per locus proportion of missing data) and add to dataframe as column
max <- max(x$sampPropNAs, na.rm = FALSE); max   

hist(x$sampPropNAs, breaks=c(seq(from = 0, to = ceiling(max), by = ceiling(max)/ceiling(max))), xlab="Proportion of missing data", ylab="number of samples")  # How many loci have what numNAs out of total loci?  You can use this to approximate how many loci you will get back if you choose a given percentage of missing data (or more importantly what proportion of loci you will be excluding by choosing that). 
# although the data still looks bimodal!!!

#tidy data
x <-subset(x, select=-sampPropNAs)
#transpose
x <- as.data.frame(t(x))   
# Export your data here to return to if needed
write.table(x, "data_Filtered1row_RC_L15MD.csv", sep=",", na = "NA", row.names = TRUE, col.names = TRUE)


