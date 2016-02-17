
# **************************************************************************** #
# ***************             Start UP Commands                *************** #
# *************************************************************************** #

rm(list=ls())
now=Sys.Date(); today=format(now, format="%d%b%y")

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #

# NOTES: THIS ANALYSIS IS RELATED TO ALPHA DIVERSITY

# Set the working directory

# Working Directory
setwd("C:/Users/djlemas/Dropbox/03_Analysis/Microbiome/ADA_Infant_Study/")
list.files()

# work directory
#--------------
dir="C:/Users/djlemas/Dropbox/03_Analysis/Microbiome/ADA_Infant_Study/"

# Meta directory
#---------------
meta.dir=paste(dir,"Metadata/SEQ127",sep="");meta.dir

# OTU Directory
#--------------
otu.dir=paste(dir,"OTU/SEQ127",sep="");otu.dir

# **************************************************************************** #
# ***************                      LIBRARY                 *************** #
# **************************************************************************** #

#library(microbes)
library(plyr)
library(ape)
library(reshape2)
library(Matrix)
library(corrplot)

# **************************************************************************** #
# ***************                 GET META Data                *************** #
# **************************************************************************** #

# SEQ metadata 
#-------------
seq.meta=read.table(file=paste(meta.dir,"/SEQ127_metadata.txt",sep=""),header=T,sep="\t")
head(seq.meta);names(seq.meta)
seq.meta$otu.lib.name=as.character(seq.meta$Lib)
head(seq.meta);dim(seq.meta) #125*22

# Alpha Diversity Metadata
#-------------------------
div.meta=read.table(file=paste(meta.dir,"/SEQ127_mom_baby_all_time_all_taxa_AlphaDiversity.txt",sep=""),header=T,sep="\t")
head(div.meta);names(div.meta)
div.meta$otu.lib.name=as.character(div.meta$Collection)
head(div.meta);dim(div.meta) #121*22

# merge metadata
#---------------
new.meta=merge(seq.meta, div.meta, x.by="otu.lib.name", y.by="otu.lib.name", all.x=)
head(new.meta);names(new.meta)
meta=new.meta
names(meta)


# create mom/baby df
meta.baby=subset(meta, mom_baby=="baby" & !study_group=="NA")
meta.mom=subset(meta, mom_baby=="mom" & !study_group=="NA")

# mom by bmi group
#-----------------
meta.mom.n=meta.mom[,c(5,17,6,26,30,34,38,42,46,58,62,70,74,78,82)];names(meta.mom.n);head(meta.mom.n)
mom.bmi.diversity=daply(meta.mom.n, .(bmi_group), colwise(mean));mom.bmi.diversity

# mom by time
#------------
mom.time.diversity=daply(meta.mom.n, .(visit_number), colwise(mean));mom.time.diversity


# baby by bmi group
#------------------
meta.baby.n=meta.baby[,c(5,17,6,26,30,34,38,42,46,58,62,70,74,78,82)];names(meta.baby.n);head(meta.baby.n)
baby.bmi.diversity=daply(meta.baby.n, .(bmi_group), colwise(mean));baby.bmi.diversity

# baby by time
#-------------
baby.time.diversity=daply(meta.baby.n, .(visit_number), colwise(mean));baby.time.diversity



# reshape data
# format: ID/time/mom_baby/variable/value
#----------------------------------------

mdata <- melt(meta[,c(1,3,5:6,17,11,26:85)], id=c("study_id","visit_number","mom_baby","bmi_group"),na.rm=TRUE);head(mdata)
unique(mdata$variable)
mdata.bmi=subset(mdata, variable=="pre_preg_bmi" & mom_baby=="mom")
mdata.bmi$value=as.numeric(mdata.bmi$value)
mdata.bmi$bmi_group=as.factor(mdata.bmi$bmi_group)

# Compare mom BMI
#----------------
fit=lm(mdata.bmi$value~mdata.bmi$bmi_group);summary(fit)
boxplot(mdata.bmi$value,mdata.bmi$bmi_group)


# Beta Diversity
#---------------
beta.div=read.table(file=paste(meta.dir,"/SEQ127_mom_baby_all_time_all_taxa_BetaDiversity_BrayCurtis.txt",sep=""),header=T,sep="\t")
head(beta.div);names(beta.div);dim(beta.div)
row.names(beta.div)=names(beta.div)
head(beta.div)

melt(beta.div, varnames=names(dimnames(beta.div)))
dimnames(beta.div)

# Formatt
# comp1 comp2 variable value

a <- array(c(1:23, NA), c(2,3,4))
melt(a)
melt(a, na.rm = TRUE)
melt(a, varnames=c("X","Y","Z"))
dimnames(a) <- lapply(dim(a), function(x) LETTERS[1:x])
melt(a)
melt(a, varnames=c("X","Y","Z"))
dimnames(a)[1] <- list(NULL)
melt(a)

band{Matrix}

# parse triangular matrix

## A random sparse matrix :
set.seed(7)
m <- matrix(0, 5, 5)
m[sample(length(m), size = 14)] <- rep(1:9, length=14)
(mm <- as(m, "CsparseMatrix"))

tril(mm)        # lower triangle
tril(mm, -1)    # strict lower triangle
triu(mm,  1)    # strict upper triangle
band(mm, -1, 2) # general band
(m5 <- Matrix(rnorm(25), nc = 5))
tril(m5)        # lower triangle
tril(m5, -1)    # strict lower triangle

b=triu(m5, 1)     # strict upper triangle
b
band(m5, -1, 2) # general band
(m65 <- Matrix(rnorm(30), nc = 5))  # not square
triu(m65)       # result in not dtrMatrix unless square
(sm5 <- crossprod(m65)) # symmetric
band(sm5, -1, 1)# symmetric band preserves symmetry property
as(band(sm5, -1, 1), "sparseMatrix")# often preferable

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

cormat=beta.div
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

# Generate short form for data
#-----------------------------
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
test=flattenCorrMatrix(beta.div)
head(test)

beta.div2=as.matrix(beta.div)

dat<-data.frame(aa=rnorm(100),bb=rnorm(100),cc=rnorm(100),dd=rnorm(100))
dat$aa<-dat$aa+dat$dd
dat$cc<-dat$cc+dat$aa
cor.matrix(dat,test=cor.test)
cor.matrix(d(aa,cc),data=dat,test=cor.test,method="kendall")
cor.matrix(d(aa,cc),d(dd,bb),data=dat,test=cor.test,method="spearman")

head(beta.div)
test1=as.matrix(beta.div.m)
test2=corrplot(test1,method = "circle",type = "upper",order = "hclust")

git config --global user.email "dominicklemas@gmail.com"


