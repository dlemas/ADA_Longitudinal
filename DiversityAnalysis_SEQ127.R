
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
meta.all=meta
meta.baby=subset(meta, mom_baby=="baby" & !study_group=="NA")
meta.mom=subset(meta, mom_baby=="mom" & !study_group=="NA")

# mom by bmi group
#-----------------
meta.mom.n=meta.mom[,c(5,17,6,26,30,34,38,42,46,58,62,70,74,78,82)];names(meta.mom.n);head(meta.mom.n)
head(meta.mom.n)
mean(meta.mom.n$ShannonH_Mean)
sd(meta.mom.n$ShannonH_Mean)
length(meta.mom.n$ShannonH_Mean)
# Parameter estimates
matage <- ddply(meta.mom.n, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(meta.mom.n$ShannonH_Mean~as.factor(meta.mom.n$bmi_group))
wilcox.test(meta.mom.n$ShannonH_Mean~as.factor(meta.mom.n$bmi_group))

# 37wk
#------------
mom.37wk=subset(meta.mom.n, visit_number==1)
mean(mom.37wk$ShannonH_Mean)
sd(mom.37wk$ShannonH_Mean)
length(mom.37wk$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.37wk, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.37wk$ShannonH_Mean~as.factor(mom.37wk$bmi_group))
wilcox.test(mom.37wk$ShannonH_Mean~as.factor(mom.37wk$bmi_group))

# 2wk
#------------
mom.2wk=subset(meta.mom.n, visit_number==3)
mean(mom.2wk$ShannonH_Mean)
sd(mom.2wk$ShannonH_Mean)
length(mom.2wk$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.2wk, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.2wk$ShannonH_Mean~as.factor(mom.2wk$bmi_group))
wilcox.test(mom.2wk$ShannonH_Mean~as.factor(mom.2wk$bmi_group))

# 4mo
#------------
mom.4mo=subset(meta.mom.n, visit_number==5)
mean(mom.4mo$ShannonH_Mean)
sd(mom.4mo$ShannonH_Mean)
length(mom.4mo$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.4mo, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.4mo$ShannonH_Mean~as.factor(mom.4mo$bmi_group))
wilcox.test(mom.4mo$ShannonH_Mean~as.factor(mom.4mo$bmi_group))

# 12mo
#------------
mom.12mo=subset(meta.mom.n, visit_number==7)
mean(mom.12mo$ShannonH_Mean)
sd(mom.12mo$ShannonH_Mean)
length(mom.12mo$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.12mo, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.12mo$ShannonH_Mean~as.factor(mom.12mo$bmi_group))
wilcox.test(mom.12mo$ShannonH_Mean~as.factor(mom.12mo$bmi_group))


# baby by bmi group
#------------------
meta.baby.n=meta.baby[,c(5,17,6,26,30,34,38,42,46,58,62,70,74,78,82)];names(meta.baby.n);head(meta.baby.n)

mean(meta.baby.n$ShannonH_Mean)
sd(meta.baby.n$ShannonH_Mean)
length(meta.baby.n$ShannonH_Mean)
# Parameter estimates
matage <- ddply(meta.baby.n, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(meta.baby.n$ShannonH_Mean~as.factor(meta.baby.n$bmi_group))
wilcox.test(meta.baby.n$ShannonH_Mean~as.factor(meta.baby.n$bmi_group))


# 2wk
#------------
mom.2wk=subset(meta.baby.n, visit_number==3)
mean(mom.2wk$ShannonH_Mean)
sd(mom.2wk$ShannonH_Mean)
length(mom.2wk$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.2wk, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.2wk$ShannonH_Mean~as.factor(mom.2wk$bmi_group))
wilcox.test(mom.2wk$ShannonH_Mean~as.factor(mom.2wk$bmi_group))

# 4mo
#------------
mom.4mo=subset(meta.baby.n, visit_number==5)
mean(mom.4mo$ShannonH_Mean)
sd(mom.4mo$ShannonH_Mean)
length(mom.4mo$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.4mo, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.4mo$ShannonH_Mean~as.factor(mom.4mo$bmi_group))
wilcox.test(mom.4mo$ShannonH_Mean~as.factor(mom.4mo$bmi_group))

# 12mo
#------------
mom.12mo=subset(meta.baby.n, visit_number==7)
mean(mom.12mo$ShannonH_Mean)
sd(mom.12mo$ShannonH_Mean)
length(mom.12mo$ShannonH_Mean)
# Parameter estimates
matage <- ddply(mom.12mo, c("bmi_group"), summarise,
                N    = length(ShannonH_Mean),
                mean = round(mean(ShannonH_Mean),digits=4),
                sd   = round(sd(ShannonH_Mean),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(mom.12mo$ShannonH_Mean~as.factor(mom.12mo$bmi_group))
wilcox.test(mom.12mo$ShannonH_Mean~as.factor(mom.12mo$bmi_group))


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

beta.div.est=read.table(file=paste(meta.dir,"/alpha_div_estimate_18Feb16.csv",sep=""),header=T,sep=",")
head(beta.div.est);names(beta.div.est);dim(beta.div.est)

library(ggplot2)

##   supp dose   len       sd
## 1   OJ  0.5 13.23 4.459709
## 2   OJ    1 22.70 3.910953
## 3   OJ    2 26.06 2.655058
## 4   VC  0.5  7.98 2.746634
## 5   VC    1 16.77 2.515309
## 6   VC    2 26.14 4.797731

test.dat=beta.div.est
# mom-baby
test2=test.dat
test2$Factor1 <- factor(test2$Factor1, c("mom", "baby"))
test2$Factor2 <- factor(test2$Factor2, c("37wk","2wk", "4mo","12mo","all"))

library(ggplot2)

# Standard deviation of the mean as error bar
p <- ggplot(test2, aes(x=Factor2, y=Mean, fill=Factor1)) 
m=p+geom_bar(stat="identity", position=position_dodge())
m + labs(title="Plot of ShannonH by Visit", x="Study Visit", y = "ShannonH")+ scale_fill_brewer(palette="Reds") + theme_minimal()

  
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(.9));m
n=m+geom_text(aes(label=N), position=position_dodge(width=0.9), vjust=16);n
n + labs(title="Plot of BC Dissimilarity  by Visit", x="Study Visit", y = "Bray-Cutis Dissimilarity")+ scale_fill_brewer(palette="Reds") + theme_minimal()




# **************************************************************************** #
# ***************                 BETA DIVERSITY               *************** #
# **************************************************************************** #
beta.div=read.table(file=paste(meta.dir,"/SEQ127_mom_baby_all_time_all_taxa_BetaDiversity_BrayCurtis.txt",sep=""),header=T,sep="\t")
head(beta.div);names(beta.div);dim(beta.div)
row.names(beta.div)=names(beta.div)
head(beta.div)

melt(beta.div, varnames=names(dimnames(beta.div)))
dimnames(beta.div)


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
test=as.data.frame(flattenCorrMatrix(beta.div))
head(test)
test$row=as.character(test$row)

# Format Beta.Div data
test$column=as.character(test$column)
drop.list=c("ADA51B2mo","ADA51B2wk","ADA51B4mo","ADA51B8mo","ADA51M4mo","NEG1","NEG2","NEG3")
test=test[-which(test$row %in% drop.list), ]


# Format Meta Data
mdata.lib=subset(mdata, variable=="otu.lib.name");head(mdata.lib)
mdata.lib=mdata.lib[-which(is.na(mdata.lib$bmi_group)), ]
mdata.lib$mom_baby=as.character(mdata.lib$mom_baby)
mdata.lib$study_id=as.character(mdata.lib$study_id)

# Index Meta Data
myIndex=length(mdata.lib$value);myIndex

# mom baby
test$row_mom_baby=NA
test$col_mom_baby=NA

# bmi group
test$row_bmi_group=NA
test$col_bmi_group=NA

# visit number
test$row_visit_number=NA
test$col_visit_number=NA

# participant ID
test$row_id=NA
test$col_id=NA
  

for (i in 1:myIndex){
 # matching
 row.match.logic=grepl(mdata.lib$value[i],test$row )
 col.match.logic=grepl(mdata.lib$value[i],test$column )
 
 # mom baby
 test$row_mom_baby=ifelse(row.match.logic==TRUE,mdata.lib$mom_baby[i],test$row_mom_baby)
 test$col_mom_baby=ifelse(col.match.logic==TRUE,mdata.lib$mom_baby[i],test$col_mom_baby)
 
 # bmi_grouvp
 test$row_bmi_group=ifelse(row.match.logic==TRUE,mdata.lib$bmi_group[i],test$row_bmi_group)
 test$col_bmi_group=ifelse(col.match.logic==TRUE,mdata.lib$bmi_group[i],test$col_bmi_group)

 # study id
 test$row_id=ifelse(row.match.logic==TRUE,mdata.lib$study_id[i],test$row_id)
 test$col_id=ifelse(col.match.logic==TRUE,mdata.lib$study_id[i],test$col_id)

 # visit number
 test$row_visit_number=ifelse(row.match.logic==TRUE,mdata.lib$visit_number[i],test$row_visit_number)
 test$col_visit_number=ifelse(col.match.logic==TRUE,mdata.lib$visit_number[i],test$col_visit_number)

}

# Reorganize data.frame
data.final=test[,c(7,1,4:6,10,2,8,9,11,3)]
# check!
head(data.final)
dim(data.final)
str(data.final)
unique(data.final$row_id)

data1=data.final[-which(is.na(data.final$row)), ]


# FINALLY COMPLETED
#------------------

# lets compute mean correaltions by groups!
#------------------------------------------

#####################
## Mom-baby Pairs 
###
#####################

# Subset to mom-baby pairs for all time-points
data2=data.final[which(data.final$row_id==data.final$col_id & data.final$row_visit_number==data.final$col_visit_number),]
dim(data2);head(data2);length(data2[,1]) # 29
all.mean=mean(data2$cor);all.mean  # 0.8701
all.sd=sd(data2$cor);all.sd        # 0.0866
data2

# Parameter estimates
matage <- ddply(data2, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8727/sd=0.0841
# Ob=0.8633/sd=0.0.0348

# significance test
t.test(data2$cor~as.factor(data2$col_bmi_group))
wilcox.test(data2$cor~as.factor(data2$col_bmi_group))
# t.test/p.value=0.8146
# wilcoxo/p.value=0.7927

# How many unique time points
unique(data2$row_visit_number) # 5, 3, 7

# 2wks Correlation
data.2wk=subset(data2, row_visit_number==3);length(data.2wk[,1])  # 10 people
all.mean.2wk=mean(data.2wk$cor);all.mean.2wk # 0.8937
all.sd.2wk=sd(data.2wk$cor);all.sd.2wk       # 0.0552

# Parameter estimates
matage <- ddply(data.2wk, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.9128/sd=0.0436
# Ob=0.8492/sd=0.0616

# significance test
t.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
wilcox.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
# t.test/p.value=0.206
# wilcoxo/p.value=0.117


# 4-month Correlation
data.4mo=subset(data2, row_visit_number==5);length(data.4mo[,1]) # 13 people
all.4mo.mean=mean(data.4mo$cor);all.4mo.mean         # 0.8783
all.4mo.sd=sd(data.4mo$cor);all.4mo.sd               # 0.0983

# Parameter estimates
matage <- ddply(data.4mo, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.8898/sd=0.0289
# Ob=0.8524/sd=0.0658

# significance test
t.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))
wilcox.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))

# t.test/p.value=0.629
# wilcoxo/p.value=0.604

# 12-months Correaltions
data.12mo=subset(data2, row_visit_number==7);length(data.12mo[,1]) # 6 people
all.12mo.mean=mean(data.12mo$cor);all.12mo.mean         # 0.8130
all.12mo.sd=sd(data.12mo$cor);all.12mo.sd               # 0.0900

# Parameter estimates
matage <- ddply(data.12mo, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.7858/sd=0.0678
# Ob=0.9488/sd=0.0

# significance test
t.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))
wilcox.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))
# t.test/p.value=NA
# wilcoxo/p.value=0.333
  


#####################
## mom Pairs 
###
#####################

# Subset to mom-baby pairs for all time-points
data3=data.final[which(data.final$row_mom_baby=="mom" & data.final$col_mom_baby=="mom"),]
dim(data3);head(data3)
unique(data3$col_visit_number)   # 5 3 1 7

length(data3$cor)    # 1128
all.mom.mean=mean(data3$cor);all.mom.mean         # 0.5222
all.mom.sd=sd(data3$cor);all.mom.sd               # 0.1357

# Parameter estimates
matage <- ddply(data3, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.5110/sd=0.1358
# Ob=0.5570/sd=0.1296

# significance test
t.test(data3$cor~as.factor(data3$col_bmi_group))
wilcox.test(data3$cor~as.factor(data3$col_bmi_group))

# t.test/p.value=<0.001
# wilcoxo/p.value=0.001

# 37 weeks
data.37wk=subset(data3, row_visit_number==1 & col_visit_number==1);length(data.37wk[,1]) # 91 people
all.mom.mean=mean(data.37wk$cor);all.mom.mean         # 0.5559
all.mom.sd=sd(data.37wk$cor);all.mom.sd               # 0.1281

# Parameter estimates
matage <- ddply(data.37wk, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.5361/sd=0.1260 66
# Ob=0.6081/sd=0.1210 25

# significance test
t.test(data.37wk$cor~as.factor(data.37wk$col_bmi_group))
wilcox.test(data.37wk$cor~as.factor(data.37wk$col_bmi_group))
# t.test/p.value=0.0159
# wilcoxo/p.value=0.016

# 2-week
data.2wk=subset(data3, row_visit_number==3 & col_visit_number==3);length(data.2wk[,1]) # 91 people
all.mom.mean=mean(data.2wk$cor);all.mom.mean         # 0.5605
all.mom.sd=sd(data.2wk$cor);all.mom.sd               # 0.1115

# Parameter estimates
matage <- ddply(data.2wk, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.5628/sd=0.1198 70
# Ob=0.5529/sd=0.0795 21


# significance test
t.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
wilcox.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
# t.test/p.value=0.660
# wilcoxo/p.value=0.966

# 4-month
data.4mo=subset(data3, row_visit_number==5 & col_visit_number==5);length(data.4mo[,1]) # 78 people
all.mom.mean=mean(data.4mo$cor);all.mom.mean         # 0.4394
all.mom.sd=sd(data.4mo$cor);all.mom.sd               # 0.1221

# Parameter estimates
matage <- ddply(data.4mo, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.4036/sd=0.0936 57
# Ob=0.5365/sd=0.1357 21

# Significance Test
t.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))
wilcox.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))
# t.test/p.value=<0.001
# wilcoxo/p.value=<0.001

# 12-month
data.12mo=subset(data3, row_visit_number==7 & col_visit_number==7);length(data.12mo[,1]) # 21 people
all.mom.mean=mean(data.12mo$cor);all.mom.mean         # 0.4287
all.mom.sd=sd(data.12mo$cor);all.mom.sd               # 0.0915

# Parameter estimates
matage <- ddply(data.12mo, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.4199/sd=0.0951 18
# Ob=0.4823/sd=0.0422 3

# significance test
t.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))
wilcox.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))
# t.test/p.value=0.106
# wilcoxo/p.value=0.153

#####################
## infant Pairs 
###
#####################

# Subset to mom-baby pairs for all time-points
data4=data.final[which(data.final$row_mom_baby=="baby" & data.final$col_mom_baby=="baby"),]
dim(data4);head(data4) 
length(data4$cor)# 2278


all.mom.mean=mean(data4$cor);all.mom.mean         # 0.7279
all.mom.sd=sd(data4$cor);all.mom.sd         # 0.1525

# Parameter estimates
matage <- ddply(data4, c("col_bmi_group"), summarise,
                N    = length(cor),
                mean = round(mean(cor),digits=4),
                sd   = round(sd(cor),digits=4),
                se   = round(sd / sqrt(N),digits=4) );matage
# NW=0.7272/sd=0.1525 1746
# Ob=0.7304/sd=0.1527 532

# signficnce test
t.test(data4$cor~as.factor(data4$col_bmi_group))
wilcox.test(data4$cor~as.factor(data4$col_bmi_group))
# NW=0.727
# ob=0.730
# t.test/p.value=0.672
# wilcoxo/p.value=0.667

# 2-week
data.2wk=subset(data4, row_visit_number==3 & col_visit_number==3);length(data.2wk[,1]) # 55 people
all.mom.mean=mean(data.2wk$cor);all.mom.mean         # 0.7465
all.mom.sd=sd(data.2wk$cor);all.mom.sd               # 0.1573

# significance test
t.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
wilcox.test(data.2wk$cor~as.factor(data.2wk$col_bmi_group))
# t.test/p.value=0.397
# wilcoxo/p.value=0.958

# 4-month
data.4mo=subset(data4, row_visit_number==5 & col_visit_number==5);length(data.4mo[,1]) # 120 people
all.mom.mean=mean(data.4mo$cor);all.mom.mean         # 0.7235
all.mom.sd=sd(data.4mo$cor);all.mom.sd               # 0.1562

# significance test
t.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))
wilcox.test(data.4mo$cor~as.factor(data.4mo$col_bmi_group))
# t.test/p.value=0.902
# wilcoxo/p.value=0.953

# 12-month
data.12mo=subset(data4, row_visit_number==7 & col_visit_number==7);length(data.12mo[,1]) # 45 people
all.mom.mean=mean(data.12mo$cor);all.mom.mean         # 0.6663
all.mom.sd=sd(data.12mo$cor);all.mom.sd         # 0.1551

# Significance test
t.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))
wilcox.test(data.12mo$cor~as.factor(data.12mo$col_bmi_group))

# t.test/p.value=0.821
# wilcoxo/p.value=0.987


# **************************************************************************** #
# ***************             IMPORT BETA DIVERSITY ESTI       *************** #
# **************************************************************************** #
beta.div.est=read.table(file=paste(meta.dir,"/beta_div_estimate_18Feb16.csv",sep=""),header=T,sep=",")
head(beta.div.est);names(beta.div.est);dim(beta.div.est)

library(ggplot2)

##   supp dose   len       sd
## 1   OJ  0.5 13.23 4.459709
## 2   OJ    1 22.70 3.910953
## 3   OJ    2 26.06 2.655058
## 4   VC  0.5  7.98 2.746634
## 5   VC    1 16.77 2.515309
## 6   VC    2 26.14 4.797731

test.dat=beta.div.est[,c(1,2,3:6)]
test.dat$Visit=as.factor(test.dat$Visit)
# mom-baby
test2=test.dat
test2$Factor1 <- factor(test2$Factor1, c("mom", "baby", "mom-baby"))
test2$Factor2 <- factor(test2$Factor2, c("37wk","2wk", "4mo","12mo","all"))


# Standard deviation of the mean as error bar
p <- ggplot(test2, aes(x=Factor2, y=Mean, fill=Factor1)) 
m=p+geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(.9));m
n=m+geom_text(aes(label=N), position=position_dodge(width=0.9), vjust=16);n
n + labs(title="Plot of BC Dissimilarity  by Visit", x="Study Visit", y = "Bray-Cutis Dissimilarity")+ scale_fill_brewer(palette="Reds") + theme_minimal()


# http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization




