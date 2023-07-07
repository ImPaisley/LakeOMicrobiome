###### BATCH CORRECTION & ASSOCIATED ANALYSES ######

##First had to go through and manually assign batches to the samples within the
##metadata file (based on mapping files)
## 10 KNOWN SEQUENCING RUNS IN TOTAL (an unknown sequence run making 11 "UNK")

###### SET WORKING DIRECTORY AND SEED ####
setwd("F:/Paise_Thesis/LakeO_Data/2019-2021_LakeO_Data/Analyses/LakeO_BatchCorrected/Analyses_Corrected")
#or setwd("/Volumes/PaiseSSD-T7/Paise_Thesis/LakeO_Data/2019-2021_LakeO_Data/Analyses/LakeO_BatchCorrected/Analyses_Corrected") for use on the lab computer
set.seed(1998)

###### Packages ######
library(vegan)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(BiocManager)
library(MMUPHin)

#updating BiocManager and installing mmuphin
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("MMUPHin")

###### Creating relative abundance data ######
set.seed(1998)
dat<-read.csv("feature_Y123_nobcmASVs-nobelow10korDupes.csv", header=TRUE, row.names = 1)
dat<-data.matrix(dat)
typeof(dat) #"integer"
dat <- t(dat)
row.names(dat) # row names should now be the sample names
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
typeof(metadata) ## "list"
dat <- as.data.frame(dat)
typeof(dat)
common.rownames <- intersect(rownames(dat), rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata))
otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund] #dominant taxa
dat.pa<-decostand(dat.dom, method ="pa") #presence/absence data
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
dat.01per<-dat.dom[,dat.otus.01per] #removed ASVs that occur less than 0.1%; 8,340 taxa present
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat.dom[,dat.otus.001per] #removed ASVs that occur less than 0.01%; 44,623 taxa present
                                      #increases the number of ASVs - includes more "microdiversity" 
dat.ra<-decostand(dat.01per, method = "total") #relative abundance of >1% taxa

###### ANOSIM by Sequencing Batch ######
set.seed(1998)
##create relative abundance table in above code
##create Bray-Curtis dissimilarity distance matrix
ra.bc.dist<-vegdist(dat.ra, method = "bray")

##betadisper calculates dispersion (variances) within each group 
dis.Batch <- betadisper(ra.bc.dist,metadata$Batch)

##permutest determines if the variances differ by groups (If differences are SIGNIFICANT - use ANOSIM
##                                                        if not use PERMANOVA (adonis))
permutest(dis.Batch, pairwise=TRUE, permutations=999)
#           Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     10 1.0196 0.101957 10.832    999  0.001 ***   SIGNIFICANT - USE ANOSIM!!
# Residuals 530 4.9886 0.009413                         
# ---

# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# ELIZA2    ELIZA23     ELIZA3       LO22      LO310     LO8382       NOAA      PAIS1      PAIS2      PAIS3   UNK
# ELIZA2             2.5800e-01 2.2800e-01 4.0000e-03 1.0000e-03 2.8000e-02 7.0300e-01 5.4400e-01 2.0000e-02 2.0000e-03 0.001
# ELIZA23 2.4836e-01            2.9000e-02 4.5000e-02 1.0000e-03 9.7000e-02 3.1600e-01 4.7000e-01 1.6100e-01 4.9000e-02 0.001
# ELIZA3  1.9314e-01 2.3483e-02            3.0000e-03 1.0000e-03 6.0000e-03 4.5400e-01 1.1300e-01 3.0000e-03 1.0000e-03 0.001
# LO22    4.2982e-03 4.4251e-02 1.9370e-03            1.1800e-01 9.9300e-01 7.5000e-02 1.0000e-02 5.1900e-01 7.1500e-01 0.005
# LO310   1.4327e-04 5.5957e-04 8.7490e-06 1.2846e-01            1.8300e-01 3.0000e-03 1.0000e-03 3.1000e-02 4.1000e-02 0.966
# LO8382  2.2967e-02 9.5642e-02 4.8732e-03 9.9098e-01 1.8031e-01            1.1600e-01 3.4000e-02 6.0500e-01 7.5900e-01 0.029
# NOAA    7.1430e-01 3.1669e-01 4.7882e-01 7.6250e-02 3.0428e-03 1.0361e-01            5.3600e-01 1.1600e-01 7.0000e-02 0.001
# PAIS1   5.5294e-01 4.9884e-01 9.5982e-02 8.2362e-03 2.4554e-04 3.9403e-02 5.4000e-01            5.3000e-02 4.0000e-03 0.001
# PAIS2   1.9845e-02 1.6321e-01 2.9914e-03 5.0465e-01 2.4943e-02 5.7184e-01 1.0857e-01 4.4072e-02            7.0500e-01 0.002
# PAIS3   1.6726e-03 4.7249e-02 1.0649e-03 6.9056e-01 3.8967e-02 7.4158e-01 7.1715e-02 3.8772e-03 7.0437e-01            0.001
# UNK     6.1433e-14 3.1874e-10 7.8632e-11 4.5319e-03 9.6747e-01 2.3258e-02 3.6162e-05 3.0384e-14 5.3552e-05 3.8664e-05      

##ANOSIM - determining if the differences between two or more groups are significant. 
## The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to
## the mean of ranked dissimilarities within groups. An R value close to “1" suggests 
## dissimilarity between groups while an R value close to “0” suggests an even distribution of
## high and low ranks within and between groups”
## the higher the R value, the more dissimilar your groups are in terms of microbial community composition.

anosim(ra.bc.dist, metadata$Batch, permutations = 999)
# ANOSIM statistic R: 0.1486
# Significance: 0.001

anosim(ra.bc.dist, metadata$Batch, permutations = 9999)
# ANOSIM statistic R: 0.1486
# Significance: 0.0001

## Conclusion? There are significantly weak differences between batches so the
## data needs to be batch corrected and ALL analyses redone.

###### BATCH CORRECTION ######
set.seed(1998)
library(MMUPHin)
library(vegan)

## Loading in feature- and metadata
dat <- read.csv("feature_Y123_nobcmASVs-nobelow10korDupes.csv", header=TRUE, row.names = 1)
dat <- data.matrix(dat)
typeof(dat) #"integer"
dat <- t(dat) #transposing data matrix
row.names(dat) # row names should now be the sample names
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
typeof(metadata) ## "list"
dat <- as.data.frame(dat)
typeof(dat)
common.rownames <- intersect(rownames(dat), rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata)) #TRUE

## Batch Correction (following Harvard tutorial)
#looking at how many samples are in each batch
table(metadata$Batch)
# ELIZA2 ELIZA23  ELIZA3   LO22   LO310   LO8382    NOAA   PAIS1   PAIS2   PAIS3   UNK 
#   62      50      11      38      20      20       6      98      40      72     124 

#Adjusting (removing) batch effect
#taxa should be rows in feature table and samples should be rows in metadata
#feature table should be a matrix while metadata should be a dataframe
fit_adjust_batch <- adjust_batch(feature_abd = t(dat), 
                                 batch = "Batch", 
                                 data = metadata)
Lake_abd_adj <- fit_adjust_batch$feature_abd_adj #now adjusted feature table MATRIX
Lake_abd_adj <- as.data.frame(Lake_abd_adj) #converting to data frame
write.csv(Lake_abd_adj, "feature_Y123_ADJUSTED.csv") #saving as csv



###### Creating a rarefaction curve on the read counts ######
library(vegan)

#load in data with NO blank samples or blank ASVs
rardat<-read.csv("feature_Y123_noblanksorbASVs.csv", header=TRUE, row.names=1, sep=',') 

#as you can see the samples are in columns and need to be in the rows so we need to flip or transpose the file
#transpose the data to rows 
trans.rardat <- t(rardat)
## check file to make sure it worked 
trans.rardat[1:5,1:5] #shows rows 1 through 5 and the samples should now be the rows
##making the transformed data matrix into main
rardat <- trans.rardat
##changing back into data frame instead of matrix (transforming the data frame turned it into a matrix)
rardat <-as.data.frame(rardat)
#check data file to make sure it looks okay 
View(rardat)

rowSums(rardat) #sums the value of each row in the data frame

#### Creating the rarefaction curve
#count the number of species within each sample
S <- specnumber(rardat)
raremax <- min(rowSums(rardat)) ## takes the sample with the lowest sample size which is 0 in this dataset

#creating color palette for curve
colors() ## lists the color names that are built into R
cc <- palette()
palette(c(cc,"purple","brown"))    ## creating the color ramp for the plot
cc <- palette()

#plotting the rarefaction curves
## auto removes samples that have no reads
pars <- expand.grid()
Hklim <- rarecurve(rardat, step = 2000, sample=raremax, col = cc, label = TRUE, main="Rarefaction Curve for Lake O read counts", 
          cex= 0.14, cex.axis= 0.7, cex.lab= 1, xlim=c(0,100000), xlab = "# of Reads", ylab = "# of ASVs", tidy = T)



#### ####



 





###### ANALYSES ON BATCH CORRECTED DATA ######
###### SET WORKING DIRECTORY AND SEED ####
setwd("F:/Paise_Thesis/LakeO_Data/2019-2021_LakeO_Data/Analyses/LakeO_BatchCorrected/Analyses_Corrected")
#or setwd("/Volumes/PaiseSSD-T7/Paise_Thesis/LakeO_Data/2019-2021_LakeO_Data/Analyses/LakeO_BatchCorrected/Analyses_Corrected")
#for use on the lab computer
set.seed(1998)

###### PACKAGES ######
library(phyloseq)
library(vegan)
library(ggplot2)
library(tidyverse)
library(RVAideMemoire)
library(DESeq2)
library(corrplot)
library(multcompView)
library(pgirmess)
library(data.table)
library(microbiome)
library(BiocManager)
library(ggthemes)
library(gplots)
library(RColorBrewer)
library(cooccur)
library(visNetwork)
library(Hmisc)
library(cowplot)
library(reshape2)
library(sjmisc)
library(MASS)
library(scales)
library(forcats)
library(leaflet)
library(eulerr)
library(microbiomeutilities)

##Installing packages
BiocManager::install("DESeq2")
BiocManager::install("lefser")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
BiocManager::install("microbiomeutilities")

##Had to install using binaries (3/9/23 on iMAC)
install.packages("tibble", type="binary")
install.packages("Hmisc", type="binary")

## Notes on packages:
# pgirmess = Kruskal-Wallis Test
# RVAideMemoire = PERMANOVA
# cowplot = making multiple plots using ggplots objects

###### GENERATING CITATIONS FOR R AND R PACKAGES ######
citation() #retrieves the citation for R 
## Citation for each package
citation("ggplot2")
citation("phyloseq")
citation("vegan")
citation("microbiome")
citation("MMUPHin")
citation("pgirmess")
citation("multcompView")
citation("RVAideMemoire")

###### Prepping data for analyses #####

### import feature-table data ### 
##change to csv or import as a tsv using read.table function
dat<-read.csv("feature_Y123_ADJUSTED.csv", header=TRUE, row.names = 1) ## do not add "header =" or "row.names =" for merging 
# 561 samples; 65294 taxa 

dat<-data.matrix(dat) ##if data is not recognized as a data.frame numeric 
typeof(dat) #"integer"
#check data file to make sure it looks okay 


#as you can see the samples are in columns and need to be in the rows so we need to flip or transpose the file
#transpose the data to rows 
trans.dat <- t(dat)

## check file to make sure it worked 
trans.dat[1:5,1:5] #shows rows 1 through 5 and the samples should now be the rows


##set transposed data to main data variable 
dat <-trans.dat
row.names(dat) # row names should now be the sample names

### import metadata ###
###(if you intend to do any statistical analyses in R)
##If not skip to refining and normalizing steps 
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)

##should read "list"
typeof(metadata) ## "list"
dat <- as.data.frame(dat) ## had to change dat back into a data frame to check for matching rows
typeof(dat) ## "list"

##check to make sure the sample names match and are correct
common.rownames <- intersect(rownames(dat), rownames(metadata))
##541 rows are in common (20 S80 samples NOT included)

##if there are any rows that do not match, they will not be included in the statistical analysis or relative abundance tables
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]

##check that all rows match
all.equal(rownames(dat),rownames(metadata)) #TRUE so yes they all match
dat[1:5,1:3] ## double-checking that everything looks good



##merging the working feature and taxonomy tables
feat <- dat
tax <- read.csv("taxonomy_Y123_edited&cleaned.csv")
feattax <- merge.data.frame(feat, tax, by= "FeatureID", all.x=TRUE, all.y = TRUE)
write.csv(feattax, "feat-tax_Y123_cleaned.csv")


## CONTINUE HERE IF YOU ARE IGNORING METADATA ###
## refining and normalizing data #
##remove singletons and doubletons -  ASVs that only show up once or twice 
##this can be modified or removed if desired. Depends on what you want to know 
library(vegan)

otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund] #46838 taxa

##all this will get rid of ASVs that appear less than a certain percent in the data 
##this is not always something that you should do depending on your question. 
dat.pa<-decostand(dat.dom, method ="pa")  #"pa" = standardization method that scales your data to presence/absence (0/1)
##remove ASVs that occur <0.01 ***
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
dat.01per<-dat.dom[,dat.otus.01per]
# 8,340 taxa
write.csv(as.data.frame(t(dat.01per)), "feature_Y123_0.01per.csv")

##remove ASVs that occur <0.001 ---> increases the number of ASVs - includes more "micro-diversity" 
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat.dom[,dat.otus.001per]
# 46,838 taxa

## relative abundance --> normalization ##
dat.ra<-decostand(dat.01per, method = "total") #"total" = standardization method that divides your data by margin total (def. margin = 1)

##export relative abundance table(s)
write.csv(dat.ra, "relative-abundance_Y123.csv")

## SHORTCUT WITH NO EXPLANATIONS
## re-creating relative abundance table
set.seed(1998)
dat<-read.csv("feature_Y123_ADJUSTED.csv", header=TRUE, row.names = 1)
dat<-data.matrix(dat)
typeof(dat) 
dat <- t(dat)
row.names(dat) 
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
typeof(metadata) 
dat <- as.data.frame(dat)
typeof(dat)
common.rownames <- intersect(rownames(dat), rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata))
otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund] 
dat.pa<-decostand(dat.dom, method ="pa")
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
dat.01per<-dat.dom[,dat.otus.01per]
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat.dom[,dat.otus.001per]
dat.ra<-decostand(dat.01per, method = "total") 


###### Merging relative abundance with taxonomy and getting averages ######
Yr1 <- read.csv("Year1_RA.csv")
Yr2 <- read.csv("Year2_RA.csv")
Yr3 <- read.csv("Year3_RA.csv")
tax <- read.csv("taxonomy_Y123_edited&cleaned.csv")
Yr1t <- merge.data.frame(Yr1,tax,by= "FeatureID", all.x = TRUE)
Yr2t <- merge.data.frame(Yr2,tax,by= "FeatureID", all.x = TRUE)
Yr3t <- merge.data.frame(Yr3,tax,by= "FeatureID", all.x = TRUE)
write.csv(Yr1t, "Year1_RA.csv")
write.csv(Yr2t, "Year2_RA.csv")
write.csv(Yr3t, "Year3_RA.csv")

### Average and St.dev abundance of each phylum in each year
library(tidyverse)

## Year 1
#first merge data with matching taxonomy and load csv
Yr1 <- read.csv("Year1_RA.csv", row.names = 1)
#Sum by phylum across samples
physumY1 <- Yr1 %>% 
  group_by(Phylum) %>% 
  summarise(across(where(is.numeric), sum))
#Average phylum across samples
Y1mean <- apply(physumY1[,-1], 1, mean, na.rm=TRUE)
#Standard deviation across samples 
Y1std <- apply(physumY1[,-1], 1, sd, na.rm=TRUE)
#merge average and st.dev with rows
Y1avsd <- as.data.frame(cbind(physumY1$Phylum,Y1mean, Y1std))
#Renaming columns and saving as csv
colnames(Y1avsd)[1] ="Phylum"
colnames(Y1avsd)[2] ="Average"
colnames(Y1avsd)[3] ="Stand.Dev"
write.csv(Y1avsd, "Year1_AvSD-UPDATED.csv")
#Extract top 10 phyla and save as csv
top101 <- names(top10phy.names.Y1)
Y1avsd10 <- filter(Y1avsd,
                   Y1avsd$Phylum %in% top101)  
write.csv(Y1avsd10, "Year1_AvSD_TOP10-UPDATED.csv")

## Year 2
Yr2 <- read.csv("Year2_RA.csv", row.names = 1)
#Sum by phylum across samples
physumY2 <- Yr2 %>% 
  group_by(Phylum) %>% 
  summarise(across(where(is.numeric), sum))
#Average phylum across samples
Y2mean <- apply(physumY2[,-1], 1, mean, na.rm=TRUE)
#Standard deviation across samples 
Y2std <- apply(physumY2[,-1], 1, sd, na.rm=TRUE)
#merge average and st.dev with rows
Y2avsd <- as.data.frame(cbind(physumY2$Phylum,Y2mean, Y2std))
#Renaming columns and saving as csv
colnames(Y2avsd)[1] ="Phylum"
colnames(Y2avsd)[2] ="Average"
colnames(Y2avsd)[3] ="Stand.Dev"
write.csv(Y2avsd, "Year2_AvSD-UPDATED.csv")
#Extract top 10 phyla and save as csv
top102 <- names(top10phy.names.Y2)
Y2avsd10 <- filter(Y2avsd,
                   Y2avsd$Phylum %in% top102)  
write.csv(Y2avsd10, "Year2_AvSD_TOP10-UPDATED.csv")

## Year 3
Yr3 <- read.csv("Year3_RA.csv", row.names = 1)
#Sum by phylum across samples
physumY3 <- Yr3 %>% 
  group_by(Phylum) %>% 
  summarise(across(where(is.numeric), sum))
#Average phylum across samples
Y3mean <- apply(physumY3[,-1], 1, mean, na.rm=TRUE)
#Standard deviation across samples 
Y3std <- apply(physumY3[,-1], 1, sd, na.rm=TRUE)
#merge average and st.dev with rows
Y3avsd <- as.data.frame(cbind(physumY3$Phylum,Y3mean, Y3std))
#Renaming columns and saving as csv
colnames(Y3avsd)[1] ="Phylum"
colnames(Y3avsd)[2] ="Average"
colnames(Y3avsd)[3] ="Stand.Dev"
write.csv(Y3avsd, "Year3_AvSD-UPDATED.csv")
#Extract top 10 phyla and save as csv
top103 <- names(top10phy.names.Y3)
Y3avsd10 <- filter(Y3avsd,
                   Y3avsd$Phylum %in% top103)  
write.csv(Y3avsd10, "Year3_AvSD_TOP10-UPDATED.csv")

# Merge all years together and save as csv
#Original lists
#put all data frames into list
Y123avstd <- list(Y1avsd, Y2avsd, Y3avsd)
#merge all data frames in list
all <- Y123avstd %>% reduce(full_join, by='Phylum')
#renaming columns
colnames(all)[2] ="Y1mean"
colnames(all)[3] ="Y1std"
colnames(all)[4] ="Y2mean"
colnames(all)[5] ="Y2std"
colnames(all)[6] ="Y3mean"
colnames(all)[7] ="Y3std"

#Top 10 lists
Y123avstd10 <- list(Y1avsd10, Y2avsd10, Y3avsd10)
top10 <- Y123avstd10 %>% reduce(full_join, by='Phylum')
colnames(top10)[2] ="Y1mean"
colnames(top10)[3] ="Y1std"
colnames(top10)[4] ="Y2mean"
colnames(top10)[5] ="Y2std"
colnames(top10)[6] ="Y3mean"
colnames(top10)[7] ="Y3std"

#Save as csvs
write.csv(all, "Year123_AvSD.csv")
write.csv(top10, "Year123_AvSD_TOP10.csv")

###### Separating feature table by Station (CSVs) ######
CLV <- as.data.frame(t(dat.ra[grep("^CLV10A", rownames(dat.ra)),]))
KISS <- as.data.frame(t(dat.ra[grep("^KISSR0.0", rownames(dat.ra)),]))
L1 <- as.data.frame(t(dat.ra[grep("^L001", rownames(dat.ra)),]))
L4 <- as.data.frame(t(dat.ra[grep("^L004", rownames(dat.ra)),]))
L5 <- as.data.frame(t(dat.ra[grep("^L005", rownames(dat.ra)),]))
L6 <- as.data.frame(t(dat.ra[grep("^L006", rownames(dat.ra)),]))
L7 <- as.data.frame(t(dat.ra[grep("^L007", rownames(dat.ra)),]))
L8 <- as.data.frame(t(dat.ra[grep("^L008", rownames(dat.ra)),]))
LZ2 <- as.data.frame(t(dat.ra[grep("^LZ2_", rownames(dat.ra)),]))
Z25A <- as.data.frame(t(dat.ra[grep("^LZ25A", rownames(dat.ra)),]))
Z30 <- as.data.frame(t(dat.ra[grep("^LZ30", rownames(dat.ra)),]))
Z40 <- as.data.frame(t(dat.ra[grep("^LZ40", rownames(dat.ra)),]))
PALM <- as.data.frame(t(dat.ra[grep("^PALMOUT", rownames(dat.ra)),]))
PEL <- as.data.frame(t(dat.ra[grep("^PELBAY3", rownames(dat.ra)),]))
POLE3S <- as.data.frame(t(dat.ra[grep("^POLE3S", rownames(dat.ra)),]))
PO <- as.data.frame(t(dat.ra[grep("^POLESOUT", rownames(dat.ra)),]))
RIT <- as.data.frame(t(dat.ra[grep("^RITTAE2", rownames(dat.ra)),]))
S308 <- as.data.frame(t(dat.ra[grep("^S308", rownames(dat.ra)),]))
S77 <- as.data.frame(t(dat.ra[grep("^S77", rownames(dat.ra)),]))
S79 <- as.data.frame(t(dat.ra[grep("^S79", rownames(dat.ra)),]))

#S80 not included in adjusted dataset

###### Separating feature table by Year then Station (CSVs) ######
dat1 <- as.data.frame(t(dat.ra[grep("_19$", rownames(dat.ra)),]))
dat2 <- as.data.frame(t(dat.ra[grep("_20$", rownames(dat.ra)),]))
dat3 <- as.data.frame(t(dat.ra[grep("_21$", rownames(dat.ra)),]))
write.csv(dat1,"feature_Y1r_ADJUSTED.csv")
write.csv(dat2,"feature_Y2r_ADJUSTED.csv")
write.csv(dat3,"feature_Y3r_ADJUSTED.csv")
dat1 <- as.data.frame(t(dat1))
dat2 <- as.data.frame(t(dat2))
dat3 <- as.data.frame(t(dat3))

#Year 1 Stations
CLV <- as.data.frame(t(dat1[grep("^CLV10A", rownames(dat1)),]))
KISS <- as.data.frame(t(dat1[grep("^KISSR0.0", rownames(dat1)),]))
L1 <- as.data.frame(t(dat1[grep("^L001", rownames(dat1)),]))
L4 <- as.data.frame(t(dat1[grep("^L004", rownames(dat1)),]))
L5 <- as.data.frame(t(dat1[grep("^L005", rownames(dat1)),]))
L6 <- as.data.frame(t(dat1[grep("^L006", rownames(dat1)),]))
L7 <- as.data.frame(t(dat1[grep("^L007", rownames(dat1)),]))
L8 <- as.data.frame(t(dat1[grep("^L008", rownames(dat1)),]))
LZ2 <- as.data.frame(t(dat1[grep("^LZ2_", rownames(dat1)),]))
Z25A <- as.data.frame(t(dat1[grep("^LZ25A", rownames(dat1)),]))
Z30 <- as.data.frame(t(dat1[grep("^LZ30", rownames(dat1)),]))
Z40 <- as.data.frame(t(dat1[grep("^LZ40", rownames(dat1)),]))
PALM <- as.data.frame(t(dat1[grep("^PALMOUT", rownames(dat1)),]))
PEL <- as.data.frame(t(dat1[grep("^PELBAY3", rownames(dat1)),]))
POLE3S <- as.data.frame(t(dat1[grep("^POLE3S", rownames(dat1)),]))
PO <- as.data.frame(t(dat1[grep("^POLESOUT", rownames(dat1)),]))
RIT <- as.data.frame(t(dat1[grep("^RITTAE2", rownames(dat1)),]))
S308 <- as.data.frame(t(dat1[grep("^S308", rownames(dat1)),]))
S77 <- as.data.frame(t(dat1[grep("^S77", rownames(dat1)),]))
S79 <- as.data.frame(t(dat1[grep("^S79", rownames(dat1)),]))

#Year 2 Stations
CLV <- as.data.frame(t(dat2[grep("^CLV10A", rownames(dat2)),]))
KISS <- as.data.frame(t(dat2[grep("^KISSR0.0", rownames(dat2)),]))
L1 <- as.data.frame(t(dat2[grep("^L001", rownames(dat2)),]))
L4 <- as.data.frame(t(dat2[grep("^L004", rownames(dat2)),]))
L5 <- as.data.frame(t(dat2[grep("^L005", rownames(dat2)),]))
L6 <- as.data.frame(t(dat2[grep("^L006", rownames(dat2)),]))
L7 <- as.data.frame(t(dat2[grep("^L007", rownames(dat2)),]))
L8 <- as.data.frame(t(dat2[grep("^L008", rownames(dat2)),]))
LZ2 <- as.data.frame(t(dat2[grep("^LZ2_", rownames(dat2)),]))
Z25A <- as.data.frame(t(dat2[grep("^LZ25A", rownames(dat2)),]))
Z30 <- as.data.frame(t(dat2[grep("^LZ30", rownames(dat2)),]))
Z40 <- as.data.frame(t(dat2[grep("^LZ40", rownames(dat2)),]))
PALM <- as.data.frame(t(dat2[grep("^PALMOUT", rownames(dat2)),]))
PEL <- as.data.frame(t(dat2[grep("^PELBAY3", rownames(dat2)),]))
POLE3S <- as.data.frame(t(dat2[grep("^POLE3S", rownames(dat2)),]))
PO <- as.data.frame(t(dat2[grep("^POLESOUT", rownames(dat2)),]))
RIT <- as.data.frame(t(dat2[grep("^RITTAE2", rownames(dat2)),]))
S308 <- as.data.frame(t(dat2[grep("^S308", rownames(dat2)),]))
S77 <- as.data.frame(t(dat2[grep("^S77", rownames(dat2)),]))
S79 <- as.data.frame(t(dat2[grep("^S79", rownames(dat2)),]))

#Year 3 Stations
CLV <- as.data.frame(t(dat3[grep("^CLV10A", rownames(dat3)),]))
KISS <- as.data.frame(t(dat3[grep("^KISSR0.0", rownames(dat3)),]))
L1 <- as.data.frame(t(dat3[grep("^L001", rownames(dat3)),]))
L4 <- as.data.frame(t(dat3[grep("^L004", rownames(dat3)),]))
L5 <- as.data.frame(t(dat3[grep("^L005", rownames(dat3)),]))
L6 <- as.data.frame(t(dat3[grep("^L006", rownames(dat3)),]))
L7 <- as.data.frame(t(dat3[grep("^L007", rownames(dat3)),]))
L8 <- as.data.frame(t(dat3[grep("^L008", rownames(dat3)),]))
LZ2 <- as.data.frame(t(dat3[grep("^LZ2_", rownames(dat3)),]))
Z25A <- as.data.frame(t(dat3[grep("^LZ25A", rownames(dat3)),]))
Z30 <- as.data.frame(t(dat3[grep("^LZ30", rownames(dat3)),]))
Z40 <- as.data.frame(t(dat3[grep("^LZ40", rownames(dat3)),]))
PALM <- as.data.frame(t(dat3[grep("^PALMOUT", rownames(dat3)),]))
PEL <- as.data.frame(t(dat3[grep("^PELBAY3", rownames(dat3)),]))
POLE3S <- as.data.frame(t(dat3[grep("^POLE3S", rownames(dat3)),]))
PO <- as.data.frame(t(dat3[grep("^POLESOUT", rownames(dat3)),]))
RIT <- as.data.frame(t(dat3[grep("^RITTAE2", rownames(dat3)),]))
S308 <- as.data.frame(t(dat3[grep("^S308", rownames(dat3)),]))
S77 <- as.data.frame(t(dat3[grep("^S77", rownames(dat3)),]))
S79 <- as.data.frame(t(dat3[grep("^S79", rownames(dat3)),]))


###### TOP 10 TAXA BAR CHART - ALL YEARS TOGETHER ######
asvdat <- as.data.frame(t(dat.ra))
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
#Merging metadata, taxonomy, and ASV tables into one phyloseq object
physeq <- phyloseq(ASV,TAX,META)
#Use transform functions from microbiome package
transform <- microbiome::transform
#Merge rare taxa in to "Other"
physeq_transform <- transform(physeq, "compositional")
ASV  # 8,340 taxa & 541 samples
TAX  # 8,340 taxa by 7 tax. ranks
META # 541 samples by 42 sample variables

### Basic stats of seq. reads 
#Check number of microbes observed in each sample
sample_sums(physeq)
##Basic stats for reads of samples
sum(sample_sums(physeq))
#Total reads = 24,093,755
mean(sample_sums(physeq))
#Mean = 44,535.59
min(sample_sums(physeq))
#Min= 10,029
max(sample_sums(physeq))
#Max = 193,655
sd(sample_sums(physeq))
#Stan.Dev = 24,782.95
ntaxa(physeq)
#Total ASVs = 65,294

physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8340 taxa and 541 samples ]
# sample_data() Sample Data:       [ 541 samples by 42 sample variables ]
# tax_table()   Taxonomy Table:    [ 8340 taxa by 7 taxonomic ranks ]

##Retrieves the unique taxonomic ranks observed in the data set
##[#] = rank (starting from Domain and onward DPCOFGS)
get_taxa_unique(physeq, taxonomic.rank=rank_names(physeq)[7], errorIfNULL=TRUE)
#Unique Domains = 4
#Unique Phyla = 56
#Unique Classes = 142
#Unique Orders = 351
#Unique Families = 508
#Unique Genera = 728
#Unique Species = 317

## make sure there is a phyloseq object which includes the data, metadata, and taxonomy ##

## Aggregating by Taxa levels
phyPhy <- aggregate_taxa(physeq, 'Phylum')
phyClass <- aggregate_taxa(physeq, 'Class')
phyOrd <- aggregate_taxa(physeq, 'Order')
phyGen <- aggregate_taxa(physeq, 'Genus')
LakeOPhy <- as.data.frame(taxa_sums(phyPhy))
LakeOClass <- as.data.frame(taxa_sums(phyClass))
LakeOOrd <- as.data.frame(taxa_sums(phyOrd))
LakeOGenus <- as.data.frame(taxa_sums(phyGen))
#Saving each table as CSV
write.csv(LakeOPhy, "LakeOPhylaTotals.csv")
write.csv(LakeOClass, "LakeOClassesTotals.csv")
write.csv(LakeOOrd, "LakeOOrdersTotals.csv")
write.csv(LakeOGenus, "LakeOGeneraTotals.csv")

## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names <- sort(tapply(taxa_sums(physeq_transform), tax_table(physeq_transform)[, "Phylum"], sum), TRUE)[1:10]
## write.csv(top10phy.names, "Top10PhylaLakeO.csv")
# Proteobacteria       Bacteroidota      Cyanobacteria   Actinobacteriota  Verrucomicrobiota    Planctomycetota    Acidobacteriota 
# 121.550676         110.168874          81.682736          57.976055          38.301827          34.610471          15.164802 
# Bdellovibrionota        Chloroflexi    Gemmatimonadota 
# 14.615002          11.278973           9.640009 
#Cut down the physeq data to only the top 10 Phyla
top10phyla <- subset_taxa(physeq_transform, Phylum %in% names(top10phy.names))

## Plotting taxa stacked bar based on Zone
LakePhylaZ <- plot_bar(top10phyla, x="Zone", y="Abundance", fill="Phylum")
LakePhylaZ <- LakePhylaZ + 
  geom_bar(aes(fill=Phylum), stat="identity", position="fill", width = 0.96)+   #width=0.96 removes any space between bars
  ggtitle("Top 10 Phyla in Lake Okeechobee by Zone - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", 
             labeller = as_labeller(c('1'='Year 1 (2019)',
                                      '2'='Year 2 (2020)',
                                      '3'='Year 3 (2021)')))+   #scales=free -> allows ggplot to change the axes for the data shown in each facet
  theme_light()+                                                #labeller -> changing the labels of the grid
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
##Changing the color (by changing the default in ggplot2 [from HELP])
LakeOTop10 <- c("#2bcaf4","#24630e","#edc427","#1f60aa","#333333",
                         "#41ea27","#806bb4","#5f421b","#f08539","#ff9eed")
                         ## listed by phyla in alphabetical order
withr::with_options(list(ggplot2.discrete.fill = LakeOTop10, ggplot2.discrete.colour = LakeOTop10),print(LakePhylaZ))

###### Top 10 phyla each year (CSVs) ######
#Subsetting original ASV table by year
Y1r <- dat.ra[grep("_19$", rownames(dat.ra)),]
Y2r <- dat.ra[grep("_20$", rownames(dat.ra)),]
Y3r <- dat.ra[grep("_21$", rownames(dat.ra)),]
write.csv(t(Y1), "Year1_RA.csv")
write.csv(t(Y2), "Year2_RA.csv")
write.csv(t(Y3), "Year3_RA.csv")

# OR

#Load in data if already exported to CSVs
Y1r <- read.csv("Year1_RA.csv", row.names = 1)
Y2r <- read.csv("Year2_RA.csv", row.names = 1)
Y3r <- read.csv("Year3_RA.csv", row.names = 1)

#Top 10 phyla in each year
##Year 1
asvdat <- Y1r
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyY1<- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyY1_transform <- transform(phyY1, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.Y1 <- sort(tapply(taxa_sums(phyY1_transform), tax_table(phyY1_transform)[, "Phylum"], sum), TRUE)[1:10]
top10phy.names.Y1
# Proteobacteria       Bacteroidota      Cyanobacteria   Actinobacteriota    Planctomycetota  Verrucomicrobiota   Bdellovibrionota 
# 37.118712             34.048403          18.633005          16.562391          11.084878          10.877789           5.230602 
# Acidobacteriota        Chloroflexi      Crenarchaeota 
# 4.481436                3.249468           2.864711 
#Cut down the physeq data to only the top 10 Phyla
top10phylaY1 <- subset_taxa(phyY1_transform, Phylum %in% names(top10phy.names.Y1))
#Saving names and proportions as a data frame then saving as csv
topphylaY1 <- as.data.frame(top10phy.names.Y1)
colnames(topphylaY1)[1] ="Abundance"
write.csv(topphylaY1, "Top10Phyla_Year1.csv")

##Year 2
asvdat <- Y2r
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyY2<- phyloseq(ASV,TAX,META)
phyY2_transform <- transform(phyY2, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.Y2 <- sort(tapply(taxa_sums(phyY2_transform), tax_table(phyY2_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaY2 <- subset_taxa(phyY2_transform, Phylum %in% names(top10phy.names.Y2))
#Saving names and proportions as a data frame then saving as csv
topphylaY2 <- as.data.frame(top10phy.names.Y2)
colnames(topphylaY2)[1] ="Abundance"
write.csv(topphylaY2, "Top10Phyla_Year2.csv")

##Year 3
asvdat <- Y3r
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyY3<- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyY3_transform <- transform(phyY3, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.Y3 <- sort(tapply(taxa_sums(phyY3_transform), tax_table(phyY3_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaY3 <- subset_taxa(phyY3_transform, Phylum %in% names(top10phy.names.Y3))
#Saving names and proportions as a data frame then saving as csv
topphylaY3 <- as.data.frame(top10phy.names.Y3)
colnames(topphylaY3)[1] ="Abundance"
write.csv(topphylaY3, "Top10Phyla_Year3.csv")


###### Top 10 by Stations (CSVs) - ALL YEARS TOGETHER ######

## Use sample name order from Metadata file to keep samples in chronological order 
#Note: psmelt() turns phyloseq object into a large dataframe that is in LONG format

## CLV10A
asvdat <- CLV
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyCLV <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyCLV_transform <- transform(phyCLV, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.CLV <- sort(tapply(taxa_sums(phyCLV_transform), tax_table(phyCLV_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaCLV <- subset_taxa(phyCLV_transform, Phylum %in% names(top10phy.names.CLV))
#Saving names and proportions as a data frame then saving as csv
topphylaCLV <- as.data.frame(top10phy.names.CLV)
colnames(topphylaCLV)[1] ="Abundance"
write.csv(topphylaCLV, "Top10Phyla_CLV.csv")


## KISSR0.0 - (Firmicutes removed-> KISSR0.0_3_20)
asvdat <- KISS
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyKISS <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyKISS_transform <- transform(phyKISS, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.KISS <- sort(tapply(taxa_sums(phyKISS_transform), tax_table(phyKISS_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaKISS <- subset_taxa(phyKISS_transform, Phylum %in% names(top10phy.names.KISS))
#Saving names and proportions as a data frame then saving as csv
topphylaKISS <- as.data.frame(top10phy.names.KISS)
colnames(topphylaKISS)[1] ="Abundance"
write.csv(topphylaKISS, "Top10Phyla_KISS.csv")


## L001 
asvdat <- L1
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL1 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL1_transform <- transform(phyL1, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L1 <- sort(tapply(taxa_sums(phyL1_transform), tax_table(phyL1_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL1 <- subset_taxa(phyL1_transform, Phylum %in% names(top10phy.names.L1))
#Saving names and proportions as a data frame then saving as csv
topphylaL1 <- as.data.frame(top10phy.names.L1)
colnames(topphylaL1)[1] ="Abundance"
write.csv(topphylaL1, "Top10Phyla_L001.csv")


## L004 
asvdat <- L4
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL4 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL4_transform <- transform(phyL4, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L4 <- sort(tapply(taxa_sums(phyL4_transform), tax_table(phyL4_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL4 <- subset_taxa(phyL4_transform, Phylum %in% names(top10phy.names.L4))
#Saving names and proportions as a data frame then saving as csv
topphylaL4 <- as.data.frame(top10phy.names.L4)
colnames(topphylaL4)[1] ="Abundance"
write.csv(topphylaL4, "Top10Phyla_L004.csv")


## L005 (Firmicutes removed-> L005_3_20)
asvdat <- L5
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL5 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL5_transform <- transform(phyL5, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L5 <- sort(tapply(taxa_sums(phyL5_transform), tax_table(phyL5_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL5 <- subset_taxa(phyL5_transform, Phylum %in% names(top10phy.names.L5))
#Saving names and proportions as a data frame then saving as csv
topphylaL5 <- as.data.frame(top10phy.names.L5)
colnames(topphylaL5)[1] ="Abundance"
write.csv(topphylaL5, "Top10Phyla_L005.csv")


## L006
asvdat <- L6
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL6 <- phyloseq(ASV,TAX,META)
phyL6_transform <- transform(phyL6, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L6 <- sort(tapply(taxa_sums(phyL6_transform), tax_table(phyL6_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL6 <- subset_taxa(phyL6_transform, Phylum %in% names(top10phy.names.L6))
#Saving names and proportions as a data frame then saving as csv
topphylaL6 <- as.data.frame(top10phy.names.L6)
colnames(topphylaL6)[1] ="Abundance"
write.csv(topphylaL6, "Top10Phyla_L006.csv")


## L007
asvdat <- L7
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL7 <- phyloseq(ASV,TAX,META)
phyL7_transform <- transform(phyL7, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L7 <- sort(tapply(taxa_sums(phyL7_transform), tax_table(phyL7_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL7 <- subset_taxa(phyL7_transform, Phylum %in% names(top10phy.names.L7))
#Saving names and proportions as a data frame then saving as csv
topphylaL7 <- as.data.frame(top10phy.names.L7)
colnames(topphylaL7)[1] ="Abundance"
write.csv(topphylaL7, "Top10Phyla_L007.csv")

## L008
asvdat <- L8
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL8 <- phyloseq(ASV,TAX,META)
phyL8_transform <- transform(phyL8, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L8 <- sort(tapply(taxa_sums(phyL8_transform), tax_table(phyL8_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL8 <- subset_taxa(phyL8_transform, Phylum %in% names(top10phy.names.L8))
#Saving names and proportions as a data frame then saving as csv
topphylaL8 <- as.data.frame(top10phy.names.L8)
colnames(topphylaL8)[1] ="Abundance"
write.csv(topphylaL8, "Top10Phyla_L008.csv")

## LZ25A 
asvdat <- Z25A
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy25A <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phy25A_transform <- transform(phy25A, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.25A <- sort(tapply(taxa_sums(phy25A_transform), tax_table(phy25A_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla25A <- subset_taxa(phy25A_transform, Phylum %in% names(top10phy.names.25A))
#Saving names and proportions as a data frame then saving as csv
topphyla25A <- as.data.frame(top10phy.names.25A)
colnames(topphyla25A)[1] ="Abundance"
write.csv(topphyla25A, "Top10Phyla_LZ25A.csv")

## LZ2 (Firmicutes contam. removed LZ2_3_20)
asvdat <- LZ2
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyLZ2 <- phyloseq(ASV,TAX,META)
phyLZ2_transform <- transform(phyLZ2, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.LZ2 <- sort(tapply(taxa_sums(phyLZ2_transform), tax_table(phyLZ2_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaLZ2 <- subset_taxa(phyLZ2_transform, Phylum %in% names(top10phy.names.LZ2))
#Saving names and proportions as a data frame then saving as csv
topphylaLZ2 <- as.data.frame(top10phy.names.LZ2)
colnames(topphylaLZ2)[1] ="Abundance"
write.csv(topphylaLZ2, "Top10Phyla_LZ2.csv")

## LZ30
asvdat <- Z30
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy30 <- phyloseq(ASV,TAX,META)
phy30_transform <- transform(phy30, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.30 <- sort(tapply(taxa_sums(phy30_transform), tax_table(phy30_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla30 <- subset_taxa(phy30_transform, Phylum %in% names(top10phy.names.30))
#Saving names and proportions as a data frame then saving as csv
topphyla30 <- as.data.frame(top10phy.names.30)
colnames(topphyla30)[1] ="Abundance"
write.csv(topphyla30, "Top10Phyla_LZ30.csv")

## LZ40 
asvdat <- Z40
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy40 <- phyloseq(ASV,TAX,META)
phy40_transform <- transform(phy40, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.40 <- sort(tapply(taxa_sums(phy40_transform), tax_table(phy40_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla40 <- subset_taxa(phy40_transform, Phylum %in% names(top10phy.names.40))
#Saving names and proportions as a data frame then saving as csv
topphyla40 <- as.data.frame(top10phy.names.40)
colnames(topphyla40)[1] ="Abundance"
write.csv(topphyla40, "Top10Phyla_LZ40.csv")

## PALMOUT (Firmicutes contam. removed PALMOUT_3_20)
asvdat <- PALM
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPALM <- phyloseq(ASV,TAX,META)
phyPALM_transform <- transform(phyPALM, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PALM <- sort(tapply(taxa_sums(phyPALM_transform), tax_table(phyPALM_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPALM <- subset_taxa(phyPALM_transform, Phylum %in% names(top10phy.names.PALM))
#Saving names and proportions as a data frame then saving as csv
topphylaPALM <- as.data.frame(top10phy.names.PALM)
colnames(topphylaPALM)[1] ="Abundance"
write.csv(topphylaPALM, "Top10Phyla_PALM.csv")

## PELBAY3 - DONE ON 11/12/22
asvdat <- PEL
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPEL <- phyloseq(ASV,TAX,META)
phyPEL_transform <- transform(phyPEL, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PEL <- sort(tapply(taxa_sums(phyPEL_transform), tax_table(phyPEL_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPEL <- subset_taxa(phyPEL_transform, Phylum %in% names(top10phy.names.PEL))
#Saving names and proportions as a data frame then saving as csv
topphylaPEL <- as.data.frame(top10phy.names.PEL)
colnames(topphylaPEL)[1] ="Abundance"
write.csv(topphylaPEL, "Top10Phyla_PEL.csv")

## POLE3S - DONE ON 11/12/22 (Firmicutes contam. removed POLE3S_3_20)
asvdat <- POLE3S
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPOLE3S <- phyloseq(ASV,TAX,META)
phyPOLE3S_transform <- transform(phyPOLE3S, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.POLE3S <- sort(tapply(taxa_sums(phyPOLE3S_transform), tax_table(phyPOLE3S_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPOLE3S <- subset_taxa(phyPOLE3S_transform, Phylum %in% names(top10phy.names.POLE3S))
#Saving names and proportions as a data frame then saving as csv
topphylaPOLE3S <- as.data.frame(top10phy.names.POLE3S)
colnames(topphylaPOLE3S)[1] ="Abundance"
write.csv(topphylaPOLE3S, "Top10Phyla_POLE3S.csv")

## POLESOUT - DONE ON 11/12/22 (Firmicutes contam. removed POLESOUT_3_20)
asvdat <- PO
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPO <- phyloseq(ASV,TAX,META)
phyPO_transform <- transform(phyPO, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PO <- sort(tapply(taxa_sums(phyPO_transform), tax_table(phyPO_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPO <- subset_taxa(phyPO_transform, Phylum %in% names(top10phy.names.PO))
#Saving names and proportions as a data frame then saving as csv
topphylaPO <- as.data.frame(top10phy.names.PO)
colnames(topphylaPO)[1] ="Abundance"
write.csv(topphylaPO, "Top10Phyla_PO.csv")

## RITTAE2 - DONE ON 11/12/22 (Firmicutes contam. removed RITTAE2_3_20)
asvdat <- RIT
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyRIT <- phyloseq(ASV,TAX,META)
phyRIT_transform <- transform(phyRIT, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.RIT <- sort(tapply(taxa_sums(phyRIT_transform), tax_table(phyRIT_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaRIT <- subset_taxa(phyRIT_transform, Phylum %in% names(top10phy.names.RIT))
#Saving names and proportions as a data frame then saving as csv
topphylaRIT <- as.data.frame(top10phy.names.RIT)
colnames(topphylaRIT)[1] ="Abundance"
write.csv(topphylaRIT, "Top10Phyla_RIT.csv")

## S308 
asvdat <- S308
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS308 <- phyloseq(ASV,TAX,META)
phyS308_transform <- transform(phyS308, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S308 <- sort(tapply(taxa_sums(phyS308_transform), tax_table(phyS308_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS308 <- subset_taxa(phyS308_transform, Phylum %in% names(top10phy.names.S308))
#Saving names and proportions as a data frame then saving as csv
topphylaS308 <- as.data.frame(top10phy.names.S308)
colnames(topphylaS308)[1] ="Abundance"
write.csv(topphylaS308, "Top10Phyla_S308.csv")

## S77  (Firmicutes contam. removed S77_3_20)
asvdat <- S77
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS77 <- phyloseq(ASV,TAX,META)
phyS77_transform <- transform(phyS77, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S77 <- sort(tapply(taxa_sums(phyS77_transform), tax_table(phyS77_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS77 <- subset_taxa(phyS77_transform, Phylum %in% names(top10phy.names.S77))
#Saving names and proportions as a data frame then saving as csv
topphylaS77 <- as.data.frame(top10phy.names.S77)
colnames(topphylaS77)[1] ="Abundance"
write.csv(topphylaS77, "Top10Phyla_S77.csv")

## S79 (Firmicutes contam. removed S79_3_20)
asvdat <- S79
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS79 <- phyloseq(ASV,TAX,META)
phyS79_transform <- transform(phyS79, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S79 <- sort(tapply(taxa_sums(phyS79_transform), tax_table(phyS79_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS79 <- subset_taxa(phyS79_transform, Phylum %in% names(top10phy.names.S79))
#Saving names and proportions as a data frame then saving as csv
topphylaS79 <- as.data.frame(top10phy.names.S79)
colnames(topphylaS79)[1] ="Abundance"
write.csv(topphylaS79, "Top10Phyla_S79.csv")

###### Plotting Taxonomy Bar plots using phyloseq - ALL YEARS TOGETHER ######
#Defining the initial plot 
CLV <- plot_bar(top10phylaCLV, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
CLV$data$Sample <- as.factor(CLV$data$Sample) #Assigning the samples as factors so I can manually put the levels in order
levels(CLV$data$Sample) #making sure each sample name is a level (should be 28 levels)
#Samples ARE NOT in chronological order here
CLV$data$Sample <- factor(CLV$data$Sample, levels=c("CLV10A_4_19","CLV10A_5_19","CLV10A_6_19","CLV10A_7_19","CLV10A_8_19",
                                                    "CLV10A_9_19","CLV10A_10_19","CLV10A_11_19","CLV10A_12_19","CLV10A_1_20",
                                                    "CLV10A_2_20","CLV10A_3_20","CLV10A_4_20","CLV10A_6_20","CLV10A_7_20",
                                                    "CLV10A_8_20","CLV10A_9_20","CLV10A_10_20","CLV10A_12_20","CLV10A_1_21",
                                                    "CLV10A_2_21","CLV10A_3_21","CLV10A_4_21","CLV10A_5_21","CLV10A_6_21",
                                                    "CLV10A_7_21","CLV10A_8_21","CLV10A_10_21"))
levels(CLV$data$Sample) #Samples ARE in chronological order now 
#Customizing the plot using ggplot2's geom_bar
CLV + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.96)+   #width=0.96 removes any space between bars
  ggtitle("Top 10 Phyla at CLV10A - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", 
             labeller = as_labeller(c('1'='Year 1 (2019)',
                                      '2'='Year 2 (2020)',
                                      '3'='Year 3 (2021)')))+   #scales=free -> allows ggplot to change the axes for the data shown in each facet
  theme_light()+                                                #labeller -> changing the labels of the grid
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
#facet_grid - splits up the graph into the variable specified
#position=fill - bars go up to 1.00, while position=stack - bar shows actual abundance (bars don't line up) 

#Defining the initial plot 
KISS <- plot_bar(top10phylaKISS, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
KISS$data$Sample <- as.factor(KISS$data$Sample)
levels(KISS$data$Sample)
KISS$data$Sample <- factor(KISS$data$Sample, levels=c("KISSR0.0_3_19","KISSR0.0_4_19","KISSR0.0_5_19","KISSR0.0_7_19","KISSR0.0_8_19","KISSR0.0_9_19",
                                                      "KISSR0.0_11_19","KISSR0.0_12_19","KISSR0.0_1_20","KISSR0.0_2_20","KISSR0.0_4_20",
                                                      "KISSR0.0_5_20","KISSR0.0_6_20","KISSR0.0_8_20","KISSR0.0_9_20","KISSR0.0_10_20","KISSR0.0_11_20",
                                                      "KISSR0.0_12_20","KISSR0.0_2_21","KISSR0.0_3_21","KISSR0.0_4_21","KISSR0.0_5_21","KISSR0.0_6_21",
                                                      "KISSR0.0_7_21","KISSR0.0_8_21","KISSR0.0_9_21","KISSR0.0_10_21"))
levels(KISS$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
KISS + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+ 
  ggtitle("Top 10 Phyla at KISSR0.0 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

#Defining the initial plot 
L1 <- plot_bar(top10phylaL1, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L1$data$Sample <- as.factor(L1$data$Sample)
levels(L1$data$Sample)
L1$data$Sample <- factor(L1$data$Sample, levels=c("L001_3_19","L001_4_19","L001_5_19","L001_6_19","L001_7_19","L001_8_19","L001_9_19",
                                                  "L001_11_19","L001_12_19","L001_1_20","L001_2_20","L001_3_20","L001_4_20",
                                                  "L001_6_20","L001_7_20","L001_8_20","L001_9_20","L001_10_20","L001_11_20",
                                                  "L001_12_20","L001_2_21","L001_3_21","L001_4_21","L001_5_21","L001_6_21",
                                                  "L001_7_21","L001_8_21","L001_9_21","L001_10_21"))
levels(L1$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L1 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L001 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

#Defining the initial plot 
L4 <- plot_bar(top10phylaL4, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L4$data$Sample <- as.factor(L4$data$Sample)
levels(L4$data$Sample)
L4$data$Sample <- factor(L4$data$Sample, levels=c("L004_3_19","L004_5_19","L004_8_19","L004_9_19",
                                                  "L004_11_19","L004_12_19","L004_1_20","L004_2_20","L004_3_20","L004_4_20",
                                                  "L004_6_20","L004_7_20","L004_8_20","L004_9_20","L004_10_20","L004_11_20",
                                                  "L004_12_20","L004_2_21","L004_3_21","L004_4_21","L004_6_21",
                                                  "L004_7_21","L004_8_21","L004_10_21"))
levels(L4$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L4 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L004 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
## Top 10 Classes - 12/01/22
#Sort Class by abundance and pick the top 10
top10class.names.L4 <- sort(tapply(taxa_sums(phyL4_transform), tax_table(phyL4_transform)[, "Class"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 classes
top10classL4 <- subset_taxa(phyL4_transform, Class %in% names(top10class.names.L4))
#Saving names and proportions as a data frame then saving as csv
topclassL4 <- as.data.frame(top10class.names.L4)
colnames(topclassL4)[1] ="Abundance"
write.csv(topclassL4, "Top10Classes_L004.csv")
### Plotting the graph -PHYLUM
#Defining the initial plot 
L4c <- plot_bar(top10classL4, x="Sample", y="Abundance", fill = "Class")
#Reordering the samples so they plot in chronological order
L4c$data$Sample <- as.factor(L4c$data$Sample)
levels(L4c$data$Sample)
L4c$data$Sample <- factor(L4c$data$Sample, levels=c("L004_3_19","L004_5_19","L004_8_19","L004_9_19",
                                                    "L004_11_19","L004_12_19","L004_1_20","L004_2_20","L004_3_20","L004_4_20",
                                                    "L004_6_20","L004_7_20","L004_8_20","L004_9_20","L004_10_20","L004_11_20",
                                                    "L004_12_20","L004_2_21","L004_3_21","L004_4_21","L004_6_21",
                                                    "L004_7_21","L004_8_21","L004_10_21"))
levels(L4c$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L4c + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Classes at L004 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

#Defining the initial plot 
L5 <- plot_bar(top10phylaL5, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L5$data$Sample <- as.factor(L5$data$Sample)
levels(L5$data$Sample)
L5$data$Sample <- factor(L5$data$Sample, levels=c("L005_3_19","L005_4_19","L005_5_19","L005_6_19","L005_7_19","L005_8_19","L005_9_19",
                                                  "L005_11_19","L005_12_19","L005_1_20","L005_2_20","L005_4_20",
                                                  "L005_6_20","L005_7_20","L005_8_20","L005_9_20","L005_10_20","L005_11_20",
                                                  "L005_12_20","L005_2_21","L005_3_21","L005_4_21","L005_5_21","L005_6_21",
                                                  "L005_7_21","L005_8_21","L005_9_21","L005_10_21"))
levels(L5$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L5 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L005 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

#Defining the initial plot 
L6 <- plot_bar(top10phylaL6, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L6$data$Sample <- as.factor(L6$data$Sample)
levels(L6$data$Sample)
L6$data$Sample <- factor(L6$data$Sample, levels=c("L006_5_19","L006_7_19","L006_8_19","L006_9_19",
                                                  "L006_11_19","L006_12_19","L006_1_20","L006_2_20","L006_3_20","L006_4_20",
                                                  "L006_5_20","L006_6_20","L006_7_20","L006_8_20","L006_9_20","L006_10_20","L006_11_20",
                                                  "L006_12_20","L006_1_21","L006_2_21","L006_3_21","L006_4_21","L006_5_21","L006_6_21",
                                                  "L006_7_21","L006_8_21","L006_10_21"))
levels(L6$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L6 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L006 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
L7 <- plot_bar(top10phylaL7, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L7$data$Sample <- as.factor(L7$data$Sample)
levels(L7$data$Sample)
L7$data$Sample <- factor(L7$data$Sample, levels=c("L007_3_19","L007_4_19","L007_5_19","L007_6_19","L007_7_19","L007_8_19","L007_9_19",
                                                  "L007_11_19","L007_12_19","L007_1_20","L007_2_20","L007_3_20","L007_4_20","L007_5_20",
                                                  "L007_6_20","L007_8_20","L007_9_20","L007_10_20","L007_11_20",
                                                  "L007_12_20","L007_1_21","L007_2_21","L007_3_21","L007_4_21","L007_5_21","L007_6_21",
                                                  "L007_7_21","L007_8_21","L007_9_21","L007_10_21"))
levels(L7$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
L7 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L007 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
L8 <- plot_bar(top10phylaL8, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
L8$data$Sample <- as.factor(L8$data$Sample)
levels(L8$data$Sample)
L8$data$Sample <- factor(L8$data$Sample, levels=c("L008_3_19","L008_5_19","L008_6_19","L008_7_19","L008_8_19","L008_9_19",
                                                  "L008_11_19","L008_12_19","L008_1_20","L008_2_20","L008_3_20","L008_4_20","L008_5_20",
                                                  "L008_6_20","L008_7_20","L008_8_20","L008_9_20","L008_10_20","L008_11_20",
                                                  "L008_12_20","L008_2_21","L008_3_21","L008_4_21","L008_5_21","L008_6_21",
                                                  "L008_7_21","L008_8_21","L008_10_21"))
levels(L8$data$Sample)
#Customizing the plot using ggplot2's geom_bar
L8 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at L008 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
Z25A <- plot_bar(top10phyla25A, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
Z25A$data$Sample <- as.factor(Z25A$data$Sample)
levels(Z25A$data$Sample)
Z25A$data$Sample <- factor(Z25A$data$Sample, levels=c("LZ25A_3_19","LZ25A_4_19","LZ25A_6_19","LZ25A_7_19","LZ25A_8_19","LZ25A_9_19",
                                                      "LZ25A_11_19","LZ25A_12_19","LZ25A_1_20","LZ25A_2_20","LZ25A_3_20","LZ25A_4_20",
                                                      "LZ25A_5_20","LZ25A_7_20","LZ25A_8_20","LZ25A_9_20","LZ25A_10_20","LZ25A_11_20",
                                                      "LZ25A_12_20","LZ25A_1_21","LZ25A_2_21","LZ25A_3_21","LZ25A_4_21","LZ25A_5_21","LZ25A_6_21",
                                                      "LZ25A_7_21","LZ25A_10_21"))
levels(Z25A$data$Sample)
#Customizing the plot using ggplot2's geom_bar
Z25A + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at LZ25A - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
LZ2 <- plot_bar(top10phylaLZ2, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
LZ2$data$Sample <- as.factor(LZ2$data$Sample)
levels(LZ2$data$Sample)
LZ2$data$Sample <- factor(LZ2$data$Sample, levels=c("LZ2_3_19","LZ2_4_19","LZ2_5_19","LZ2_6_19","LZ2_8_19","LZ2_9_19",
                                                    "LZ2_11_19","LZ2_12_19","LZ2_1_20","LZ2_2_20","LZ2_4_20",
                                                    "LZ2_5_20","LZ2_6_20","LZ2_7_20","LZ2_8_20","LZ2_9_20","LZ2_10_20","LZ2_11_20",
                                                    "LZ2_12_20","LZ2_2_21","LZ2_3_21","LZ2_4_21","LZ2_5_21","LZ2_6_21",
                                                    "LZ2_7_21","LZ2_8_21","LZ2_10_21"))
levels(LZ2$data$Sample)
#Customizing the plot using ggplot2's geom_bar
LZ2 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at LZ2 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))


### Plotting the graph
#Defining the initial plot 
Z30 <- plot_bar(top10phyla30, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
Z30$data$Sample <- as.factor(Z30$data$Sample)
levels(Z30$data$Sample)
Z30$data$Sample <- factor(Z30$data$Sample, levels=c("LZ30_4_19","LZ30_5_19","LZ30_6_19","LZ30_7_19","LZ30_8_19","LZ30_9_19","LZ30_10_19",
                                                    "LZ30_11_19","LZ30_12_19","LZ30_1_20","LZ30_2_20","LZ30_3_20","LZ30_4_20","LZ30_5_20",
                                                    "LZ30_6_20","LZ30_7_20","LZ30_8_20","LZ30_9_20","LZ30_10_20","LZ30_11_20",
                                                    "LZ30_12_20","LZ30_1_21","LZ30_2_21","LZ30_3_21","LZ30_4_21","LZ30_5_21","LZ30_6_21",
                                                    "LZ30_7_21","LZ30_8_21","LZ30_10_21"))
levels(Z30$data$Sample)
#Customizing the plot using ggplot2's geom_bar
Z30 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at LZ30 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
Z40 <- plot_bar(top10phyla40, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
Z40$data$Sample <- as.factor(Z40$data$Sample)
levels(Z40$data$Sample)
Z40$data$Sample <- factor(Z40$data$Sample, levels=c("LZ40_3_19","LZ40_4_19","LZ40_5_19","LZ40_6_19","LZ40_7_19","LZ40_8_19","LZ40_9_19",
                                                    "LZ40_11_19","LZ40_12_19","LZ40_1_20","LZ40_2_20","LZ40_3_20","LZ40_4_20",
                                                    "LZ40_5_20","LZ40_6_20","LZ40_7_20","LZ40_8_20","LZ40_9_20","LZ40_10_20",
                                                    "LZ40_12_20","LZ40_1_21","LZ40_2_21","LZ40_3_21","LZ40_4_21","LZ40_5_21","LZ40_6_21",
                                                    "LZ40_7_21","LZ40_8_21","LZ40_9_21","LZ40_10_21"))
levels(Z40$data$Sample)
#Customizing the plot using ggplot2's geom_bar
Z40 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at LZ40 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
PALM <- plot_bar(top10phylaPALM, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
PALM$data$Sample <- as.factor(PALM$data$Sample)
levels(PALM$data$Sample)
PALM$data$Sample <- factor(PALM$data$Sample, levels=c("PALMOUT_3_19","PALMOUT_4_19","PALMOUT_6_19","PALMOUT_7_19","PALMOUT_8_19",
                                                      "PALMOUT_11_19","PALMOUT_12_19","PALMOUT_1_20","PALMOUT_2_20","PALMOUT_4_20",
                                                      "PALMOUT_5_20","PALMOUT_6_20","PALMOUT_7_20","PALMOUT_8_20","PALMOUT_9_20","PALMOUT_10_20",
                                                      "PALMOUT_12_20","PALMOUT_1_21","PALMOUT_2_21","PALMOUT_3_21","PALMOUT_4_21","PALMOUT_5_21","PALMOUT_6_21",
                                                      "PALMOUT_7_21","PALMOUT_8_21","PALMOUT_9_21","PALMOUT_10_21"))
levels(PALM$data$Sample)

#Customizing the plot using ggplot2's geom_bar
PALM + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at PALMOUT - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
PEL <- plot_bar(top10phylaPEL, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
PEL$data$Sample <- as.factor(PEL$data$Sample)
levels(PEL$data$Sample)
PEL$data$Sample <- factor(PEL$data$Sample, levels=c("PELBAY3_3_19","PELBAY3_5_19","PELBAY3_6_19","PELBAY3_7_19","PELBAY3_8_19","PELBAY3_9_19",
                                                    "PELBAY3_11_19","PELBAY3_12_19","PELBAY3_1_20","PELBAY3_2_20","PELBAY3_4_20","PELBAY3_5_20",
                                                    "PELBAY3_6_20","PELBAY3_7_20","PELBAY3_8_20","PELBAY3_9_20","PELBAY3_10_20","PELBAY3_11_20",
                                                    "PELBAY3_12_20","PELBAY3_1_21","PELBAY3_2_21","PELBAY3_3_21","PELBAY3_4_21","PELBAY3_5_21","PELBAY3_6_21",
                                                    "PELBAY3_7_21","PELBAY3_8_21","PELBAY3_10_21"))
levels(PEL$data$Sample) 
#Customizing the plot using ggplot2's geom_bar
PEL + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at PELBAY3 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))

### Plotting the graph
#Defining the initial plot 
POLE3S <- plot_bar(top10phylaPOLE3S, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
POLE3S$data$Sample <- as.factor(POLE3S$data$Sample)
levels(POLE3S$data$Sample)
POLE3S$data$Sample <- factor(POLE3S$data$Sample, levels=c("POLE3S_3_19","POLE3S_5_19","POLE3S_6_19","POLE3S_7_19","POLE3S_8_19",
                                                          "POLE3S_12_19","POLE3S_1_20","POLE3S_2_20","POLE3S_4_20",
                                                          "POLE3S_7_20","POLE3S_8_20","POLE3S_9_20","POLE3S_10_20","POLE3S_11_20",
                                                          "POLE3S_12_20","POLE3S_1_21","POLE3S_2_21","POLE3S_3_21","POLE3S_4_21","POLE3S_5_21","POLE3S_6_21",
                                                          "POLE3S_7_21","POLE3S_8_21","POLE3S_10_21"))
levels(POLE3S$data$Sample) 
#Customizing the plot using ggplot2's geom_bar and exporting as PNG file
png("Top10PhylaPOLE3S.png", width = 885, height = 575) # creates a named png file in your working directory
POLE3S + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at POLE3S - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off() #stops writing to the png file and saves it

### Plotting the graph
#Defining the initial plot 
PO <- plot_bar(top10phylaPO, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
PO$data$Sample <- as.factor(PO$data$Sample)
levels(PO$data$Sample)
PO$data$Sample <- factor(PO$data$Sample, levels=c("POLESOUT_3_19","POLESOUT_4_19","POLESOUT_5_19","POLESOUT_6_19","POLESOUT_7_19","POLESOUT_8_19",
                                                  "POLESOUT_11_19","POLESOUT_1_20","POLESOUT_2_20","POLESOUT_4_20",
                                                  "POLESOUT_6_20","POLESOUT_7_20","POLESOUT_8_20","POLESOUT_9_20","POLESOUT_10_20","POLESOUT_11_20",
                                                  "POLESOUT_12_20","POLESOUT_2_21","POLESOUT_3_21","POLESOUT_4_21","POLESOUT_5_21","POLESOUT_6_21",
                                                  "POLESOUT_7_21","POLESOUT_8_21","POLESOUT_9_21","POLESOUT_10_21"))
levels(PO$data$Sample)
#Customizing the plot using ggplot2's geom_bar
png("Top10PhylaPOLESOUT.png", width = 885, height = 575)
PO + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at POLESOUT - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off()

### Plotting the graph
#Defining the initial plot 
RIT <- plot_bar(top10phylaRIT, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
RIT$data$Sample <- as.factor(RIT$data$Sample)
levels(RIT$data$Sample)
RIT$data$Sample <- factor(RIT$data$Sample, levels=c("RITTAE2_3_19","RITTAE2_4_19","RITTAE2_6_19","RITTAE2_7_19","RITTAE2_8_19",
                                                    "RITTAE2_11_19","RITTAE2_12_19","RITTAE2_1_20","RITTAE2_2_20","RITTAE2_4_20",
                                                    "RITTAE2_8_20","RITTAE2_9_20","RITTAE2_10_20","RITTAE2_11_20",
                                                    "RITTAE2_12_20","RITTAE2_1_21","RITTAE2_2_21","RITTAE2_3_21","RITTAE2_4_21","RITTAE2_5_21","RITTAE2_6_21",
                                                    "RITTAE2_7_21","RITTAE2_8_21","RITTAE2_10_21"))
levels(RIT$data$Sample)
#Customizing the plot using ggplot2's geom_bar
png("Top10PhylaRITTAE2.png", width = 885, height = 575)
RIT + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at RITTAE2 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off()

### Plotting the graph
#Defining the initial plot 
S308 <- plot_bar(top10phylaS308, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
S308$data$Sample <- as.factor(S308$data$Sample)
levels(S308$data$Sample)
S308$data$Sample <- factor(S308$data$Sample, levels=c("S308_3_19","S308_4_19","S308_5_19","S308_6_19","S308_7_19","S308_9_19","S308_10_19",
                                                      "S308_11_19","S308_12_19","S308_1_20","S308_2_20","S308_3_20","S308_4_20","S308_5_20",
                                                      "S308_6_20","S308_7_20","S308_8_20","S308_9_20","S308_10_20","S308_11_20",
                                                      "S308_12_20","S308_1_21","S308_2_21","S308_3_21","S308_4_21","S308_5_21","S308_6_21"))
levels(S308$data$Sample)
#Customizing the plot using ggplot2's geom_bar
png("Top10PhylaS308.png", width = 885, height = 575)
S308 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at S308 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off()

### Plotting the graph
#Defining the initial plot 
S77 <- plot_bar(top10phylaS77, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
S77$data$Sample <- as.factor(S77$data$Sample)
levels(S77$data$Sample)
S77$data$Sample <- factor(S77$data$Sample, levels=c("S77_3_19","S77_4_19","S77_5_19","S77_6_19","S77_7_19","S77_8_19","S77_10_19",
                                                    "S77_11_19","S77_12_19","S77_1_20","S77_2_20","S77_4_20",
                                                    "S77_6_20","S77_7_20","S77_8_20","S77_9_20","S77_10_20",
                                                    "S77_12_20","S77_2_21","S77_3_21","S77_4_21","S77_5_21","S77_6_21",
                                                    "S77_7_21","S77_8_21","S77_9_21","S77_10_21"))
levels(S77$data$Sample)
#Customizing the plot using ggplot2's geom_bar
png("Top10PhylaS77.png", width = 885, height = 575)
S77 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at S77 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off()

### Plotting the graph
#Defining the initial plot 
S79 <- plot_bar(top10phylaS79, x="Sample", y="Abundance", fill = "Phylum")
#Reordering the samples so they plot in chronological order
S79$data$Sample <- as.factor(S79$data$Sample)
levels(S79$data$Sample)
S79$data$Sample <- factor(S79$data$Sample, levels=c("S79_3_19","S79_4_19","S79_6_19","S79_7_19","S79_8_19","S79_12_19",
                                                    "S79_1_20","S79_2_20","S79_4_20",
                                                    "S79_7_20","S79_9_20","S79_10_20","S79_11_20",
                                                    "S79_12_20","S79_2_21","S79_3_21","S79_4_21","S79_5_21","S79_6_21",
                                                    "S79_7_21","S79_8_21","S79_10_21"))
levels(S79$data$Sample)
#Customizing the plot using ggplot2's geom_bar
png("Top10PhylaS79.png", width = 885, height = 575)
S79 + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="fill", width = 0.98)+
  ggtitle("Top 10 Phyla at S79 - March 2019 to October 2021")+
  facet_grid(.~Year, scales = "free", labeller = as_labeller(c('1'='Year 1 (2019)',
                                                               '2'='Year 2 (2020)',
                                                               '3'='Year 3 (2021)')))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
dev.off()


###### Top 10 by Stations each Year - exporting CSVs ######
#### Year 1
## CLV10A
asvdat <- CLV
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyCLV <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyCLV_transform <- transform(phyCLV, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.CLV <- sort(tapply(taxa_sums(phyCLV_transform), tax_table(phyCLV_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaCLV <- subset_taxa(phyCLV_transform, Phylum %in% names(top10phy.names.CLV))
#Saving names and proportions as a data frame then saving as csv
topphylaCLV <- as.data.frame(top10phy.names.CLV)
colnames(topphylaCLV)[1] ="Abundance"
write.csv(topphylaCLV, "Top10Phyla_CLV_Y1.csv")

## KISSR0.0 - (Firmicutes removed-> KISSR0.0_3_20)
asvdat <- KISS
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyKISS <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyKISS_transform <- transform(phyKISS, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.KISS <- sort(tapply(taxa_sums(phyKISS_transform), tax_table(phyKISS_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaKISS <- subset_taxa(phyKISS_transform, Phylum %in% names(top10phy.names.KISS))
#Saving names and proportions as a data frame then saving as csv
topphylaKISS <- as.data.frame(top10phy.names.KISS)
colnames(topphylaKISS)[1] ="Abundance"
write.csv(topphylaKISS, "Top10Phyla_KISS_Y1.csv")

## L001 
asvdat <- L1
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL1 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL1_transform <- transform(phyL1, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L1 <- sort(tapply(taxa_sums(phyL1_transform), tax_table(phyL1_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL1 <- subset_taxa(phyL1_transform, Phylum %in% names(top10phy.names.L1))
#Saving names and proportions as a data frame then saving as csv
topphylaL1 <- as.data.frame(top10phy.names.L1)
colnames(topphylaL1)[1] ="Abundance"
write.csv(topphylaL1, "Top10Phyla_L001_Y1.csv")

## L004 
asvdat <- L4
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL4 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL4_transform <- transform(phyL4, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L4 <- sort(tapply(taxa_sums(phyL4_transform), tax_table(phyL4_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL4 <- subset_taxa(phyL4_transform, Phylum %in% names(top10phy.names.L4))
#Saving names and proportions as a data frame then saving as csv
topphylaL4 <- as.data.frame(top10phy.names.L4)
colnames(topphylaL4)[1] ="Abundance"
write.csv(topphylaL4, "Top10Phyla_L004_Y1.csv")

## L005 (Firmicutes removed-> L005_3_20)
asvdat <- L5
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL5 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL5_transform <- transform(phyL5, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L5 <- sort(tapply(taxa_sums(phyL5_transform), tax_table(phyL5_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL5 <- subset_taxa(phyL5_transform, Phylum %in% names(top10phy.names.L5))
#Saving names and proportions as a data frame then saving as csv
topphylaL5 <- as.data.frame(top10phy.names.L5)
colnames(topphylaL5)[1] ="Abundance"
write.csv(topphylaL5, "Top10Phyla_L005_Y1.csv")

## L006
asvdat <- L6
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL6 <- phyloseq(ASV,TAX,META)
phyL6_transform <- transform(phyL6, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L6 <- sort(tapply(taxa_sums(phyL6_transform), tax_table(phyL6_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL6 <- subset_taxa(phyL6_transform, Phylum %in% names(top10phy.names.L6))
#Saving names and proportions as a data frame then saving as csv
topphylaL6 <- as.data.frame(top10phy.names.L6)
colnames(topphylaL6)[1] ="Abundance"
write.csv(topphylaL6, "Top10Phyla_L006_Y1.csv")

## L007
asvdat <- L7
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL7 <- phyloseq(ASV,TAX,META)
phyL7_transform <- transform(phyL7, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L7 <- sort(tapply(taxa_sums(phyL7_transform), tax_table(phyL7_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL7 <- subset_taxa(phyL7_transform, Phylum %in% names(top10phy.names.L7))
#Saving names and proportions as a data frame then saving as csv
topphylaL7 <- as.data.frame(top10phy.names.L7)
colnames(topphylaL7)[1] ="Abundance"
write.csv(topphylaL7, "Top10Phyla_L007_Y1.csv")

## L008
asvdat <- L8
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL8 <- phyloseq(ASV,TAX,META)
phyL8_transform <- transform(phyL8, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L8 <- sort(tapply(taxa_sums(phyL8_transform), tax_table(phyL8_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL8 <- subset_taxa(phyL8_transform, Phylum %in% names(top10phy.names.L8))
#Saving names and proportions as a data frame then saving as csv
topphylaL8 <- as.data.frame(top10phy.names.L8)
colnames(topphylaL8)[1] ="Abundance"
write.csv(topphylaL8, "Top10Phyla_L008_Y1.csv")

## LZ25A 
asvdat <- Z25A
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy25A <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phy25A_transform <- transform(phy25A, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.25A <- sort(tapply(taxa_sums(phy25A_transform), tax_table(phy25A_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla25A <- subset_taxa(phy25A_transform, Phylum %in% names(top10phy.names.25A))
#Saving names and proportions as a data frame then saving as csv
topphyla25A <- as.data.frame(top10phy.names.25A)
colnames(topphyla25A)[1] ="Abundance"
write.csv(topphyla25A, "Top10Phyla_LZ25A_Y1.csv")

## LZ2 (Firmicutes contam. removed LZ2_3_20)
asvdat <- LZ2
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyLZ2 <- phyloseq(ASV,TAX,META)
phyLZ2_transform <- transform(phyLZ2, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.LZ2 <- sort(tapply(taxa_sums(phyLZ2_transform), tax_table(phyLZ2_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaLZ2 <- subset_taxa(phyLZ2_transform, Phylum %in% names(top10phy.names.LZ2))
#Saving names and proportions as a data frame then saving as csv
topphylaLZ2 <- as.data.frame(top10phy.names.LZ2)
colnames(topphylaLZ2)[1] ="Abundance"
write.csv(topphylaLZ2, "Top10Phyla_LZ2_Y1.csv")

## LZ30
asvdat <- Z30
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy30 <- phyloseq(ASV,TAX,META)
phy30_transform <- transform(phy30, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.30 <- sort(tapply(taxa_sums(phy30_transform), tax_table(phy30_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla30 <- subset_taxa(phy30_transform, Phylum %in% names(top10phy.names.30))
#Saving names and proportions as a data frame then saving as csv
topphyla30 <- as.data.frame(top10phy.names.30)
colnames(topphyla30)[1] ="Abundance"
write.csv(topphyla30, "Top10Phyla_LZ30_Y1.csv")

## LZ40 
asvdat <- Z40
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy40 <- phyloseq(ASV,TAX,META)
phy40_transform <- transform(phy40, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.40 <- sort(tapply(taxa_sums(phy40_transform), tax_table(phy40_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla40 <- subset_taxa(phy40_transform, Phylum %in% names(top10phy.names.40))
#Saving names and proportions as a data frame then saving as csv
topphyla40 <- as.data.frame(top10phy.names.40)
colnames(topphyla40)[1] ="Abundance"
write.csv(topphyla40, "Top10Phyla_LZ40_Y1.csv")

## PALMOUT (Firmicutes contam. removed PALMOUT_3_20)
asvdat <- PALM
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPALM <- phyloseq(ASV,TAX,META)
phyPALM_transform <- transform(phyPALM, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PALM <- sort(tapply(taxa_sums(phyPALM_transform), tax_table(phyPALM_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPALM <- subset_taxa(phyPALM_transform, Phylum %in% names(top10phy.names.PALM))
#Saving names and proportions as a data frame then saving as csv
topphylaPALM <- as.data.frame(top10phy.names.PALM)
colnames(topphylaPALM)[1] ="Abundance"
write.csv(topphylaPALM, "Top10Phyla_PALM_Y1.csv")

## PELBAY3 - DONE ON 11/12/22
asvdat <- PEL
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPEL <- phyloseq(ASV,TAX,META)
phyPEL_transform <- transform(phyPEL, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PEL <- sort(tapply(taxa_sums(phyPEL_transform), tax_table(phyPEL_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPEL <- subset_taxa(phyPEL_transform, Phylum %in% names(top10phy.names.PEL))
#Saving names and proportions as a data frame then saving as csv
topphylaPEL <- as.data.frame(top10phy.names.PEL)
colnames(topphylaPEL)[1] ="Abundance"
write.csv(topphylaPEL, "Top10Phyla_PEL_Y1.csv")

## POLE3S (Firmicutes contam. removed POLE3S_3_20)
asvdat <- POLE3S
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPOLE3S <- phyloseq(ASV,TAX,META)
phyPOLE3S_transform <- transform(phyPOLE3S, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.POLE3S <- sort(tapply(taxa_sums(phyPOLE3S_transform), tax_table(phyPOLE3S_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPOLE3S <- subset_taxa(phyPOLE3S_transform, Phylum %in% names(top10phy.names.POLE3S))
#Saving names and proportions as a data frame then saving as csv
topphylaPOLE3S <- as.data.frame(top10phy.names.POLE3S)
colnames(topphylaPOLE3S)[1] ="Abundance"
write.csv(topphylaPOLE3S, "Top10Phyla_POLE3S_Y1.csv")

## POLESOUT (Firmicutes contam. removed POLESOUT_3_20)
asvdat <- PO
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPO <- phyloseq(ASV,TAX,META)
phyPO_transform <- transform(phyPO, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PO <- sort(tapply(taxa_sums(phyPO_transform), tax_table(phyPO_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPO <- subset_taxa(phyPO_transform, Phylum %in% names(top10phy.names.PO))
#Saving names and proportions as a data frame then saving as csv
topphylaPO <- as.data.frame(top10phy.names.PO)
colnames(topphylaPO)[1] ="Abundance"
write.csv(topphylaPO, "Top10Phyla_PO_Y1.csv")

## RITTAE2 (Firmicutes contam. removed RITTAE2_3_20)
asvdat <- RIT
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyRIT <- phyloseq(ASV,TAX,META)
phyRIT_transform <- transform(phyRIT, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.RIT <- sort(tapply(taxa_sums(phyRIT_transform), tax_table(phyRIT_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaRIT <- subset_taxa(phyRIT_transform, Phylum %in% names(top10phy.names.RIT))
#Saving names and proportions as a data frame then saving as csv
topphylaRIT <- as.data.frame(top10phy.names.RIT)
colnames(topphylaRIT)[1] ="Abundance"
write.csv(topphylaRIT, "Top10Phyla_RIT_Y1.csv")

## S308 
asvdat <- S308
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS308 <- phyloseq(ASV,TAX,META)
phyS308_transform <- transform(phyS308, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S308 <- sort(tapply(taxa_sums(phyS308_transform), tax_table(phyS308_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS308 <- subset_taxa(phyS308_transform, Phylum %in% names(top10phy.names.S308))
#Saving names and proportions as a data frame then saving as csv
topphylaS308 <- as.data.frame(top10phy.names.S308)
colnames(topphylaS308)[1] ="Abundance"
write.csv(topphylaS308, "Top10Phyla_S308_Y1.csv")

## S77  (Firmicutes contam. removed S77_3_20)
asvdat <- S77
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS77 <- phyloseq(ASV,TAX,META)
phyS77_transform <- transform(phyS77, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S77 <- sort(tapply(taxa_sums(phyS77_transform), tax_table(phyS77_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS77 <- subset_taxa(phyS77_transform, Phylum %in% names(top10phy.names.S77))
#Saving names and proportions as a data frame then saving as csv
topphylaS77 <- as.data.frame(top10phy.names.S77)
colnames(topphylaS77)[1] ="Abundance"
write.csv(topphylaS77, "Top10Phyla_S77_Y1.csv")

## S79 (Firmicutes contam. removed S79_3_20)
asvdat <- S79
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS79 <- phyloseq(ASV,TAX,META)
phyS79_transform <- transform(phyS79, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S79 <- sort(tapply(taxa_sums(phyS79_transform), tax_table(phyS79_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS79 <- subset_taxa(phyS79_transform, Phylum %in% names(top10phy.names.S79))
#Saving names and proportions as a data frame then saving as csv
topphylaS79 <- as.data.frame(top10phy.names.S79)
colnames(topphylaS79)[1] ="Abundance"
write.csv(topphylaS79, "Top10Phyla_S79_Y1.csv")

#### Year 2
## CLV10A
asvdat <- CLV
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyCLV <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyCLV_transform <- transform(phyCLV, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.CLV <- sort(tapply(taxa_sums(phyCLV_transform), tax_table(phyCLV_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaCLV <- subset_taxa(phyCLV_transform, Phylum %in% names(top10phy.names.CLV))
#Saving names and proportions as a data frame then saving as csv
topphylaCLV <- as.data.frame(top10phy.names.CLV)
colnames(topphylaCLV)[1] ="Abundance"
write.csv(topphylaCLV, "Top10Phyla_CLV_Y2.csv")

## KISSR0.0 - (Firmicutes removed-> KISSR0.0_3_20)
asvdat <- KISS
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyKISS <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyKISS_transform <- transform(phyKISS, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.KISS <- sort(tapply(taxa_sums(phyKISS_transform), tax_table(phyKISS_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaKISS <- subset_taxa(phyKISS_transform, Phylum %in% names(top10phy.names.KISS))
#Saving names and proportions as a data frame then saving as csv
topphylaKISS <- as.data.frame(top10phy.names.KISS)
colnames(topphylaKISS)[1] ="Abundance"
write.csv(topphylaKISS, "Top10Phyla_KISS_Y2.csv")

## L001 
asvdat <- L1
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL1 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL1_transform <- transform(phyL1, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L1 <- sort(tapply(taxa_sums(phyL1_transform), tax_table(phyL1_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL1 <- subset_taxa(phyL1_transform, Phylum %in% names(top10phy.names.L1))
#Saving names and proportions as a data frame then saving as csv
topphylaL1 <- as.data.frame(top10phy.names.L1)
colnames(topphylaL1)[1] ="Abundance"
write.csv(topphylaL1, "Top10Phyla_L001_Y2.csv")

## L004 
asvdat <- L4
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL4 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL4_transform <- transform(phyL4, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L4 <- sort(tapply(taxa_sums(phyL4_transform), tax_table(phyL4_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL4 <- subset_taxa(phyL4_transform, Phylum %in% names(top10phy.names.L4))
#Saving names and proportions as a data frame then saving as csv
topphylaL4 <- as.data.frame(top10phy.names.L4)
colnames(topphylaL4)[1] ="Abundance"
write.csv(topphylaL4, "Top10Phyla_L004_Y2.csv")

## L005 (Firmicutes removed-> L005_3_20)
asvdat <- L5
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL5 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL5_transform <- transform(phyL5, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L5 <- sort(tapply(taxa_sums(phyL5_transform), tax_table(phyL5_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL5 <- subset_taxa(phyL5_transform, Phylum %in% names(top10phy.names.L5))
#Saving names and proportions as a data frame then saving as csv
topphylaL5 <- as.data.frame(top10phy.names.L5)
colnames(topphylaL5)[1] ="Abundance"
write.csv(topphylaL5, "Top10Phyla_L005_Y2.csv")

## L006
asvdat <- L6
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL6 <- phyloseq(ASV,TAX,META)
phyL6_transform <- transform(phyL6, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L6 <- sort(tapply(taxa_sums(phyL6_transform), tax_table(phyL6_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL6 <- subset_taxa(phyL6_transform, Phylum %in% names(top10phy.names.L6))
#Saving names and proportions as a data frame then saving as csv
topphylaL6 <- as.data.frame(top10phy.names.L6)
colnames(topphylaL6)[1] ="Abundance"
write.csv(topphylaL6, "Top10Phyla_L006_Y2.csv")

## L007
asvdat <- L7
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL7 <- phyloseq(ASV,TAX,META)
phyL7_transform <- transform(phyL7, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L7 <- sort(tapply(taxa_sums(phyL7_transform), tax_table(phyL7_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL7 <- subset_taxa(phyL7_transform, Phylum %in% names(top10phy.names.L7))
#Saving names and proportions as a data frame then saving as csv
topphylaL7 <- as.data.frame(top10phy.names.L7)
colnames(topphylaL7)[1] ="Abundance"
write.csv(topphylaL7, "Top10Phyla_L007_Y2.csv")

## L008
asvdat <- L8
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL8 <- phyloseq(ASV,TAX,META)
phyL8_transform <- transform(phyL8, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L8 <- sort(tapply(taxa_sums(phyL8_transform), tax_table(phyL8_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL8 <- subset_taxa(phyL8_transform, Phylum %in% names(top10phy.names.L8))
#Saving names and proportions as a data frame then saving as csv
topphylaL8 <- as.data.frame(top10phy.names.L8)
colnames(topphylaL8)[1] ="Abundance"
write.csv(topphylaL8, "Top10Phyla_L008_Y2.csv")

## LZ25A 
asvdat <- Z25A
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy25A <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phy25A_transform <- transform(phy25A, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.25A <- sort(tapply(taxa_sums(phy25A_transform), tax_table(phy25A_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla25A <- subset_taxa(phy25A_transform, Phylum %in% names(top10phy.names.25A))
#Saving names and proportions as a data frame then saving as csv
topphyla25A <- as.data.frame(top10phy.names.25A)
colnames(topphyla25A)[1] ="Abundance"
write.csv(topphyla25A, "Top10Phyla_LZ25A_Y2.csv")

## LZ2 (Firmicutes contam. removed LZ2_3_20)
asvdat <- LZ2
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyLZ2 <- phyloseq(ASV,TAX,META)
phyLZ2_transform <- transform(phyLZ2, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.LZ2 <- sort(tapply(taxa_sums(phyLZ2_transform), tax_table(phyLZ2_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaLZ2 <- subset_taxa(phyLZ2_transform, Phylum %in% names(top10phy.names.LZ2))
#Saving names and proportions as a data frame then saving as csv
topphylaLZ2 <- as.data.frame(top10phy.names.LZ2)
colnames(topphylaLZ2)[1] ="Abundance"
write.csv(topphylaLZ2, "Top10Phyla_LZ2_Y2.csv")

## LZ30
asvdat <- Z30
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy30 <- phyloseq(ASV,TAX,META)
phy30_transform <- transform(phy30, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.30 <- sort(tapply(taxa_sums(phy30_transform), tax_table(phy30_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla30 <- subset_taxa(phy30_transform, Phylum %in% names(top10phy.names.30))
#Saving names and proportions as a data frame then saving as csv
topphyla30 <- as.data.frame(top10phy.names.30)
colnames(topphyla30)[1] ="Abundance"
write.csv(topphyla30, "Top10Phyla_LZ30_Y2.csv")

## LZ40 
asvdat <- Z40
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy40 <- phyloseq(ASV,TAX,META)
phy40_transform <- transform(phy40, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.40 <- sort(tapply(taxa_sums(phy40_transform), tax_table(phy40_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla40 <- subset_taxa(phy40_transform, Phylum %in% names(top10phy.names.40))
#Saving names and proportions as a data frame then saving as csv
topphyla40 <- as.data.frame(top10phy.names.40)
colnames(topphyla40)[1] ="Abundance"
write.csv(topphyla40, "Top10Phyla_LZ40_Y2.csv")

## PALMOUT (Firmicutes contam. removed PALMOUT_3_20)
asvdat <- PALM
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPALM <- phyloseq(ASV,TAX,META)
phyPALM_transform <- transform(phyPALM, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PALM <- sort(tapply(taxa_sums(phyPALM_transform), tax_table(phyPALM_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPALM <- subset_taxa(phyPALM_transform, Phylum %in% names(top10phy.names.PALM))
#Saving names and proportions as a data frame then saving as csv
topphylaPALM <- as.data.frame(top10phy.names.PALM)
colnames(topphylaPALM)[1] ="Abundance"
write.csv(topphylaPALM, "Top10Phyla_PALM_Y2.csv")

## PELBAY3 - DONE ON 11/12/22
asvdat <- PEL
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPEL <- phyloseq(ASV,TAX,META)
phyPEL_transform <- transform(phyPEL, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PEL <- sort(tapply(taxa_sums(phyPEL_transform), tax_table(phyPEL_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPEL <- subset_taxa(phyPEL_transform, Phylum %in% names(top10phy.names.PEL))
#Saving names and proportions as a data frame then saving as csv
topphylaPEL <- as.data.frame(top10phy.names.PEL)
colnames(topphylaPEL)[1] ="Abundance"
write.csv(topphylaPEL, "Top10Phyla_PEL_Y2.csv")

## POLE3S (Firmicutes contam. removed POLE3S_3_20)
asvdat <- POLE3S
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPOLE3S <- phyloseq(ASV,TAX,META)
phyPOLE3S_transform <- transform(phyPOLE3S, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.POLE3S <- sort(tapply(taxa_sums(phyPOLE3S_transform), tax_table(phyPOLE3S_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPOLE3S <- subset_taxa(phyPOLE3S_transform, Phylum %in% names(top10phy.names.POLE3S))
#Saving names and proportions as a data frame then saving as csv
topphylaPOLE3S <- as.data.frame(top10phy.names.POLE3S)
colnames(topphylaPOLE3S)[1] ="Abundance"
write.csv(topphylaPOLE3S, "Top10Phyla_POLE3S_Y2.csv")

## POLESOUT (Firmicutes contam. removed POLESOUT_3_20)
asvdat <- PO
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPO <- phyloseq(ASV,TAX,META)
phyPO_transform <- transform(phyPO, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PO <- sort(tapply(taxa_sums(phyPO_transform), tax_table(phyPO_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPO <- subset_taxa(phyPO_transform, Phylum %in% names(top10phy.names.PO))
#Saving names and proportions as a data frame then saving as csv
topphylaPO <- as.data.frame(top10phy.names.PO)
colnames(topphylaPO)[1] ="Abundance"
write.csv(topphylaPO, "Top10Phyla_PO_Y2.csv")

## RITTAE2 (Firmicutes contam. removed RITTAE2_3_20)
asvdat <- RIT
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyRIT <- phyloseq(ASV,TAX,META)
phyRIT_transform <- transform(phyRIT, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.RIT <- sort(tapply(taxa_sums(phyRIT_transform), tax_table(phyRIT_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaRIT <- subset_taxa(phyRIT_transform, Phylum %in% names(top10phy.names.RIT))
#Saving names and proportions as a data frame then saving as csv
topphylaRIT <- as.data.frame(top10phy.names.RIT)
colnames(topphylaRIT)[1] ="Abundance"
write.csv(topphylaRIT, "Top10Phyla_RIT_Y2.csv")

## S308 
asvdat <- S308
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS308 <- phyloseq(ASV,TAX,META)
phyS308_transform <- transform(phyS308, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S308 <- sort(tapply(taxa_sums(phyS308_transform), tax_table(phyS308_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS308 <- subset_taxa(phyS308_transform, Phylum %in% names(top10phy.names.S308))
#Saving names and proportions as a data frame then saving as csv
topphylaS308 <- as.data.frame(top10phy.names.S308)
colnames(topphylaS308)[1] ="Abundance"
write.csv(topphylaS308, "Top10Phyla_S308_Y2.csv")

## S77  (Firmicutes contam. removed S77_3_20)
asvdat <- S77
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS77 <- phyloseq(ASV,TAX,META)
phyS77_transform <- transform(phyS77, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S77 <- sort(tapply(taxa_sums(phyS77_transform), tax_table(phyS77_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS77 <- subset_taxa(phyS77_transform, Phylum %in% names(top10phy.names.S77))
#Saving names and proportions as a data frame then saving as csv
topphylaS77 <- as.data.frame(top10phy.names.S77)
colnames(topphylaS77)[1] ="Abundance"
write.csv(topphylaS77, "Top10Phyla_S77_Y2.csv")

## S79 (Firmicutes contam. removed S79_3_20)
asvdat <- S79
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS79 <- phyloseq(ASV,TAX,META)
phyS79_transform <- transform(phyS79, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S79 <- sort(tapply(taxa_sums(phyS79_transform), tax_table(phyS79_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS79 <- subset_taxa(phyS79_transform, Phylum %in% names(top10phy.names.S79))
#Saving names and proportions as a data frame then saving as csv
topphylaS79 <- as.data.frame(top10phy.names.S79)
colnames(topphylaS79)[1] ="Abundance"
write.csv(topphylaS79, "Top10Phyla_S79_Y2.csv")

#### Year 3
## CLV10A
asvdat <- CLV
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyCLV <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyCLV_transform <- transform(phyCLV, "compositional")
### Assigning Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.CLV <- sort(tapply(taxa_sums(phyCLV_transform), tax_table(phyCLV_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaCLV <- subset_taxa(phyCLV_transform, Phylum %in% names(top10phy.names.CLV))
#Saving names and proportions as a data frame then saving as csv
topphylaCLV <- as.data.frame(top10phy.names.CLV)
colnames(topphylaCLV)[1] ="Abundance"
write.csv(topphylaCLV, "Top10Phyla_CLV_Y3.csv")

## KISSR0.0 - (Firmicutes removed-> KISSR0.0_3_20)
asvdat <- KISS
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyKISS <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyKISS_transform <- transform(phyKISS, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.KISS <- sort(tapply(taxa_sums(phyKISS_transform), tax_table(phyKISS_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaKISS <- subset_taxa(phyKISS_transform, Phylum %in% names(top10phy.names.KISS))
#Saving names and proportions as a data frame then saving as csv
topphylaKISS <- as.data.frame(top10phy.names.KISS)
colnames(topphylaKISS)[1] ="Abundance"
write.csv(topphylaKISS, "Top10Phyla_KISS_Y3.csv")

## L001 
asvdat <- L1
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL1 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL1_transform <- transform(phyL1, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L1 <- sort(tapply(taxa_sums(phyL1_transform), tax_table(phyL1_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL1 <- subset_taxa(phyL1_transform, Phylum %in% names(top10phy.names.L1))
#Saving names and proportions as a data frame then saving as csv
topphylaL1 <- as.data.frame(top10phy.names.L1)
colnames(topphylaL1)[1] ="Abundance"
write.csv(topphylaL1, "Top10Phyla_L001_Y3.csv")

## L004 
asvdat <- L4
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL4 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL4_transform <- transform(phyL4, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L4 <- sort(tapply(taxa_sums(phyL4_transform), tax_table(phyL4_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL4 <- subset_taxa(phyL4_transform, Phylum %in% names(top10phy.names.L4))
#Saving names and proportions as a data frame then saving as csv
topphylaL4 <- as.data.frame(top10phy.names.L4)
colnames(topphylaL4)[1] ="Abundance"
write.csv(topphylaL4, "Top10Phyla_L004_Y3.csv")

## L005 (Firmicutes removed-> L005_3_20)
asvdat <- L5
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL5 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phyL5_transform <- transform(phyL5, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L5 <- sort(tapply(taxa_sums(phyL5_transform), tax_table(phyL5_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL5 <- subset_taxa(phyL5_transform, Phylum %in% names(top10phy.names.L5))
#Saving names and proportions as a data frame then saving as csv
topphylaL5 <- as.data.frame(top10phy.names.L5)
colnames(topphylaL5)[1] ="Abundance"
write.csv(topphylaL5, "Top10Phyla_L005_Y3.csv")

## L006
asvdat <- L6
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL6 <- phyloseq(ASV,TAX,META)
phyL6_transform <- transform(phyL6, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L6 <- sort(tapply(taxa_sums(phyL6_transform), tax_table(phyL6_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL6 <- subset_taxa(phyL6_transform, Phylum %in% names(top10phy.names.L6))
#Saving names and proportions as a data frame then saving as csv
topphylaL6 <- as.data.frame(top10phy.names.L6)
colnames(topphylaL6)[1] ="Abundance"
write.csv(topphylaL6, "Top10Phyla_L006_Y3.csv")

## L007
asvdat <- L7
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL7 <- phyloseq(ASV,TAX,META)
phyL7_transform <- transform(phyL7, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L7 <- sort(tapply(taxa_sums(phyL7_transform), tax_table(phyL7_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL7 <- subset_taxa(phyL7_transform, Phylum %in% names(top10phy.names.L7))
#Saving names and proportions as a data frame then saving as csv
topphylaL7 <- as.data.frame(top10phy.names.L7)
colnames(topphylaL7)[1] ="Abundance"
write.csv(topphylaL7, "Top10Phyla_L007_Y3.csv")

## L008
asvdat <- L8
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyL8 <- phyloseq(ASV,TAX,META)
phyL8_transform <- transform(phyL8, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.L8 <- sort(tapply(taxa_sums(phyL8_transform), tax_table(phyL8_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaL8 <- subset_taxa(phyL8_transform, Phylum %in% names(top10phy.names.L8))
#Saving names and proportions as a data frame then saving as csv
topphylaL8 <- as.data.frame(top10phy.names.L8)
colnames(topphylaL8)[1] ="Abundance"
write.csv(topphylaL8, "Top10Phyla_L008_Y3.csv")

## LZ25A 
asvdat <- Z25A
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy25A <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
phy25A_transform <- transform(phy25A, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.25A <- sort(tapply(taxa_sums(phy25A_transform), tax_table(phy25A_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla25A <- subset_taxa(phy25A_transform, Phylum %in% names(top10phy.names.25A))
#Saving names and proportions as a data frame then saving as csv
topphyla25A <- as.data.frame(top10phy.names.25A)
colnames(topphyla25A)[1] ="Abundance"
write.csv(topphyla25A, "Top10Phyla_LZ25A_Y3.csv")

## LZ2 (Firmicutes contam. removed LZ2_3_20)
asvdat <- LZ2
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyLZ2 <- phyloseq(ASV,TAX,META)
phyLZ2_transform <- transform(phyLZ2, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.LZ2 <- sort(tapply(taxa_sums(phyLZ2_transform), tax_table(phyLZ2_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaLZ2 <- subset_taxa(phyLZ2_transform, Phylum %in% names(top10phy.names.LZ2))
#Saving names and proportions as a data frame then saving as csv
topphylaLZ2 <- as.data.frame(top10phy.names.LZ2)
colnames(topphylaLZ2)[1] ="Abundance"
write.csv(topphylaLZ2, "Top10Phyla_LZ2_Y3.csv")

## LZ30
asvdat <- Z30
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy30 <- phyloseq(ASV,TAX,META)
phy30_transform <- transform(phy30, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.30 <- sort(tapply(taxa_sums(phy30_transform), tax_table(phy30_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla30 <- subset_taxa(phy30_transform, Phylum %in% names(top10phy.names.30))
#Saving names and proportions as a data frame then saving as csv
topphyla30 <- as.data.frame(top10phy.names.30)
colnames(topphyla30)[1] ="Abundance"
write.csv(topphyla30, "Top10Phyla_LZ30_Y3.csv")

## LZ40 
asvdat <- Z40
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phy40 <- phyloseq(ASV,TAX,META)
phy40_transform <- transform(phy40, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.40 <- sort(tapply(taxa_sums(phy40_transform), tax_table(phy40_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phyla40 <- subset_taxa(phy40_transform, Phylum %in% names(top10phy.names.40))
#Saving names and proportions as a data frame then saving as csv
topphyla40 <- as.data.frame(top10phy.names.40)
colnames(topphyla40)[1] ="Abundance"
write.csv(topphyla40, "Top10Phyla_LZ40_Y3.csv")

## PALMOUT (Firmicutes contam. removed PALMOUT_3_20)
asvdat <- PALM
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPALM <- phyloseq(ASV,TAX,META)
phyPALM_transform <- transform(phyPALM, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PALM <- sort(tapply(taxa_sums(phyPALM_transform), tax_table(phyPALM_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPALM <- subset_taxa(phyPALM_transform, Phylum %in% names(top10phy.names.PALM))
#Saving names and proportions as a data frame then saving as csv
topphylaPALM <- as.data.frame(top10phy.names.PALM)
colnames(topphylaPALM)[1] ="Abundance"
write.csv(topphylaPALM, "Top10Phyla_PALM_Y3.csv")

## PELBAY3 - DONE ON 11/12/22
asvdat <- PEL
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPEL <- phyloseq(ASV,TAX,META)
phyPEL_transform <- transform(phyPEL, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PEL <- sort(tapply(taxa_sums(phyPEL_transform), tax_table(phyPEL_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPEL <- subset_taxa(phyPEL_transform, Phylum %in% names(top10phy.names.PEL))
#Saving names and proportions as a data frame then saving as csv
topphylaPEL <- as.data.frame(top10phy.names.PEL)
colnames(topphylaPEL)[1] ="Abundance"
write.csv(topphylaPEL, "Top10Phyla_PEL_Y3.csv")

## POLE3S (Firmicutes contam. removed POLE3S_3_20)
asvdat <- POLE3S
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPOLE3S <- phyloseq(ASV,TAX,META)
phyPOLE3S_transform <- transform(phyPOLE3S, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.POLE3S <- sort(tapply(taxa_sums(phyPOLE3S_transform), tax_table(phyPOLE3S_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPOLE3S <- subset_taxa(phyPOLE3S_transform, Phylum %in% names(top10phy.names.POLE3S))
#Saving names and proportions as a data frame then saving as csv
topphylaPOLE3S <- as.data.frame(top10phy.names.POLE3S)
colnames(topphylaPOLE3S)[1] ="Abundance"
write.csv(topphylaPOLE3S, "Top10Phyla_POLE3S_Y3.csv")

## POLESOUT (Firmicutes contam. removed POLESOUT_3_20)
asvdat <- PO
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyPO <- phyloseq(ASV,TAX,META)
phyPO_transform <- transform(phyPO, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.PO <- sort(tapply(taxa_sums(phyPO_transform), tax_table(phyPO_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaPO <- subset_taxa(phyPO_transform, Phylum %in% names(top10phy.names.PO))
#Saving names and proportions as a data frame then saving as csv
topphylaPO <- as.data.frame(top10phy.names.PO)
colnames(topphylaPO)[1] ="Abundance"
write.csv(topphylaPO, "Top10Phyla_PO_Y3.csv")

## RITTAE2 (Firmicutes contam. removed RITTAE2_3_20)
asvdat <- RIT
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyRIT <- phyloseq(ASV,TAX,META)
phyRIT_transform <- transform(phyRIT, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.RIT <- sort(tapply(taxa_sums(phyRIT_transform), tax_table(phyRIT_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaRIT <- subset_taxa(phyRIT_transform, Phylum %in% names(top10phy.names.RIT))
#Saving names and proportions as a data frame then saving as csv
topphylaRIT <- as.data.frame(top10phy.names.RIT)
colnames(topphylaRIT)[1] ="Abundance"
write.csv(topphylaRIT, "Top10Phyla_RIT_Y3.csv")

## S308 
asvdat <- S308
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS308 <- phyloseq(ASV,TAX,META)
phyS308_transform <- transform(phyS308, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S308 <- sort(tapply(taxa_sums(phyS308_transform), tax_table(phyS308_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS308 <- subset_taxa(phyS308_transform, Phylum %in% names(top10phy.names.S308))
#Saving names and proportions as a data frame then saving as csv
topphylaS308 <- as.data.frame(top10phy.names.S308)
colnames(topphylaS308)[1] ="Abundance"
write.csv(topphylaS308, "Top10Phyla_S308_Y3.csv")

## S77  (Firmicutes contam. removed S77_3_20)
asvdat <- S77
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS77 <- phyloseq(ASV,TAX,META)
phyS77_transform <- transform(phyS77, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S77 <- sort(tapply(taxa_sums(phyS77_transform), tax_table(phyS77_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS77 <- subset_taxa(phyS77_transform, Phylum %in% names(top10phy.names.S77))
#Saving names and proportions as a data frame then saving as csv
topphylaS77 <- as.data.frame(top10phy.names.S77)
colnames(topphylaS77)[1] ="Abundance"
write.csv(topphylaS77, "Top10Phyla_S77_Y3.csv")

## S79 (Firmicutes contam. removed S79_3_20)
asvdat <- S79
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
phyS79 <- phyloseq(ASV,TAX,META)
phyS79_transform <- transform(phyS79, "compositional")
## Top 10 Phyla
#Sort Phylum by abundance and pick the top 10
top10phy.names.S79 <- sort(tapply(taxa_sums(phyS79_transform), tax_table(phyS79_transform)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq data to only the top 10 Phyla
top10phylaS79 <- subset_taxa(phyS79_transform, Phylum %in% names(top10phy.names.S79))
#Saving names and proportions as a data frame then saving as csv
topphylaS79 <- as.data.frame(top10phy.names.S79)
colnames(topphylaS79)[1] ="Abundance"
write.csv(topphylaS79, "Top10Phyla_S79_Y3.csv")

###### Top 10 Phyla by Station - ALL 3 YEARS ######
### You need tidyverse package in order to do this

## Loading in each station on their own (make sure the two columns in the csv is labeled 'Phylum' 'Station Name')
CLV <- read.csv("Top10Phyla_CLV.csv")
KISS <- read.csv("Top10Phyla_KISS.csv")
L1 <- read.csv("Top10Phyla_L001.csv")
L4 <- read.csv("Top10Phyla_L004.csv")
L5 <- read.csv("Top10Phyla_L005.csv")
L6 <- read.csv("Top10Phyla_L006.csv")
L7 <- read.csv("Top10Phyla_L007.csv")
L8 <- read.csv("Top10Phyla_L008.csv")
LZ2 <- read.csv("Top10Phyla_LZ2.csv")
Z25A <- read.csv("Top10Phyla_LZ25A.csv")
Z30 <- read.csv("Top10Phyla_LZ30.csv")
Z40 <- read.csv("Top10Phyla_LZ40.csv")
PALM <- read.csv("Top10Phyla_PALM.csv")
PEL <- read.csv("Top10Phyla_PEL.csv")
POLE3S <- read.csv("Top10Phyla_POLE3S.csv")
PO <- read.csv("Top10Phyla_PO.csv")
RIT <- read.csv("Top10Phyla_RIT.csv")
S308 <- read.csv("Top10Phyla_S308.csv")
S77 <- read.csv("Top10Phyla_S77.csv")
S79 <- read.csv("Top10Phyla_S79.csv")
## Creating a list of the stations
Stations <- list(CLV, KISS, L1, L4, L5, L6, L7, L8, LZ2, Z25A, Z30, Z40, PALM,
                 PEL, POLE3S, PO, RIT, S308, S77, S79)
## Merging all of the data frames in the list (USES TIDYVERSE)
Station_merge <- Stations %>% reduce(full_join, by= "Phylum")
Station_merge[is.na(Station_merge)] = 0  #replacing the NAs with zeros
## Saving merged data frame as CSV
write.csv(Station_merge, "Top10Phyla-Stations.csv")


## Testing to see if I can create a stacked bar chart using the merged station data frame
## Converting the data frame into long format (which converts it into a tibble)
S_tibble <-Station_merge %>% pivot_longer(cols=c(2:21),names_to= "Station",values_to= "Abundance")
write.csv(S_tibble, "StationPhyla_long.csv")
# ## Reloading in previous data frame (went into excel and replaced NA with 0)
# StationPhyla <- read.csv("StationPhyla_long.csv", header = T) or SKIP AND GO TO NEXT LINE!!
StationPhyla <- S_tibble

## Plotting using custom colors
Top10Station <- ggplot(StationPhyla, aes(fill=Phylum, x=Abundance, y=Station)) + 
  geom_bar(position='fill', stat='identity')+      #position="fill" creates a stacked bar plot with abundance as a percentage
  theme_minimal()+
  labs(x='Abundance', y='Stations', title='Top Phyla Found in Lake Okeechobee by Station')+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
Top10Stat <- c("#2bcaf4","#24630e","#000080","#edc427","#1f60aa","#333333","#d841ad",
                        "#41ea27","red3","#806bb4","#cbcc8f","#5f421b","#f08539","#ff9eed")
                        ## listed by phyla in alphabetical order
withr::with_options(list(ggplot2.discrete.fill = Top10Stat),print(Top10Station))

###### Top 10 Phyla by Station - EACH YEAR ######
### You need tidyverse package in order to do this

##Year 1
CLV <- read.csv("Top10Phyla_CLV_Y1.csv")
KISS <- read.csv("Top10Phyla_KISS_Y1.csv")
L1 <- read.csv("Top10Phyla_L001_Y1.csv")
L4 <- read.csv("Top10Phyla_L004_Y1.csv")
L5 <- read.csv("Top10Phyla_L005_Y1.csv")
L6 <- read.csv("Top10Phyla_L006_Y1.csv")
L7 <- read.csv("Top10Phyla_L007_Y1.csv")
L8 <- read.csv("Top10Phyla_L008_Y1.csv")
LZ2 <- read.csv("Top10Phyla_LZ2_Y1.csv")
Z25A <- read.csv("Top10Phyla_LZ25A_Y1.csv")
Z30 <- read.csv("Top10Phyla_LZ30_Y1.csv")
Z40 <- read.csv("Top10Phyla_LZ40_Y1.csv")
PALM <- read.csv("Top10Phyla_PALM_Y1.csv")
PEL <- read.csv("Top10Phyla_PEL_Y1.csv")
POLE3S <- read.csv("Top10Phyla_POLE3S_Y1.csv")
PO <- read.csv("Top10Phyla_PO_Y1.csv")
RIT <- read.csv("Top10Phyla_RIT_Y1.csv")
S308 <- read.csv("Top10Phyla_S308_Y1.csv")
S77 <- read.csv("Top10Phyla_S77_Y1.csv")
S79 <- read.csv("Top10Phyla_S79_Y1.csv")
## Creating a list of the stations
Stations <- list(CLV, KISS, L1, L4, L5, L6, L7, L8, LZ2, Z25A, Z30, Z40, PALM,
                 PEL, POLE3S, PO, RIT, S308, S77, S79)
## Merging all of the data frames in the list (USES TIDYVERSE)
Station_merge <- Stations %>% reduce(full_join, by= "Phylum")
Station_merge[is.na(Station_merge)] = 0  #replacing the NAs with zeros
## Saving merged data frame as CSV
write.csv(Station_merge, "Top10Phyla-Stations_Y1.csv")


## Testing to see if I can create a stacked bar chart using the merged station data frame
## Converting the data frame into long format (which converts it into a tibble)
S_tibble <-Station_merge %>% pivot_longer(cols=c(2:21),names_to= "Station",values_to= "Abundance")
write.csv(S_tibble, "StationPhyla_long_Y1.csv")
# StationPhyla <- read.csv("StationPhyla_long_Y1.csv", header = T) or SKIP AND GO TO NEXT LINE!!
StationPhyla <- S_tibble

## Plotting using custom colors
Top10Station <- ggplot(StationPhyla, aes(fill=Phylum, x=Abundance, y=Station)) + 
  geom_bar(position='fill', stat='identity')+      #position="fill" creates a stacked bar plot with abundance as a percentage
  theme_minimal()+
  labs(x='Abundance', y='Stations', title='Top Phyla Found in Lake Okeechobee by Station - Year 1')+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
Top10Stat <- c("#2bcaf4","#24630e","#edc427","#1f60aa","#333333","#d841ad",
                        "#41ea27","red3","#806bb4","#cbcc8f","#5f421b","#f08539","purple4","#ff9eed")
                        ## listed by phyla in alphabetical order
withr::with_options(list(ggplot2.discrete.fill = Top10Stat),print(Top10Station))

##Year 2
CLV <- read.csv("Top10Phyla_CLV_Y2.csv")
KISS <- read.csv("Top10Phyla_KISS_Y2.csv")
L1 <- read.csv("Top10Phyla_L001_Y2.csv")
L4 <- read.csv("Top10Phyla_L004_Y2.csv")
L5 <- read.csv("Top10Phyla_L005_Y2.csv")
L6 <- read.csv("Top10Phyla_L006_Y2.csv")
L7 <- read.csv("Top10Phyla_L007_Y2.csv")
L8 <- read.csv("Top10Phyla_L008_Y2.csv")
LZ2 <- read.csv("Top10Phyla_LZ2_Y2.csv")
Z25A <- read.csv("Top10Phyla_LZ25A_Y2.csv")
Z30 <- read.csv("Top10Phyla_LZ30_Y2.csv")
Z40 <- read.csv("Top10Phyla_LZ40_Y2.csv")
PALM <- read.csv("Top10Phyla_PALM_Y2.csv")
PEL <- read.csv("Top10Phyla_PEL_Y2.csv")
POLE3S <- read.csv("Top10Phyla_POLE3S_Y2.csv")
PO <- read.csv("Top10Phyla_PO_Y2.csv")
RIT <- read.csv("Top10Phyla_RIT_Y2.csv")
S308 <- read.csv("Top10Phyla_S308_Y2.csv")
S77 <- read.csv("Top10Phyla_S77_Y2.csv")
S79 <- read.csv("Top10Phyla_S79_Y2.csv")
## Creating a list of the stations
Stations <- list(CLV, KISS, L1, L4, L5, L6, L7, L8, LZ2, Z25A, Z30, Z40, PALM,
                 PEL, POLE3S, PO, RIT, S308, S77, S79)
## Merging all of the data frames in the list (USES TIDYVERSE)
Station_merge <- Stations %>% reduce(full_join, by= "Phylum")
Station_merge[is.na(Station_merge)] = 0  #replacing the NAs with zeros
## Saving merged data frame as CSV
write.csv(Station_merge, "Top10Phyla-Stations_Y2.csv")


## Testing to see if I can create a stacked bar chart using the merged station data frame
## Converting the data frame into long format (which converts it into a tibble)
S_tibble <-Station_merge %>% pivot_longer(cols=c(2:21),names_to= "Station",values_to= "Abundance")
write.csv(S_tibble, "StationPhyla_long_Y2.csv")
# StationPhyla <- read.csv("StationPhyla_long_Y2.csv", header = T) or SKIP AND GO TO NEXT LINE!!
StationPhyla <- S_tibble

## Plotting using custom colors
Top10Station <- ggplot(StationPhyla, aes(fill=Phylum, x=Abundance, y=Station)) + 
  geom_bar(position='fill', stat='identity')+      #position="fill" creates a stacked bar plot with abundance as a percentage
  theme_minimal()+
  labs(x='Abundance', y='Stations', title='Top Phyla Found in Lake Okeechobee by Station - Year 2')+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
Top10Stat <- c("#2bcaf4","#24630e","#000080","#edc427","#1f60aa","#333333","#d841ad",
                        "#41ea27","red3","#806bb4","#5f421b","#f08539","#ff9eed")
                        ## listed by phyla in alphabetical order
withr::with_options(list(ggplot2.discrete.fill = Top10Stat),print(Top10Station))

##Year 3
CLV <- read.csv("Top10Phyla_CLV_Y3.csv")
KISS <- read.csv("Top10Phyla_KISS_Y3.csv")
L1 <- read.csv("Top10Phyla_L001_Y3.csv")
L4 <- read.csv("Top10Phyla_L004_Y3.csv")
L5 <- read.csv("Top10Phyla_L005_Y3.csv")
L6 <- read.csv("Top10Phyla_L006_Y3.csv")
L7 <- read.csv("Top10Phyla_L007_Y3.csv")
L8 <- read.csv("Top10Phyla_L008_Y3.csv")
LZ2 <- read.csv("Top10Phyla_LZ2_Y3.csv")
Z25A <- read.csv("Top10Phyla_LZ25A_Y3.csv")
Z30 <- read.csv("Top10Phyla_LZ30_Y3.csv")
Z40 <- read.csv("Top10Phyla_LZ40_Y3.csv")
PALM <- read.csv("Top10Phyla_PALM_Y3.csv")
PEL <- read.csv("Top10Phyla_PEL_Y3.csv")
POLE3S <- read.csv("Top10Phyla_POLE3S_Y3.csv")
PO <- read.csv("Top10Phyla_PO_Y3.csv")
RIT <- read.csv("Top10Phyla_RIT_Y3.csv")
S308 <- read.csv("Top10Phyla_S308_Y3.csv")
S77 <- read.csv("Top10Phyla_S77_Y3.csv")
S79 <- read.csv("Top10Phyla_S79_Y3.csv")
## Creating a list of the stations
Stations <- list(CLV, KISS, L1, L4, L5, L6, L7, L8, LZ2, Z25A, Z30, Z40, PALM,
                 PEL, POLE3S, PO, RIT, S308, S77, S79)
## Merging all of the data frames in the list (USES TIDYVERSE)
Station_merge <- Stations %>% reduce(full_join, by= "Phylum")
Station_merge[is.na(Station_merge)] = 0  #replacing the NAs with zeros
## Saving merged data frame as CSV
write.csv(Station_merge, "Top10Phyla-Stations_Y3.csv")


## Testing to see if I can create a stacked bar chart using the merged station data frame
## Converting the data frame into long format (which converts it into a tibble)
S_tibble <-Station_merge %>% pivot_longer(cols=c(2:21),names_to= "Station",values_to= "Abundance")
write.csv(S_tibble, "StationPhyla_long_Y3.csv")
# StationPhyla <- read.csv("StationPhyla_long_Y3.csv", header = T) or SKIP AND GO TO NEXT LINE!!
StationPhyla <- S_tibble

## Plotting using custom colors
Top10Station <- ggplot(StationPhyla, aes(fill=Phylum, x=Abundance, y=Station)) + 
  geom_bar(position='fill', stat='identity')+      #position="fill" creates a stacked bar plot with abundance as a percentage
  theme_minimal()+
  labs(x='Abundance', y='Stations', title='Top Phyla Found in Lake Okeechobee by Station - Year 3')+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
Top10Stat <- c("#2bcaf4","#24630e","#edc427","#1f60aa","#333333","#d841ad",
                        "#41ea27","red3","#806bb4","#cbcc8f","#a97548","#5f421b","#f08539","#ff9eed")
                        ## listed by phyla in alphabetical order
withr::with_options(list(ggplot2.discrete.fill = Top10Stat),print(Top10Station))


###### Top 15 Orders in Year 3 by Station ######
##Merge feature table with taxonomy and save
dat.Y3 <- read.csv("feature_Y3r_ADJUSTED.csv")
tax <- read.csv("taxonomy_Y123_edited&cleaned.csv")
Yr3t <- merge.data.frame(dat.Y3,tax,by= "FeatureID", all.x = TRUE)
write.csv(Yr3t, "feature_Y3r_ADJUSTED_tax.csv")

##Load feature/tax table
dat.Y3 <- as.data.frame(t(read.csv("feature_Y3r_ADJUSTED_tax.csv", row.names = 1)))

##Separate Station and create master list of top 15
CLV <- as.data.frame(t(dat.Y3[grep("^CLV10A", rownames(dat.Y3)),]))
KISS <- as.data.frame(t(dat.Y3[grep("^KISSR0.0", rownames(dat.Y3)),]))
L1 <- as.data.frame(t(dat.Y3[grep("^L001", rownames(dat.Y3)),]))
L4 <- as.data.frame(t(dat.Y3[grep("^L004", rownames(dat.Y3)),]))
L5 <- as.data.frame(t(dat.Y3[grep("^L005", rownames(dat.Y3)),]))
L6 <- as.data.frame(t(dat.Y3[grep("^L006", rownames(dat.Y3)),]))
L7 <- as.data.frame(t(dat.Y3[grep("^L007", rownames(dat.Y3)),]))
L8 <- as.data.frame(t(dat.Y3[grep("^L008", rownames(dat.Y3)),]))
LZ2 <- as.data.frame(t(dat.Y3[grep("^LZ2_", rownames(dat.Y3)),]))
Z25A <- as.data.frame(t(dat.Y3[grep("^LZ25A", rownames(dat.Y3)),]))
Z30 <- as.data.frame(t(dat.Y3[grep("^LZ30", rownames(dat.Y3)),]))
Z40 <- as.data.frame(t(dat.Y3[grep("^LZ40", rownames(dat.Y3)),]))
PALM <- as.data.frame(t(dat.Y3[grep("^PALMOUT", rownames(dat.Y3)),]))
PEL <- as.data.frame(t(dat.Y3[grep("^PELBAY3", rownames(dat.Y3)),]))
POLE3S <- as.data.frame(t(dat.Y3[grep("^POLE3S", rownames(dat.Y3)),]))
PO <- as.data.frame(t(dat.Y3[grep("^POLESOUT", rownames(dat.Y3)),]))
RIT <- as.data.frame(t(dat.Y3[grep("^RITTAE2", rownames(dat.Y3)),]))
S308 <- as.data.frame(t(dat.Y3[grep("^S308", rownames(dat.Y3)),]))
S77 <- as.data.frame(t(dat.Y3[grep("^S77", rownames(dat.Y3)),]))
S79 <- as.data.frame(t(dat.Y3[grep("^S79", rownames(dat.Y3)),]))

##Assigning top 15 orders by Station
## CLV10A
asvdat <- CLV
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordCLV <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ordCLV_transform <- transform(ordCLV, "compositional")
### Assigning Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.CLV <- sort(tapply(taxa_sums(ordCLV_transform), tax_table(ordCLV_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordCLV <- subset_taxa(ordCLV_transform, Order %in% names(top15ord.names.CLV))
#Saving names and proportions as a data frame then saving as csv
topordCLV <- as.data.frame(top15ord.names.CLV)
colnames(topordCLV)[1] ="Abundance"
write.csv(topordCLV, "Top15Ord_CLV.csv")


## KISSR0.0 - (Firmicutes removed-> KISSR0.0_3_20)
asvdat <- KISS
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordKISS <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ordKISS_transform <- transform(ordKISS, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.KISS <- sort(tapply(taxa_sums(ordKISS_transform), tax_table(ordKISS_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordKISS <- subset_taxa(ordKISS_transform, Order %in% names(top15ord.names.KISS))
#Saving names and proportions as a data frame then saving as csv
topordKISS <- as.data.frame(top15ord.names.KISS)
colnames(topordKISS)[1] ="Abundance"
write.csv(topordKISS, "Top15Ord_KISS.csv")


## L001 
asvdat <- L1
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL1 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ordL1_transform <- transform(ordL1, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L1 <- sort(tapply(taxa_sums(ordL1_transform), tax_table(ordL1_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL1 <- subset_taxa(ordL1_transform, Order %in% names(top15ord.names.L1))
#Saving names and proportions as a data frame then saving as csv
topordL1 <- as.data.frame(top15ord.names.L1)
colnames(topordL1)[1] ="Abundance"
write.csv(topordL1, "Top15Ord_L001.csv")


## L004 
asvdat <- L4
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL4 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ordL4_transform <- transform(ordL4, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L4 <- sort(tapply(taxa_sums(ordL4_transform), tax_table(ordL4_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL4 <- subset_taxa(ordL4_transform, Order %in% names(top15ord.names.L4))
#Saving names and proportions as a data frame then saving as csv
topordL4 <- as.data.frame(top15ord.names.L4)
colnames(topordL4)[1] ="Abundance"
write.csv(topordL4, "Top15Ord_L004.csv")


## L005 (Firmicutes removed-> L005_3_20)
asvdat <- L5
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL5 <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ordL5_transform <- transform(ordL5, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L5 <- sort(tapply(taxa_sums(ordL5_transform), tax_table(ordL5_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL5 <- subset_taxa(ordL5_transform, Order %in% names(top15ord.names.L5))
#Saving names and proportions as a data frame then saving as csv
topordL5 <- as.data.frame(top15ord.names.L5)
colnames(topordL5)[1] ="Abundance"
write.csv(topordL5, "Top15Ord_L005.csv")


## L006
asvdat <- L6
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL6 <- phyloseq(ASV,TAX,META)
ordL6_transform <- transform(ordL6, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L6 <- sort(tapply(taxa_sums(ordL6_transform), tax_table(ordL6_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL6 <- subset_taxa(ordL6_transform, Order %in% names(top15ord.names.L6))
#Saving names and proportions as a data frame then saving as csv
topordL6 <- as.data.frame(top15ord.names.L6)
colnames(topordL6)[1] ="Abundance"
write.csv(topordL6, "Top15Ord_L006.csv")


## L007
asvdat <- L7
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL7 <- phyloseq(ASV,TAX,META)
ordL7_transform <- transform(ordL7, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L7 <- sort(tapply(taxa_sums(ordL7_transform), tax_table(ordL7_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL7 <- subset_taxa(ordL7_transform, Order %in% names(top15ord.names.L7))
#Saving names and proportions as a data frame then saving as csv
topordL7 <- as.data.frame(top15ord.names.L7)
colnames(topordL7)[1] ="Abundance"
write.csv(topordL7, "Top15Ord_L007.csv")

## L008
asvdat <- L8
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordL8 <- phyloseq(ASV,TAX,META)
ordL8_transform <- transform(ordL8, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.L8 <- sort(tapply(taxa_sums(ordL8_transform), tax_table(ordL8_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordL8 <- subset_taxa(ordL8_transform, Order %in% names(top15ord.names.L8))
#Saving names and proportions as a data frame then saving as csv
topordL8 <- as.data.frame(top15ord.names.L8)
colnames(topordL8)[1] ="Abundance"
write.csv(topordL8, "Top15Ord_L008.csv")

## LZ25A 
asvdat <- Z25A
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ord25A <- phyloseq(ASV,TAX,META)
transform <- microbiome::transform
ord25A_transform <- transform(ord25A, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.25A <- sort(tapply(taxa_sums(ord25A_transform), tax_table(ord25A_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ord25A <- subset_taxa(ord25A_transform, Order %in% names(top15ord.names.25A))
#Saving names and proportions as a data frame then saving as csv
topord25A <- as.data.frame(top15ord.names.25A)
colnames(topord25A)[1] ="Abundance"
write.csv(topord25A, "Top15Ord_LZ25A.csv")

## LZ2 (Firmicutes contam. removed LZ2_3_20)
asvdat <- LZ2
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordLZ2 <- phyloseq(ASV,TAX,META)
ordLZ2_transform <- transform(ordLZ2, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.LZ2 <- sort(tapply(taxa_sums(ordLZ2_transform), tax_table(ordLZ2_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordLZ2 <- subset_taxa(ordLZ2_transform, Order %in% names(top15ord.names.LZ2))
#Saving names and proportions as a data frame then saving as csv
topordLZ2 <- as.data.frame(top15ord.names.LZ2)
colnames(topordLZ2)[1] ="Abundance"
write.csv(topordLZ2, "Top15Ord_LZ2.csv")

## LZ30
asvdat <- Z30
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ord30 <- phyloseq(ASV,TAX,META)
ord30_transform <- transform(ord30, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.30 <- sort(tapply(taxa_sums(ord30_transform), tax_table(ord30_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ord30 <- subset_taxa(ord30_transform, Order %in% names(top15ord.names.30))
#Saving names and proportions as a data frame then saving as csv
topord30 <- as.data.frame(top15ord.names.30)
colnames(topord30)[1] ="Abundance"
write.csv(topord30, "Top15Ord_LZ30.csv")

## LZ40 
asvdat <- Z40
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ord40 <- phyloseq(ASV,TAX,META)
ord40_transform <- transform(ord40, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.40 <- sort(tapply(taxa_sums(ord40_transform), tax_table(ord40_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ord40 <- subset_taxa(ord40_transform, Order %in% names(top15ord.names.40))
#Saving names and proportions as a data frame then saving as csv
topord40 <- as.data.frame(top15ord.names.40)
colnames(topord40)[1] ="Abundance"
write.csv(topord40, "Top15Ord_LZ40.csv")

## PALMOUT (Firmicutes contam. removed PALMOUT_3_20)
asvdat <- PALM
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordPALM <- phyloseq(ASV,TAX,META)
ordPALM_transform <- transform(ordPALM, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.PALM <- sort(tapply(taxa_sums(ordPALM_transform), tax_table(ordPALM_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordPALM <- subset_taxa(ordPALM_transform, Order %in% names(top15ord.names.PALM))
#Saving names and proportions as a data frame then saving as csv
topordPALM <- as.data.frame(top15ord.names.PALM)
colnames(topordPALM)[1] ="Abundance"
write.csv(topordPALM, "Top15Ord_PALM.csv")

## PELBAY3 - DONE ON 11/12/22
asvdat <- PEL
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordPEL <- phyloseq(ASV,TAX,META)
ordPEL_transform <- transform(ordPEL, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.PEL <- sort(tapply(taxa_sums(ordPEL_transform), tax_table(ordPEL_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordPEL <- subset_taxa(ordPEL_transform, Order %in% names(top15ord.names.PEL))
#Saving names and proportions as a data frame then saving as csv
topordPEL <- as.data.frame(top15ord.names.PEL)
colnames(topordPEL)[1] ="Abundance"
write.csv(topordPEL, "Top15Ord_PEL.csv")

## POLE3S - DONE ON 11/12/22 (Firmicutes contam. removed POLE3S_3_20)
asvdat <- POLE3S
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordPOLE3S <- phyloseq(ASV,TAX,META)
ordPOLE3S_transform <- transform(ordPOLE3S, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.POLE3S <- sort(tapply(taxa_sums(ordPOLE3S_transform), tax_table(ordPOLE3S_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordPOLE3S <- subset_taxa(ordPOLE3S_transform, Order %in% names(top15ord.names.POLE3S))
#Saving names and proportions as a data frame then saving as csv
topordPOLE3S <- as.data.frame(top15ord.names.POLE3S)
colnames(topordPOLE3S)[1] ="Abundance"
write.csv(topordPOLE3S, "Top15Ord_POLE3S.csv")

## POLESOUT - DONE ON 11/12/22 (Firmicutes contam. removed POLESOUT_3_20)
asvdat <- PO
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordPO <- phyloseq(ASV,TAX,META)
ordPO_transform <- transform(ordPO, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.PO <- sort(tapply(taxa_sums(ordPO_transform), tax_table(ordPO_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordPO <- subset_taxa(ordPO_transform, Order %in% names(top15ord.names.PO))
#Saving names and proportions as a data frame then saving as csv
topordPO <- as.data.frame(top15ord.names.PO)
colnames(topordPO)[1] ="Abundance"
write.csv(topordPO, "Top15Ord_PO.csv")

## RITTAE2 - DONE ON 11/12/22 (Firmicutes contam. removed RITTAE2_3_20)
asvdat <- RIT
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordRIT <- phyloseq(ASV,TAX,META)
ordRIT_transform <- transform(ordRIT, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.RIT <- sort(tapply(taxa_sums(ordRIT_transform), tax_table(ordRIT_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordRIT <- subset_taxa(ordRIT_transform, Order %in% names(top15ord.names.RIT))
#Saving names and proportions as a data frame then saving as csv
topordRIT <- as.data.frame(top15ord.names.RIT)
colnames(topordRIT)[1] ="Abundance"
write.csv(topordRIT, "Top15Ord_RIT.csv")

## S308 
asvdat <- S308
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordS308 <- phyloseq(ASV,TAX,META)
ordS308_transform <- transform(ordS308, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.S308 <- sort(tapply(taxa_sums(ordS308_transform), tax_table(ordS308_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordS308 <- subset_taxa(ordS308_transform, Order %in% names(top15ord.names.S308))
#Saving names and proportions as a data frame then saving as csv
topordS308 <- as.data.frame(top15ord.names.S308)
colnames(topordS308)[1] ="Abundance"
write.csv(topordS308, "Top15Ord_S308.csv")

## S77  (Firmicutes contam. removed S77_3_20)
asvdat <- S77
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordS77 <- phyloseq(ASV,TAX,META)
ordS77_transform <- transform(ordS77, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.S77 <- sort(tapply(taxa_sums(ordS77_transform), tax_table(ordS77_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordS77 <- subset_taxa(ordS77_transform, Order %in% names(top15ord.names.S77))
#Saving names and proportions as a data frame then saving as csv
topordS77 <- as.data.frame(top15ord.names.S77)
colnames(topordS77)[1] ="Abundance"
write.csv(topordS77, "Top15Ord_S77.csv")

## S79 (Firmicutes contam. removed S79_3_20)
asvdat <- S79
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
ordS79 <- phyloseq(ASV,TAX,META)
ordS79_transform <- transform(ordS79, "compositional")
## Top 15 ord
#Sort Order by abundance and pick the top 15
top15ord.names.S79 <- sort(tapply(taxa_sums(ordS79_transform), tax_table(ordS79_transform)[, "Order"], sum), TRUE)[1:15]
#Cut down the phyloseq data to only the top 15 ord
top15ordS79 <- subset_taxa(ordS79_transform, Order %in% names(top15ord.names.S79))
#Saving names and proportions as a data frame then saving as csv
topordS79 <- as.data.frame(top15ord.names.S79)
colnames(topordS79)[1] ="Abundance"
write.csv(topordS79, "Top15Ord_S79.csv")


## Creating a list of the stations
CLV <- read.csv("Top15Ord_CLV.csv")
KISS <- read.csv("Top15Ord_KISS.csv")
L1 <- read.csv("Top15Ord_L001.csv")
L4 <- read.csv("Top15Ord_L004.csv")
L5 <- read.csv("Top15Ord_L005.csv")
L6 <- read.csv("Top15Ord_L006.csv")
L7 <- read.csv("Top15Ord_L007.csv")
L8 <- read.csv("Top15Ord_L008.csv")
LZ2 <- read.csv("Top15Ord_LZ2.csv")
Z25A <- read.csv("Top15Ord_LZ25A.csv")
Z30 <- read.csv("Top15Ord_LZ30.csv")
Z40 <- read.csv("Top15Ord_LZ40.csv")
PALM <- read.csv("Top15Ord_PALM.csv")
PEL <- read.csv("Top15Ord_PEL.csv")
POLE3S <- read.csv("Top15Ord_POLE3S.csv")
PO <- read.csv("Top15Ord_PO.csv")
RIT <- read.csv("Top15Ord_RIT.csv")
S308 <- read.csv("Top15Ord_S308.csv")
S77 <- read.csv("Top15Ord_S77.csv")
S79 <- read.csv("Top15Ord_S79.csv")
## Creating a list of the stations (fix in Excel before moving on!)
Stations <- list(CLV, KISS, L1, L4, L5, L6, L7, L8, LZ2, Z25A, Z30, Z40, PALM,
                 PEL, POLE3S, PO, RIT, S308, S77, S79)
## Merging all of the data frames in the list (USES TIDYVERSE)
Station_merge <- Stations %>% reduce(full_join, by= "Order")
Station_merge[is.na(Station_merge)] = 0  #replacing the NAs with zeros
Station_merge[5,1] <- "NA" #renaming a cell in the dataframe
## Saving merged data frame as CSV
write.csv(Station_merge, "Top15Order-Stations_Y3.csv")


## Testing to see if I can create a stacked bar chart using the merged station data frame
## Converting the data frame into long format (which converts it into a tibble)
S_tibble <-Station_merge %>% pivot_longer(cols=c(2:21),names_to= "Station",values_to= "Abundance")
write.csv(S_tibble, "StationOrders_long_Y3.csv")
# StationOrd <- read.csv("StationPhyla_long_Y3.csv", header = T) or SKIP AND GO TO NEXT LINE!!
StationOrd <- S_tibble

## Plotting using custom colors
Top15Station <- ggplot(StationOrd, aes(fill=Order, x=Abundance, y=Station)) + 
  geom_bar(position='fill', stat='identity')+      #position="fill" creates a stacked bar plot with abundance as a percentage
  theme_minimal()+
  labs(x='Abundance', y='Stations', title='Top Orders Found in Lake Okeechobee by Station - Year 3')+
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5))+
  theme(legend.title = element_text(face="italic"))
Top15Stat <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                               "#78C675","#006ddb","#b66dff","#6db6ff","#b6dbff",
                               "#920000","#924900","#db6d00","navy","#ffff6d", 
                        "antiquewhite2", "#1D91C0", "#67005F", "khaki3", "#CB181D", 
                        "#A6D854", "#F46D43", "#A6CEE3", "#FD8D3C", "#490092", "#999999")
## 15-color palette, colorlblind friendly
withr::with_options(list(ggplot2.discrete.fill = Top15Stat),print(Top15Station))


###### Environmental variable - Scatter plots by Year ######
library(ggplot2)
library(cowplot)
#Loading in metadata
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
#Subsetting metadata table by year
met1 <- metadata[grep("_19$", rownames(metadata)),]
met2 <- metadata[grep("_20$", rownames(metadata)),]
met3 <- metadata[grep("_21$", rownames(metadata)),]
write.csv(met1, "Metadata_BATCH_Y1.csv")
write.csv(met2, "Metadata_BATCH_Y2.csv")
write.csv(met3, "Metadata_BATCH_Y3.csv")

### PLOTTING 
#Chlorophyll a
ch1 <- ggplot(met1, aes(x = as.factor(Month), y = Chlorophyll.a)) +
  geom_jitter(size = 2, color = "green4", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Chlorophyll a (ug/L)")+
  ylim(-25, 150)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

ch2 <- ggplot(met2, aes(x = as.factor(Month), y = Chlorophyll.a)) +
  geom_jitter(size = 2, color = "green4", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-25, 150)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

ch3 <- ggplot(met3, aes(x = as.factor(Month), y = Chlorophyll.a)) +
  geom_jitter(size = 2, color = "green4", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-25, 150)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving as png
png(file="Chla_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(ch1, ch2, ch3, ncol = 3, labels = "AUTO")
graphics.off()


#Total Phosphorus 
tp1 <- ggplot(met1, aes(x = as.factor(Month), y = Phosphate.Total)) +
  geom_jitter(size = 2, color = "darkred", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Total Phosphorus (mg/L)")+
  ylim(0, 0.5)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tp2 <- ggplot(met2, aes(x = as.factor(Month), y = Phosphate.Total)) +
  geom_jitter(size = 2, color = "darkred", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 0.5)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tp3 <- ggplot(met3, aes(x = as.factor(Month), y = Phosphate.Total)) +
  geom_jitter(size = 2, color = "darkred", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 0.5)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving png
png(file="TP_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(tp1, tp2, tp3, ncol = 3, labels = "AUTO")
graphics.off()

#Nitrate + Nitrite
tn1 <- ggplot(met1, aes(x = as.factor(Month), y = Nitrate.Nitrite)) +
  geom_jitter(size = 2, color = "dodgerblue2", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Nitrate + Nitrite (mg/L)")+
  ylim(-0.2, 0.6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tn2 <- ggplot(met2, aes(x = as.factor(Month), y = Nitrate.Nitrite)) +
  geom_jitter(size = 2, color = "dodgerblue2", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.2, 0.6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tn3 <- ggplot(met3, aes(x = as.factor(Month), y = Nitrate.Nitrite)) +
  geom_jitter(size = 2, color = "dodgerblue2", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.2, 0.6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving png
png(file="Nit_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(tn1, tn2, tn3, ncol = 3, labels = "AUTO")
graphics.off()

#Ammonia
a1 <- ggplot(met1, aes(x = as.factor(Month), y = Ammonia)) +
  geom_jitter(size = 2, color = "mediumpurple3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Ammonia (mg/L)")+
  ylim(-0.2, 0.8)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

a2 <- ggplot(met2, aes(x = as.factor(Month), y = Ammonia)) +
  geom_jitter(size = 2, color = "mediumpurple3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.2, 0.8)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

a3 <- ggplot(met3, aes(x = as.factor(Month), y = Ammonia)) +
  geom_jitter(size = 2, color = "mediumpurple3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.2, 0.8)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving png
png(file="Ammonia_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(a1, a2, a3, ncol = 3, labels = "AUTO")
graphics.off()

#Temperature
t1 <- ggplot(met1, aes(x = as.factor(Month), y = Temperature)) +
  geom_jitter(size = 2, color = "sienna", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Temperature (°C)")+
  ylim(0, 35)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

t2 <- ggplot(met2, aes(x = as.factor(Month), y = Temperature)) +
  geom_jitter(size = 2, color = "sienna", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 35)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

t3 <- ggplot(met3, aes(x = as.factor(Month), y = Temperature)) +
  geom_jitter(size = 2, color = "sienna", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 35)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph 
png(file="Temp_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(t1, t2, t3, ncol = 3, labels = "AUTO")
graphics.off()

#Microcystin.LR
m1 <- ggplot(met1, aes(x = as.factor(Month), y = Microcystin.LR)) +
  geom_jitter(size = 2, color = "hotpink3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Microcystin (ug/L)")+
  ylim(-10, 55)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

m2 <- ggplot(met2, aes(x = as.factor(Month), y = Microcystin.LR)) +
  geom_jitter(size = 2, color = "hotpink3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-10, 55)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

m3 <- ggplot(met3, aes(x = as.factor(Month), y = Microcystin.LR)) +
  geom_jitter(size = 2, color = "hotpink3", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-10, 55)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving png
png(file="MicrocystinLR_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(m1, m2, m3, ncol = 3, labels = "AUTO")
graphics.off()

#pH
p1 <- ggplot(met1, aes(x = as.factor(Month), y = pH)) +
  geom_jitter(size = 2, color = "darkorange", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "pH")+
  ylim(0, 11)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

p2 <- ggplot(met2, aes(x = as.factor(Month), y = pH)) +
  geom_jitter(size = 2, color = "darkorange", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 11)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

p3 <- ggplot(met3, aes(x = as.factor(Month), y = pH)) +
  geom_jitter(size = 2, color = "darkorange", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 11)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph 
png(file="PH_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")
graphics.off()

#Total Nitrogen
tn4 <- ggplot(met1, aes(x = as.factor(Month), y = Total.Nitrogen)) +
  geom_jitter(size = 2, color = "navy", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Total Nitrogen (mg/L)")+
  ylim(0, 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tn5 <- ggplot(met2, aes(x = as.factor(Month), y = Total.Nitrogen)) +
  geom_jitter(size = 2, color = "navy", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tn6 <- ggplot(met3, aes(x = as.factor(Month), y = Total.Nitrogen)) +
  geom_jitter(size = 2, color = "navy", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 4)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph 
png(file="TotN_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(tn4, tn5, tn6, ncol = 3, labels = "AUTO")
graphics.off()

#TN:TP
np1 <- ggplot(met1, aes(x = as.factor(Month), y = TN.TP.ratio)) +
  geom_jitter(size = 2, color = "lightsalmon2", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "TN : TP")+
  ylim(0, 46)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

np2 <- ggplot(met2, aes(x = as.factor(Month), y = TN.TP.ratio)) +
  geom_jitter(size = 2, color = "lightsalmon2", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 46)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

np3 <- ggplot(met3, aes(x = as.factor(Month), y = TN.TP.ratio)) +
  geom_jitter(size = 2, color = "lightsalmon2", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 46)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph 
png(file="TNTP_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(np1, np2, np3, ncol = 3, labels = "AUTO")
graphics.off()

#Total Depth
d1 <- ggplot(met1, aes(x = as.factor(Month), y = TotalDepth)) +
  geom_jitter(size = 2, color = "cornsilk4", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Total Depth (m)")+
  ylim(0, 6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

d2 <- ggplot(met2, aes(x = as.factor(Month), y = TotalDepth)) +
  geom_jitter(size = 2, color = "cornsilk4", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

d3 <- ggplot(met3, aes(x = as.factor(Month), y = TotalDepth)) +
  geom_jitter(size = 2, color = "cornsilk4", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(0, 6)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving as png
png(file="Depth_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(d1, d2, d3, ncol = 3, labels = "AUTO")
graphics.off()

#Total Phosphate 
tph1 <- ggplot(met1, aes(x = as.factor(Month), y = Phosphate.Ortho)) +
  geom_jitter(size = 2, color = "grey35", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Total Phosphate (mg/L)")+
  ylim(-0.01, 0.25)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tph2 <- ggplot(met2, aes(x = as.factor(Month), y = Phosphate.Ortho)) +
  geom_jitter(size = 2, color = "grey35", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.01, 0.25)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

tph3 <- ggplot(met3, aes(x = as.factor(Month), y = Phosphate.Ortho)) +
  geom_jitter(size = 2, color = "grey35", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.01, 0.25)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph and saving png
png(file="TPhos_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(tph1, tph2, tph3, ncol = 3, labels = "AUTO")
graphics.off()

###### Viewing Microcystis RA over time ######
#Loading in metadata
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
#Subsetting metadata table by year
met1 <- metadata[grep("_19$", rownames(metadata)),]
met2 <- metadata[grep("_20$", rownames(metadata)),]
met3 <- metadata[grep("_21$", rownames(metadata)),]

#Plotting
mc1 <- ggplot(met1, aes(x = as.factor(Month), y = Microcystis.Abundance)) +
  geom_jitter(size = 1.8, color = "darkcyan", width = 0.25)+
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = "Microcystis Relative Abundance")+
  ylim(-0.01, 0.1)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 1 - 2019") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

mc2 <- ggplot(met2, aes(x = as.factor(Month), y = Microcystis.Abundance)) +
  geom_jitter(size = 1.8, color = "darkcyan", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.01, 0.1)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 2 - 2020") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

mc3 <- ggplot(met3, aes(x = as.factor(Month), y = Microcystis.Abundance)) +
  geom_jitter(size = 1.8, color = "darkcyan", width = 0.25) +
  stat_summary(fun=mean, aes(group=1), geom="line",
               colour="black", linewidth= 0.7)+
  theme_grey()+
  labs(x = "Month", y = NULL)+
  ylim(-0.01, 0.1)+
  theme(legend.position="none")+
  theme(axis.title = element_text(size = 15,face = "bold"))+
  theme(axis.text = element_text(size = 14))+
  labs(title = "Year 3 - 2021") + 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

#Viewing all plots in one graph 
png(file="Microcystis_scatter.png", width=1406, height=573, bg="transparent")
plot_grid(mc1, mc2, mc3, ncol = 3, labels = "AUTO")
graphics.off()

###### Alpha Diversity - Measures ###### 
#### alpha diversity: the species richness that occurs within a given area within a region
#### that is smaller than the entire distribution of the species (Moore, 2013)
#### uses the relative abundance data

###Diversity by Sample (MAKE SURE YOU ONLY HAVE vegan INSTALLED!!)
# Species richness:
S <- as.data.frame(specnumber(dat.01per))
colnames(S)[1] ="Species Richness"
## Species richness: the number of species within a region (Moore, 2013)
#No. individuals:
N <- as.data.frame(rowSums(dat.01per))
colnames(N)[1] ="No. of Individuals"
#Shannon-Weiner Diversity:
H <- as.data.frame(diversity(dat.ra), index="shannon")
colnames(H)[1] ="Shannon Diverisity Index"
## Shannon index: a measure of the information content of a community rather than of the particular species
##               that is present (Moore, 2013) [species richness index]
## strongly influenced by species richness and by rare species (so sample size is negligible)
#Pielou's Evenness:
J = H/log(S)
colnames(J)[1] ="Species Evenness"
## Pielou's evenness: an index that measures diversity along with the species richness
## Formula - J = H/log(S) (aka Shannon evenness index)
## evenness = the count of individuals of each species in an area; 0 is no evenness & 1 is complete evenness 
#Simpson's Diversity (1/D) (inverse):
inv.D <- as.data.frame(diversity(dat.ra, index="inv"))
colnames(inv.D)[1] ="inverse Simpson Diversity Index"
## gives the Simpson index the property of increasing as diversity increases (the dominance of
## a few species decreases)

#Combine data together into a single new data frame, export as CSV
diversitybysample <- cbind(S, N, H, J,inv.D)
write.csv(diversitybysample, "AlphaDiversityBATCH.csv")

#merging with metadata table and export as csv (edited OUTSIDE of R in Excel)
diversitybysample <- read.csv("AlphaDiversityBATCH.csv", row.names = 1)
met <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)
adivmet <- cbind(diversitybysample,met)
write.csv(adivmet,"Metadata-Diversity_BATCH.csv")

###### Alpha Diversity Stats. - ALL YEARS TOGETHER ######
# Packages Used
library(vegan)
library(stats)
library(ggplot2)
library(ggfortify)

#### Alpha Diversities analyses
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)

#### Testing Statistical Significance
## Normality - Shapiro Test (only done on NUMERIC data)
## p <= 0.05 = H0 REJECTED -> DATA IS NOT NORMAL 
## p > 0.05 = H0 ACCEPTED -> DATA IS NORMAL
## Attempted to transform twice using log and sqrt

#Alpha Diversity Variables
shapiro.test(metadata$S) #NOT NORMAL
#W = 0.97921, p-value = 5.777e-07
shapiro.test(metadata$N) #NOT NORMAL
#W = 0.91551, p-value < 2.2e-16
shapiro.test(metadata$H) #NOT NORMAL
#W = 0.96059, p-value = 7.456e-11
shapiro.test(metadata$J) #NOT NORMAL
#W = 0.72606, p-value < 2.2e-16
shapiro.test(metadata$inv.D) #NOT NORMAL
#W = 0.9247, p-value = 8.049e-16



## NOT NORMAL -> Transformations also didn't work -> Non-parametric test (KRUSKAL-WALLIS)
library(pgirmess)
library(multcompView)

#### Hypothesis 1 Comparisons (Diversity & Year)
# Kruskal Wallis: Nonparametric Data (not normal)
## Pairwise Wilcox Test - calculate pairwise comparisons between group levels 
##                        with corrections for multiple testing (non-parametric)

kruskal.test(metadata$S ~ metadata$Year)
#Kruskal-Wallis chi-squared = 13.385, df = 2, p-value = 0.00124 (< 0.05; reject null - significant)
pairwise.wilcox.test(metadata$S, metadata$Year, p.adjust.method = "fdr") #Difference between year 1 and 3 & year 2 and 3
kmc <- kruskalmc(metadata$S ~ metadata$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#  1   2   3 
# "a" "a" "b" 

kruskal.test(metadata$N ~ metadata$Year)
#Kruskal-Wallis chi-squared = 19.73, df = 2, p-value = 5.196e-05 (< 0.05; reject null - significant)
pairwise.wilcox.test(metadata$N, metadata$Year, p.adjust.method = "fdr") #Difference between year 1 and 3 & year 2 and 3
kmc <- kruskalmc(metadata$N ~ metadata$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#  1   2   3 
# "a" "a" "b" 


kruskal.test(metadata$H ~ metadata$Year)
#Kruskal-Wallis chi-squared = 8.5305, df = 2, p-value = 0.01405 (< 0.05; reject null - significant)
pairwise.wilcox.test(metadata$H, metadata$Year, p.adjust.method = "fdr") #Difference between year 2 and 3
kmc <- kruskalmc(metadata$H ~ metadata$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#  1   2   3 
# "ab" "a" "b" 

kruskal.test(metadata$J ~ metadata$Year)
#Kruskal-Wallis chi-squared = 16.987, df = 2, p-value = 0.0002048 (< 0.05; reject null - significant)
pairwise.wilcox.test(metadata$J, metadata$Year, p.adjust.method = "fdr") #Difference between year 1 and 2 & 1 and 3
kmc <- kruskalmc(metadata$J ~ metadata$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#  1   2   3 
# "a" "ab" "b" 

kruskal.test(metadata$inv.D ~ metadata$Year)
#Kruskal-Wallis chi-squared = 16.987, df = 2, p-value = 0.0002048 (< 0.05; reinv.Dect null - significant)
pairwise.wilcox.test(metadata$inv.D, metadata$Year, p.adjust.method = "fdr") #Difference between year 1 and 2 & 1 and 3
kmc <- kruskalmc(metadata$inv.D ~ metadata$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#  1   2   3 
# "a" "b" "a" 


## Plotting boxplots of alpha diversity by year
# Creating pdf for the plots to populate
pdf("AlphaDiverisityPlots.pdf")
par(mfrow=c(2,2))
par(mar=c(5,6,2,2)+0.1)
# plot each boxplot on its own page
boxplot(S~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Richness (S)", line=4.25, cex.lab=1.15)
text(y=1500, x=3, labels="b", col="blue", cex=1.2)
text(y=1420, x=2, labels="a", col="red", cex=1.2)        # labeling which groups are significantly different than the other 
text(y=1585, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,4.5,2,2)+0.1)
boxplot(H~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Shannon Diversity Index (H)", line=2.8, cex.lab=1.15)
text(y=4, x=3, labels="b", col="blue", cex=1.2)
text(y=3.4, x=2, labels="a", col="red", cex=1.2)        
text(y=3.6, x=1, labels="ab", col="purple", cex=1.2)

boxplot(J~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="Species Evenness (J)", line=3, cex.lab=1.15)
text(y=0.73, x=3, labels="b", col="blue", cex=1.2)
text(y=0.685, x=2, labels="ab", col="purple", cex=1.2)        
text(y=0.73, x=1, labels="a", col="red", cex=1.2)

par(mar=c(5,6,2,2)+0.1)
boxplot(inv.D~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="inverse Simpson Diversity Index (inv.D)", line=3.6, cex.lab=1.15)
text(y=440, x=3, labels="a", col="red", cex=1.2)
text(y=420, x=2, labels="b", col="blue", cex=1.2)        
text(y=340, x=1, labels="a", col="red", cex=1.2)

boxplot(N~Year, data=metadata, horizontal = F, las=1, ylab = "", xlab = "")
title(xlab="Year", line = 3, cex.lab=1.15)
title(ylab="No. of Individuals (N)", line=4.25, cex.lab=1.15)
text(y=90000, x=3, labels="b", col="blue", cex=1.2)
text(y=130000, x=2, labels="a", col="red", cex=1.2)        
text(y=160000, x=1, labels="a", col="red", cex=1.2)

# stop saving to pdf 
dev.off()

###### Alpha Diversity by Year ######
Y1 <- metadata[grep("_19$", rownames(metadata)),]
Y2 <- metadata[grep("_20$", rownames(metadata)),]
Y3 <- metadata[grep("_21$", rownames(metadata)),]
## Packages
library(pgirmess)
library(multcompView)
library(vegan)

###### Differences by ZONE - Richness, Shannon, inv. Simpson, Evenness ####
# Boxplot colors by zone (4 different zones so 4 different colors)
Zones <- c("palegreen3","wheat2","rosybrown1","violetred2")


#Year 1
kruskal.test(Y1$S ~ Y1$Zone)
#Kruskal-Wallis chi-squared = 12.026, df = 3, p-value = 0.007295
pairwise.wilcox.test(Y1$S, Y1$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$S ~ Y1$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow  Nearshore   Pelagic     S79 
#  "ab"       "a"       "b"      "ab"

kruskal.test(Y1$H ~ Y1$Zone)
#Kruskal-Wallis chi-squared = 11.77, df = 3, p-value = 0.008214
pairwise.wilcox.test(Y1$H, Y1$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$H ~ Y1$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow Nearshore   Pelagic       S79 
# "ab"       "a"       "b"      "ab" 

kruskal.test(Y1$inv.D ~ Y1$Zone)
#Kruskal-Wallis chi-squared = 8.5961, df = 3, p-value = 0.03517
pairwise.wilcox.test(Y1$inv.D, Y1$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$inv.D ~ Y1$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow Nearshore   Pelagic   S79 
# "a"       "a"       "a"      "a"  -> NO DIFFERENCES

kruskal.test(Y1$J ~ Y1$Zone)
#Kruskal-Wallis chi-squared = 13.726, df = 3, p-value = 0.003303
pairwise.wilcox.test(Y1$J, Y1$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$J ~ Y1$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow Nearshore   Pelagic    S79 
# "a"       "b"      "ab"      "ab" 

## Plotting all the Year 1 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Zone, data=Y1, las=1, col= Zones, ylab = "Species Richness")
boxplot(H~Zone, data=Y1, las=1,col= Zones, ylab = "Shannon Diversity Index")
boxplot(inv.D~Zone, data=Y1, las=1,col= Zones, ylab = "inverse Simpson Diversity Index")
boxplot(J~Zone, data=Y1, las=1,col= Zones, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Zone - Year 1", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)


#Year 2 - NO SIGINIFICANT DIFFERENCES!
kruskal.test(Y2$S ~ Y2$Zone)
#Kruskal-Wallis chi-squared = 2.1354, df = 3, p-value = 0.5448
pairwise.wilcox.test(Y2$S, Y2$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$S ~ Y2$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$H ~ Y2$Zone)
#Kruskal-Wallis chi-squared = 0.90469, df = 3, p-value = 0.8243
pairwise.wilcox.test(Y2$H, Y2$Zone, p.adjust.method = "fdr")
kmc <- kruskalmc(Y2$H ~ Y2$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$inv.D ~ Y2$Zone)
#Kruskal-Wallis chi-squared = 2.1509, df = 3, p-value = 0.5417
pairwise.wilcox.test(Y2$inv.D, Y2$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$inv.D ~ Y2$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$J ~ Y2$Zone)
#Kruskal-Wallis chi-squared = 6.2334, df = 3, p-value = 0.1008
pairwise.wilcox.test(Y2$J, Y2$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$J ~ Y2$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

## Plotting all the Year 2 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Zone, data=Y2, las=1, col= Zones, ylab = "Species Richness")
boxplot(H~Zone, data=Y2, las=1,col= Zones, ylab = "Shannon Diversity Index")
boxplot(inv.D~Zone, data=Y2, las=1,col= Zones, ylab = "inverse Simpson Diversity Index")
boxplot(J~Zone, data=Y2, las=1,col= Zones, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Zone - Year 2", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)


#Year 3
kruskal.test(Y3$S ~ Y3$Zone)
#Kruskal-Wallis chi-squared = 18.21, df = 3, p-value = 0.0003981
pairwise.wilcox.test(Y3$S, Y3$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$S ~ Y3$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow Nearshore   Pelagic    S79 
# "a"       "b"       "a"       "b" 

kruskal.test(Y3$H ~ Y3$Zone)
#Kruskal-Wallis chi-squared = 14.781, df = 3, p-value = 0.002014
pairwise.wilcox.test(Y3$H, Y3$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$H ~ Y3$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# Inflow Nearshore   Pelagic    S79 
#  "a"       "b"      "ab"      "ab" 

kruskal.test(Y3$inv.D ~ Y3$Zone)
#Kruskal-Wallis chi-squared = 13.68, df = 3, p-value = 0.003374
pairwise.wilcox.test(Y3$inv.D, Y3$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$inv.D ~ Y3$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#Inflow Nearshore   Pelagic    S79 
# "a"       "b"       "b"      "ab" 

kruskal.test(Y3$J ~ Y3$Zone)
#Kruskal-Wallis chi-squared = 15.472, df = 3, p-value = 0.001454
pairwise.wilcox.test(Y3$J, Y3$Zone, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$J ~ Y3$Zone) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#Inflow Nearshore   Pelagic   S79 
#"a"       "b"       "b"      "ab" 

## Plotting all the Year 3 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Zone, data=Y3, las=1, col= Zones, ylab = "Species Richness")
boxplot(H~Zone, data=Y3, las=1,col= Zones, ylab = "Shannon Diversity Index")
boxplot(inv.D~Zone, data=Y3, las=1,col= Zones, ylab = "inverse Simpson Diversity Index")
boxplot(J~Zone, data=Y3, las=1,col= Zones, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Zone - Year 3", side = 3, line = -2.4, outer = TRUE, cex = 1.4)

###### Differences by SEASON - Richness, Shannon, inv. Simpson, Evenness ####
# Boxplot colors by season (2 seasons so 2 different colors)
Seasons <- c("lemonchiffon2","royalblue1")

#Year 1 - NO SIGNIFICANT DIFFERENCES ALL AROUND!
kruskal.test(Y1$S ~ Y1$Season)
#Kruskal-Wallis chi-squared = 0.10935, df = 1, p-value = 0.7409
pairwise.wilcox.test(Y1$S, Y1$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$S ~ Y1$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y1$H ~ Y1$Season)
#Kruskal-Wallis chi-squared = 0.18617, df = 1, p-value = 0.6661
pairwise.wilcox.test(Y1$H, Y1$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$H ~ Y1$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y1$inv.D ~ Y1$Season)
#Kruskal-Wallis chi-squared = 0.16256, df = 1, p-value = 0.6868
pairwise.wilcox.test(Y1$inv.D, Y1$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$inv.D ~ Y1$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y1$J ~ Y1$Season)
#Kruskal-Wallis chi-squared = 1.5322, df = 1, p-value = 0.2158
pairwise.wilcox.test(Y1$J, Y1$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$J ~ Y1$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

## Plotting all the Year 1 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Season, data=Y1, las=1, col= Seasons, ylab = "Species Richness")
boxplot(H~Season, data=Y1, las=1,col= Seasons, ylab = "Shannon Diversity Index")
boxplot(inv.D~Season, data=Y1, las=1,col= Seasons, ylab = "inverse Simpson Diversity Index")
boxplot(J~Season, data=Y1, las=1,col= Seasons, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Season - Year 1", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)

#Year 2 - Difference found in evenness
kruskal.test(Y2$S ~ Y2$Season)
#Kruskal-Wallis chi-squared = 0.0066879, df = 1, p-value = 0.9348
pairwise.wilcox.test(Y2$S, Y2$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$S ~ Y2$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$H ~ Y2$Season)
#Kruskal-Wallis chi-squared = 0.018269, df = 1, p-value = 0.8925
pairwise.wilcox.test(Y2$H, Y2$Season, p.adjust.method = "fdr")
kmc <- kruskalmc(Y2$H ~ Y2$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$inv.D ~ Y2$Season)
#Kruskal-Wallis chi-squared = 0.17949, df = 1, p-value = 0.6718
pairwise.wilcox.test(Y2$inv.D, Y2$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$inv.D ~ Y2$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$J ~ Y2$Season)
#Kruskal-Wallis chi-squared = 11.159, df = 1, p-value = 0.0008365
pairwise.wilcox.test(Y2$J, Y2$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$J ~ Y2$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# dry  wet
# "a"  "b"

## Plotting all the Year 2 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Season, data=Y2, las=1, col= Seasons, ylab = "Species Richness")
boxplot(H~Season, data=Y2, las=1,col= Seasons, ylab = "Shannon Diversity Index")
boxplot(inv.D~Season, data=Y2, las=1,col= Seasons, ylab = "inverse Simpson Diversity Index")
boxplot(J~Season, data=Y2, las=1,col= Seasons, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Season - Year 2", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)


#Year 3 - Differences found in evenness
kruskal.test(Y3$S ~ Y3$Season)
#Kruskal-Wallis chi-squared = 2.0537, df = 1, p-value = 0.1518
pairwise.wilcox.test(Y3$S, Y3$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$S ~ Y3$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!! 

kruskal.test(Y3$H ~ Y3$Season)
#Kruskal-Wallis chi-squared = 0.075109, df = 1, p-value = 0.784
pairwise.wilcox.test(Y3$H, Y3$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$H ~ Y3$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!! 

kruskal.test(Y3$inv.D ~ Y3$Season)
#Kruskal-Wallis chi-squared = 0.41548, df = 1, p-value = 0.5192
pairwise.wilcox.test(Y3$inv.D, Y3$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$inv.D ~ Y3$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGINIFICANT DIFFERENCES FOUND!! 

kruskal.test(Y3$J ~ Y3$Season)
#Kruskal-Wallis chi-squared = 4.3677, df = 1, p-value = 0.03663
pairwise.wilcox.test(Y3$J, Y3$Season, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$J ~ Y3$Season) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# dry  wet
# "a"  "b"  

## Plotting all the Year 3 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Season, data=Y3, las=1, col= Seasons, ylab = "Species Richness")
boxplot(H~Season, data=Y3, las=1,col= Seasons, ylab = "Shannon Diversity Index")
boxplot(inv.D~Season, data=Y3, las=1,col= Seasons, ylab = "inverse Simpson Diversity Index")
boxplot(J~Season, data=Y3, las=1,col= Seasons, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Season - Year 3", side = 3, line = -2.4, outer = TRUE, cex = 1.4)

###### Differences by STATION - Richness, Shannon, inv. Simpson, Evenness ####
# Boxplot colors by station
## Expanding the color palette using color ramp
library(RColorBrewer)
nb.cols <- 20 #defines the number of colors you want
Stations <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols) #now the color ramp has 20 colors


#Year 1
kruskal.test(Y1$S ~ Y1$Station)
#Kruskal-Wallis chi-squared = 38.321, df = 19, p-value = 0.0054
pairwise.wilcox.test(Y1$S, Y1$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$S ~ Y1$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0    L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"     "ab"      "a"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"      "b" 
# POLESOUT  RITTAE2     S308      S77      S79 
#     "ab"     "ab"     "ab"     "ab"     "ab" 

kruskal.test(Y1$H ~ Y1$Station)
#Kruskal-Wallis chi-squared = 40.886, df = 19, p-value = 0.002499
pairwise.wilcox.test(Y1$H, Y1$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$H ~ Y1$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0    L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"     "ab"      "a"     "ab"     "ab"     "ab"      "b"     "ab"     "ab"      "b"     "ab"     "ab"     "ab"     "ab"      "b" 
# POLESOUT  RITTAE2     S308      S77      S79 
# "ab"      "b"         "ab"     "ab"     "ab" 

kruskal.test(Y1$inv.D ~ Y1$Station)
#Kruskal-Wallis chi-squared = 40.482, df = 19, p-value = 0.002827
pairwise.wilcox.test(Y1$inv.D, Y1$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$inv.D ~ Y1$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0   L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"     "ab"      "a"     "ab"     "ab"     "b"     "b"     "ab"     "ab"      "b"     "ab"     "ab"     "ab"     "ab"      "b" 
# POLESOUT  RITTAE2     S308      S77      S79 
# "ab"     "b"         "ab"     "ab"     "ab" 

kruskal.test(Y1$J ~ Y1$Station)
#Kruskal-Wallis chi-squared = 34.478, df = 19, p-value = 0.01613
pairwise.wilcox.test(Y1$J, Y1$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$J ~ Y1$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

## Plotting all the Year 1 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Station, data=Y1, las=2, col= Stations, ylab = "Species Richness", xlab = "", cex.axis = 0.88)
boxplot(H~Station, data=Y1, las=2,col= Stations, ylab = "Shannon Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(inv.D~Station, data=Y1, las=2,col= Stations, ylab = "inverse Simpson Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(J~Station, data=Y1, las=2,col= Stations, ylab = "Evenness", xlab = "", cex.axis = 0.88)
#Creating main title
mtext("Alpha Diversity by Station - Year 1", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)


#Year 2 - Differences found in evenness
kruskal.test(Y2$S ~ Y2$Station)
#Kruskal-Wallis chi-squared = 7.7969, df = 19, p-value = 0.9886
pairwise.wilcox.test(Y2$S, Y2$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$S ~ Y2$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$H ~ Y2$Station)
#Kruskal-Wallis chi-squared = 12.192, df = 19, p-value = 0.8772
pairwise.wilcox.test(Y2$H, Y2$Station, p.adjust.method = "fdr")
kmc <- kruskalmc(Y2$H ~ Y2$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$inv.D ~ Y2$Station)
#Kruskal-Wallis chi-squared = 21.503, df = 19, p-value = 0.3097
pairwise.wilcox.test(Y2$inv.D, Y2$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$inv.D ~ Y2$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y2$J ~ Y2$Station)
#Kruskal-Wallis chi-squared = 36.956, df = 19, p-value = 0.008036
pairwise.wilcox.test(Y2$J, Y2$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$J ~ Y2$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0     L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30 
# "ab"     "ab"      "a"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 
# LZ40  PALMOUT  PELBAY3   POLE3S POLESOUT  RITTAE2     S308      S77      S79 
# "b"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 

## Plotting all the Year 2 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Station, data=Y2, las=2, col= Stations, ylab = "Species Richness", xlab = "", cex.axis = 0.88)
boxplot(H~Station, data=Y2, las=2,col= Stations, ylab = "Shannon Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(inv.D~Station, data=Y2, las=2,col= Stations, ylab = "inverse Simpson Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(J~Station, data=Y2, las=2,col= Stations, ylab = "Evenness", xlab = "", cex.axis = 0.88)
#Creating main title
mtext("Alpha Diversity by Station - Year 2", side = 3, line = -2.4, outer = TRUE, cex = 1.4)



#Year 3
kruskal.test(Y3$S ~ Y3$Station)
#Kruskal-Wallis chi-squared = 36.513, df = 19, p-value = 0.009123
pairwise.wilcox.test(Y3$S, Y3$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$S ~ Y3$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0    L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"     "ab"      "a"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 
# POLESOUT  RITTAE2     S308      S77      S79 
# "ab"     "ab"         "ab"     "ab"      "b" 

kruskal.test(Y3$H ~ Y3$Station)
#Kruskal-Wallis chi-squared = 37.551, df = 19, p-value = 0.006766
pairwise.wilcox.test(Y3$H, Y3$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$H ~ Y3$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0     L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"     "ab"        "a"     "ab"     "ab"     "ab"      "b"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 
# POLESOUT  RITTAE2     S308      S77      S79 
# "ab"     "ab"         "ab"     "ab"     "ab" 

kruskal.test(Y3$inv.D ~ Y3$Station)
#Kruskal-Wallis chi-squared = 42.098, df = 19, p-value = 0.001719
pairwise.wilcox.test(Y3$inv.D, Y3$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$inv.D ~ Y3$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0     L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30     LZ40  PALMOUT  PELBAY3   POLE3S 
# "ab"      "ab"        "a"     "ab"     "ab"     "ab"      "b"   "ab"      "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 
# POLESOUT  RITTAE2     S308      S77      S79 
# "ab"     "ab"         "ab"     "ab"     "ab" 

kruskal.test(Y3$J ~ Y3$Station)
#Kruskal-Wallis chi-squared = 42.614, df = 19, p-value = 0.001463
pairwise.wilcox.test(Y3$J, Y3$Station, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$J ~ Y3$Station) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# CLV10A KISSR0.0     L001     L004     L005     L006     L007     L008      LZ2    LZ25A     LZ30 
# "ab"     "ab"      "a"     "ab"     "ab"     "ab"      "b"     "ab"     "ab"     "ab"     "ab" 
# LZ40  PALMOUT  PELBAY3   POLE3S POLESOUT  RITTAE2     S308      S77      S79 
# "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab"     "ab" 

## Plotting all the Year 3 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Station, data=Y3, las=2, col= Stations, ylab = "Species Richness", xlab = "", cex.axis = 0.88)
boxplot(H~Station, data=Y3, las=2,col= Stations, ylab = "Shannon Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(inv.D~Station, data=Y3, las=2,col= Stations, ylab = "inverse Simpson Diversity Index", xlab = "", cex.axis = 0.88)
boxplot(J~Station, data=Y3, las=2,col= Stations, ylab = "Evenness", xlab = "", cex.axis = 0.88)
#Creating main title
mtext("Alpha Diversity by Station - Year 3", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)

###### Differences by MONTH - Richness, Shannon, inv. Simpson, Evenness ####
# Boxplot colors by month (different for each year)
Year1col <- c("lightgoldenrod1","goldenrod1","green3","cadetblue2","dodgerblue2",
                               "mediumpurple2","lightpink1","tan","sienna","seashell3")
Year2col <- c("firebrick2","darkorange1","lightgoldenrod1","goldenrod1","green3",
                          "cadetblue2","dodgerblue2","mediumpurple2","lightpink1",
                          "tan","sienna","seashell3")
Year3col <- c("firebrick2","darkorange1","lightgoldenrod1","goldenrod1","green3",
                          "cadetblue2","dodgerblue2","mediumpurple2","lightpink1",
                          "tan")
#Year 1
kruskal.test(Y1$S ~ Y1$Month)
#Kruskal-Wallis chi-squared = 26.535, df = 9, p-value = 0.001669
pairwise.wilcox.test(Y1$S, Y1$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$S ~ Y1$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   3     4     5     6     7     8     9    10    11    12 
# "ab" "abc" "abc"   "a" "abc"   "c" "abc" "abc"  "bc" "abc" 

kruskal.test(Y1$H ~ Y1$Month)
#Kruskal-Wallis chi-squared = 25.593, df = 9, p-value = 0.002381
pairwise.wilcox.test(Y1$H, Y1$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$H ~ Y1$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   3    4    5    6    7    8    9   10   11   12 
# "ab" "ab" "ab"  "a" "ab" "ab" "ab" "ab"  "b" "ab" 

kruskal.test(Y1$inv.D ~ Y1$Month)
#Kruskal-Wallis chi-squared = 18.778, df = 9, p-value = 0.02715
pairwise.wilcox.test(Y1$inv.D, Y1$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$inv.D ~ Y1$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

kruskal.test(Y1$J ~ Y1$Month)
#Kruskal-Wallis chi-squared = 13.89, df = 9, p-value = 0.1263
pairwise.wilcox.test(Y1$J, Y1$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y1$J ~ Y1$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#NO SIGNIFICANT DIFFERENCES FOUND!!

## Plotting all the Year 1 boxplots on one graph
#defining plotting area as one row and 4 columns
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Month, data=Y1, las=1, col= Year1col, ylab = "Species Richness")
boxplot(H~Month, data=Y1, las=1,col= Year1col, ylab = "Shannon Diversity Index")
boxplot(inv.D~Month, data=Y1, las=1,col= Year1col, ylab = "inverse Simpson Diversity Index")
boxplot(J~Month, data=Y1, las=1,col= Year1col, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Month - Year 1", side = 3, line = -2.4, outer = TRUE, cex = 1.4)


#Year 2
kruskal.test(Y2$S ~ Y2$Month)
#Kruskal-Wallis chi-squared = 144.03, df = 11, p-value < 2.2e-16
pairwise.wilcox.test(Y2$S, Y2$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$S ~ Y2$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# 1      2      3      4      5      6      7      8      9     10     11     12 
# "abc"  "abc"    "d"   "de"    "d"    "d"  "ade" "abce"   "bc"    "b"    "b" "acde" 

kruskal.test(Y2$H ~ Y2$Month)
#Kruskal-Wallis chi-squared = 131.82, df = 11, p-value < 2.2e-16
pairwise.wilcox.test(Y2$H, Y2$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$H ~ Y2$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# 1     2     3     4     5     6     7     8     9    10    11    12 
# "ab"  "ab"   "c"  "cd"  "cd"  "cd" "acd"  "ab"  "ab"  "ab"   "b"  "ad" 

kruskal.test(Y2$inv.D ~ Y2$Month)
#Kruskal-Wallis chi-squared = 104.87, df = 11, p-value < 2.2e-16
pairwise.wilcox.test(Y2$inv.D, Y2$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$inv.D ~ Y2$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
# 1     2     3     4     5     6     7     8     9    10    11    12 
# "a"   "a"   "b"  "bc"  "bc"  "bc" "abc"   "a"   "a"   "a"   "a"  "ac" 

kruskal.test(Y2$J ~ Y2$Month)
#Kruskal-Wallis chi-squared = 34.984, df = 11, p-value = 0.0002494
pairwise.wilcox.test(Y2$J, Y2$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y2$J ~ Y2$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   1    2    3    4    5    6    7    8    9   10   11   12 
# "ab"  "a" "ab" "ab" "ab"  "b" "ab" "ab" "ab" "ab"  "b" "ab" 

## Plotting all the Year 2 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Month, data=Y2, las=1, col= Year2col, ylab = "Species Richness")
boxplot(H~Month, data=Y2, las=1,col= Year2col, ylab = "Shannon Diversity Index")
boxplot(inv.D~Month, data=Y2, las=1,col= Year2col, ylab = "inverse Simpson Diversity Index")
boxplot(J~Month, data=Y2, las=1,col= Year2col, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Month - Year 2", side = 3, line = - 2.4, outer = TRUE, cex = 1.4)


#Year 3
kruskal.test(Y3$S ~ Y3$Month)
#Kruskal-Wallis chi-squared = 50.462, df = 9, p-value = 8.819e-08
pairwise.wilcox.test(Y3$S, Y3$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$S ~ Y3$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   1     2     3     4     5     6     7     8     9    10 
# "ab"   "c"   "a"  "ab" "abc"  "bc"  "bc"  "bc" "abc"  "ab" 

kruskal.test(Y3$H ~ Y3$Month)
#Kruskal-Wallis chi-squared = 45.298, df = 9, p-value = 8.126e-07
pairwise.wilcox.test(Y3$H, Y3$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$H ~ Y3$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#    1     2     3     4     5     6     7     8     9    10 
# "abc"   "a"   "b" "abc"  "ac"   "a" "abc" "abc" "abc"  "bc"   

kruskal.test(Y3$inv.D ~ Y3$Month)
#Kruskal-Wallis chi-squared = 38.56, df = 9, p-value = 1.383e-05
pairwise.wilcox.test(Y3$inv.D, Y3$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$inv.D ~ Y3$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   1    2    3    4    5    6    7    8    9   10 
# "ab"  "a"  "b" "ab" "ab"  "a" "ab" "ab" "ab"  "b" 

kruskal.test(Y3$J ~ Y3$Month)
#Kruskal-Wallis chi-squared = 36.807, df = 9, p-value = 2.848e-05
pairwise.wilcox.test(Y3$J, Y3$Month, p.adjust.method = "fdr") 
kmc <- kruskalmc(Y3$J ~ Y3$Month) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$difference # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   1    2    3    4    5    6    7    8    9   10 
# "ab"  "a"  "b" "ab" "ab" "ab" "ab"  "b" "ab"  "b" 

## Plotting all the Year 3 boxplots on one graph
par(mfrow = c(1,4))
#plotting the boxplots for each alpha diversity variable
boxplot(S~Month, data=Y3, las=1, col= Year3col, ylab = "Species Richness")
boxplot(H~Month, data=Y3, las=1,col= Year3col, ylab = "Shannon Diversity Index")
boxplot(inv.D~Month, data=Y3, las=1,col= Year3col, ylab = "inverse Simpson Diversity Index")
boxplot(J~Month, data=Y3, las=1,col= Year3col, ylab = "Evenness")
#Creating main title
mtext("Alpha Diversity by Month - Year 3", side = 3, line = -2.4, outer = TRUE, cex = 1.4)



###### Correlation of alpha diversity measures and chlorophyll a ######
metadata <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)
par(mfrow=c(2,2))
#Shannon vs. Chl.a
#calculating correlation (-1 to 0 to +1; negatively correlated to none to positively correlated)
cor.test(metadata$Chlorophyll.a, metadata$H, method ="pearson")
#t = -0.74435, df = 539, p-value = 0.457, Pearson coeff. = -0.03204502 <- NOT SIGNIFICANT
#plotting them against each other
plot(metadata$Chlorophyll.a, metadata$H, pch = 19, col = "gray52", xlab = "", ylab = "")
# Adding text
title(main="Shannon Diversity vs Chlorophyll-a Correlation",
      xlab = "Chlorophyll a (ug/L)",
      ylab = "Shannon Diversity Index")
#inv.Simpson vs. Chl.a 
cor.test(metadata$Chlorophyll.a, metadata$inv.D, method ="pearson")
# t = 1.1217, df = 539, p-value = 0.2625, Pearson coeff. = 0.04825728 <- NOT SIGNIFICANT
plot(metadata$Chlorophyll.a, metadata$inv.D, pch = 19, col = "gray52", xlab = "", ylab = "")
title(main="inverse Simpson Diversity vs Chlorophyll-a Correlation",
      xlab = "Chlorophyll a (ug/L)",
      ylab = "inverse Simpson Diversity Index")
#Richness vs. Chl.a
cor.test(metadata$Chlorophyll.a, metadata$S, method ="pearson")
# t = 0.49649, df = 539, p-value = 0.6198, Pearson coeff. = 0.0213804 <- NOT SIGNIFICANT
plot(metadata$Chlorophyll.a, metadata$S, pch = 19, col = "gray52", xlab = "", ylab = "")
title(main="Species Richness vs Chlorophyll-a Correlation",
      xlab = "Chlorophyll a (ug/L)",
      ylab = "Species Richness")
#Evenness vs. Chl.a
cor.test(metadata$Chlorophyll.a, metadata$J, method ="pearson")
# t = -1.9153, df = 539, p-value = 0.05599, Pearson coeff. = -0.08221648   <- NOT SIGNIFICANT
plot(metadata$Chlorophyll.a, metadata$J, pch = 19, col = "gray52", xlab = "", ylab = "")
title(main="Evenness vs Chlorophyll-a Correlation",
      xlab = "Chlorophyll a (ug/L)",
      ylab = "Evenness")

###### Correlation of Microcystis  vs. Chl a (and Microcystin LR) ######
metadata <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)
## Chl a
#calculating correlation (-1 to 0 to +1; negatively correlated to no correlation to positively correlated)
cor.test(metadata$Chlorophyll.a, metadata$Microcystis.Abundance, method ="pearson")
# t = 5.4696, df = 539, p-value = 6.914e-08, Pearson coeff. = 0.229314    -> weakly positive (SIGNIFICANT)       
#plotting them against each other
plot(metadata$Microcystis.Abundance, metadata$Chlorophyll.a, pch = 19, xlab = "", ylab = "")
lines(lowess(metadata$Microcystis.Abundance, metadata$Chlorophyll.a), col = 2, lwd = 2)
# Adding text
title(main="Microcystis Relative Abundance vs Chlorophyll-a Correlation",
      xlab = "Microcystis Relative Abundance",
      ylab = "Chlorophyll a (ug/L)")
text(0.063,136,"Pearson R: 0.23", cex=1.05)
## Microcystin
cor.test(metadata$Microcystin.LR, metadata$Microcystis.Abundance, method ="pearson")
# t = 17.318, df = 539, p-value < 2.2e-16, Pearson coeff. = 0.5979055     -> positive (SIGNIFICANT)       
#plotting them against each other
plot(metadata$Microcystis.Abundance, metadata$Microcystin.LR, pch = 19, xlab = "", ylab = "")
lines(lowess(metadata$Microcystis.Abundance, metadata$Microcystin.LR), col = 2, lwd = 2)
# Adding text
title(main="Microcystis Relative Abundance vs Microcystin (ug/L) Correlation",
      xlab = "Microcystis Relative Abundance",
      ylab = "Microcystin (ug/L)")
text(0.063,45,"Pearson R: 0.60", cex=1.05)

###### Alpha Diversity vs Microcystis Abundance ######
metadata <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)
##Scatter plots
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Microcystis.Abundance, metadata$S, pch = 19, xlab= "Microcystis Relative Abundance",
     ylab = "Species Richness")
lines(lowess(metadata$Microcystis.Abundance, metadata$S), col = 2, lwd = 2)
title(main="Species richness vs Microcystis Relative Abundance", cex.main = 1)

plot(metadata$Microcystis.Abundance, metadata$H, pch = 19, xlab= "Microcystis Relative Abundance",
     ylab = "Shannon Diversity Index")
lines(lowess(metadata$Microcystis.Abundance, metadata$H), col = 2, lwd = 2)
text(0.062,6,"Pearson's r = -0.23", cex=0.9)
title(main="Shannon Diversity vs Microcystis Relative Abundance", cex.main = 1)

plot(metadata$Microcystis.Abundance, metadata$J, pch = 19, xlab= "Microcystis Relative Abundance",
     ylab = "Evenness")
lines(lowess(metadata$Microcystis.Abundance, metadata$J), col = 2, lwd = 2)
text(0.062,0.88,"Pearson's r = -0.72", cex=0.9)
title(main="Species Evenness vs Microcystis Relative Abundance", cex.main = 1)

plot(metadata$Microcystis.Abundance, metadata$inv.D, pch = 19, xlab= "Microcystis Relative Abundance",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Microcystis.Abundance, metadata$inv.D), col = 2, lwd = 2)
text(0.062,400,"Pearson's r = -0.22", cex=0.9)
title(main="inverse Simpson Diversity vs Microcystis Relative Abundance", cex.main = 1)


## Looking at the correlations
cor.test(metadata$Microcystis.Abundance, metadata$S, method ="pearson")
#t = 1.4678, df = 539, Pearson coeff. = 0.0630954 , p-value = 0.1427 -> NOT SIGNIFICANT (NO CORRELATION)
cor.test(metadata$Microcystis.Abundance, metadata$H, method ="pearson")
#t = -5.5028, df = 539, Pearson coeff. = -0.2306343, p-value = 5.785e-08 -> SIGNIFICANT (NEG. CORRELATION)
cor.test(metadata$Microcystis.Abundance, metadata$J, method ="pearson")
#t = -24.34, df = 539, Pearson coeff. = -0.7236151, p-value < 2.2e-16 -> SIGNIFICANT (NEG. CORRELATION)
cor.test(metadata$Microcystis.Abundance, metadata$inv.D, method ="pearson")
#t = -5.3297, df = 539, Pearson coeff. = -0.2237471, p-value = 1.448e-07 -> SIGNIFICANT (NEG. CORRELATION)

###### Alpha Diversity vs Environmental Variables - Scatter plots ######
metadata <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)
##Scatter plots
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph

#Chlorophyll a
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Chlorophyll.a, metadata$S, pch = 19, xlab= "Chlorophyll a (ug/L)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Chlorophyll a (ug/L)", cex.main = 1)
plot(metadata$Chlorophyll.a, metadata$H, pch = 19, xlab= "Chlorophyll a (ug/L)",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs Chlorophyll a (ug/L)", cex.main = 1)
plot(metadata$Chlorophyll.a, metadata$J, pch = 19, xlab= "Chlorophyll a (ug/L)",
     ylab = "Evenness", col="grey54")
title(main="Species Evenness vs Chlorophyll a (ug/L)", cex.main = 1)
plot(metadata$Chlorophyll.a, metadata$inv.D, pch = 19, xlab= "Chlorophyll a (ug/L)",
     ylab = "inverse Simpson Diversity Index", col="grey54")
title(main="inverse Simpson Diversity vs Chlorophyll a (ug/L)", cex.main = 1)

#Ammonia
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Ammonia, metadata$S, pch = 19, xlab= "Ammonia (mg/L)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Ammonia (mg/L)", cex.main = 1)
plot(metadata$Ammonia, metadata$H, pch = 19, xlab= "Ammonia (mg/L)",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs Ammonia (mg/L)", cex.main = 1)
plot(metadata$Ammonia, metadata$J, pch = 19, xlab= "Ammonia (mg/L)",
     ylab = "Evenness")
lines(lowess(metadata$Ammonia, metadata$J), col = 2, lwd = 2)
text(0.68,0.8,"Pearson's r = 0.11", cex=0.9)
title(main="Species Evenness vs Ammonia (mg/L)", cex.main = 1)
plot(metadata$Ammonia, metadata$inv.D, pch = 19, xlab= "Ammonia (mg/L)",
     ylab = "inverse Simpson Diversity Index", col="grey54")
title(main="inverse Simpson Diversity vs Ammonia (mg/L)", cex.main = 1)

#Nitrate(ite)
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Nitrate.Nitrite, metadata$S, pch = 19, xlab= "Nitrate + Nitrite (mg/L)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Nitrate + Nitrite (mg/L)", cex.main = 1)
plot(metadata$Nitrate.Nitrite, metadata$H, pch = 19, xlab= "Nitrate + Nitrite (mg/L)",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs Nitrate + Nitrite (mg/L)", cex.main = 1)
plot(metadata$Nitrate.Nitrite, metadata$J, pch = 19, xlab= "Nitrate + Nitrite (mg/L)",
     ylab = "Evenness")
lines(lowess(metadata$Nitrate.Nitrite, metadata$J), col = 2, lwd = 2)
text(0.5,0.7,"Pearson's r = -0.10", cex=0.9)
title(main="Species Evenness vs Nitrate + Nitrite (mg/L)", cex.main = 1)
plot(metadata$Nitrate.Nitrite, metadata$inv.D, pch = 19, xlab= "Nitrate + Nitrite (mg/L)",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Nitrate.Nitrite, metadata$inv.D), col = 2, lwd = 2)
text(0.52,400,"Pearson's r = -0.10", cex=0.9)
title(main="inverse Simpson Diversity vs Nitrate + Nitrite (mg/L)", cex.main = 1)

#Total Phosphorus
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Phosphate.Total, metadata$S, pch = 19, xlab= "Total Phosphorus (mg/L)",
     ylab = "Species Richness")
lines(lowess(metadata$Phosphate.Total, metadata$S), col = 2, lwd = 2)
text(0.45,1750,"Pearson's r = 0.18", cex=0.9)
title(main="Species richness vs Total Phosphorus (mg/L)", cex.main = 1)
plot(metadata$Phosphate.Total, metadata$H, pch = 19, xlab= "Total Phosphorus (mg/L)",
     ylab = "Shannon Diversity Index")
lines(lowess(metadata$Phosphate.Total, metadata$H), col = 2, lwd = 2)
text(0.44,6.4,"Pearson's r = 0.06", cex=0.9)
title(main="Shannon Diversity vs Total Phosphorus (mg/L)", cex.main = 1)
plot(metadata$Phosphate.Total, metadata$J, pch = 19, xlab= "Total Phosphorus (mg/L)",
     ylab = "Evenness", col="grey54")
title(main="Species Evenness vs Total Phosphorus (mg/L)", cex.main = 1)
plot(metadata$Phosphate.Total, metadata$inv.D, pch = 19, xlab= "Total Phosphorus (mg/L)",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Phosphate.Total, metadata$inv.D), col = 2, lwd = 2)
text(0.44,400,"Pearson's r = 0.10", cex=0.9)
title(main="inverse Simpson Diversity vs Total Phosphorus (mg/L)", cex.main = 1)

#Microcystin
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Microcystin.LR, metadata$S, pch = 19, xlab= "Microcystin (ug/L)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Microcystin (ug/L)", cex.main = 1)
plot(metadata$Microcystin.LR, metadata$H, pch = 19, xlab= "Microcystin (ug/L)",
     ylab = "Shannon Diversity Index")
lines(lowess(metadata$Microcystin.LR, metadata$H), col = 2, lwd = 2)
text(48,6.2,"Pearson's r = -0.23", cex=0.9)
title(main="Shannon Diversity vs Microcystin (ug/L)", cex.main = 1)
plot(metadata$Microcystin.LR, metadata$J, pch = 19, xlab= "Microcystin (ug/L)",
     ylab = "Evenness")
lines(lowess(metadata$Microcystin.LR, metadata$J), col = 2, lwd = 2)
text(48,0.88,"Pearson's r = -0.49", cex=0.9)
title(main="Species Evenness vs Microcystin (ug/L)", cex.main = 1)
plot(metadata$Microcystin.LR, metadata$inv.D, pch = 19, xlab= "Microcystin (ug/L)",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Microcystin.LR, metadata$inv.D), col = 2, lwd = 2)
text(46,375,"Pearson's r = -0.20", cex=0.9)
title(main="inverse Simpson Diversity vs Microcystin (ug/L)", cex.main = 1)

#Temperature
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Temperature, metadata$S, pch = 19, xlab= "Temperature (°C)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Temperature (°C)", cex.main = 1)
plot(metadata$Temperature, metadata$H, pch = 19, xlab= "Temperature (°C)",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs Temperature (°C)", cex.main = 1)
plot(metadata$Temperature, metadata$J, pch = 19, xlab= "Temperature (°C)",
     ylab = "Evenness", col="grey54")
title(main="Species Evenness vs Temperature (°C)", cex.main = 1)
plot(metadata$Temperature, metadata$inv.D, pch = 19, xlab= "Temperature (°C)",
     ylab = "inverse Simpson Diversity Index", col="grey54")
title(main="inverse Simpson Diversity vs Temperature (°C)", cex.main = 1)

#Total Nitrogen
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Total.Nitrogen, metadata$S, pch = 19, xlab= "Total Nitrogen (mg/L)",
     ylab = "Species Richness")
lines(lowess(metadata$Total.Nitrogen, metadata$S), col = 2, lwd = 2)
text(3.15,1750,"Pearson's r = 0.17", cex=0.9)
title(main="Species richness vs Total Nitrogen (mg/L)", cex.main = 1)
plot(metadata$Total.Nitrogen, metadata$H, pch = 19, xlab= "Total Nitrogen (mg/L)",
     ylab = "Shannon Diversity Index")
lines(lowess(metadata$Total.Nitrogen, metadata$H), col = 2, lwd = 2)
text(3,6.8,"Pearson's r = 0.13", cex=0.9)
title(main="Shannon Diversity vs Total Nitrogen (mg/L)", cex.main = 1)
plot(metadata$Total.Nitrogen, metadata$J, pch = 19, xlab= "Total Nitrogen (mg/L)",
     ylab = "Evenness", col="grey54")
title(main="Species Evenness vs Total Nitrogen (mg/L)", cex.main = 1)
plot(metadata$Total.Nitrogen, metadata$inv.D, pch = 19, xlab= "Total Nitrogen (mg/L)",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Total.Nitrogen, metadata$inv.D), col = 2, lwd = 2)
text(3.1,440,"Pearson's r = 0.17", cex=0.9)
title(main="inverse Simpson Diversity vs Total Nitrogen (mg/L)", cex.main = 1)

#pH
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$pH, metadata$S, pch = 19, xlab= "pH",
     ylab = "Species Richness")
lines(lowess(metadata$pH, metadata$S), col = 2, lwd = 2)
text(2,1730,"Pearson's r = -0.13", cex=0.9)
title(main="Species richness vs pH", cex.main = 1)
plot(metadata$pH, metadata$H, pch = 19, xlab= "pH",
     ylab = "Shannon Diversity Index")
lines(lowess(metadata$pH, metadata$H), col = 2, lwd = 2)
text(2,6.4,"Pearson's r = -0.15", cex=0.9)
title(main="Shannon Diversity vs pH", cex.main = 1)
plot(metadata$pH, metadata$J, pch = 19, xlab= "pH",
     ylab = "Evenness")
lines(lowess(metadata$pH, metadata$J), col = 2, lwd = 2)
text(2,0.74,"Pearson's r = -0.11", cex=0.9)
title(main="Species Evenness vs pH", cex.main = 1)
plot(metadata$pH, metadata$inv.D, pch = 19, xlab= "pH",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$pH, metadata$inv.D), col = 2, lwd = 2)
text(2,400,"Pearson's r = -0.16", cex=0.9)
title(main="inverse Simpson Diversity vs pH", cex.main = 1)

#TN:TP
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$TN.TP.ratio, metadata$S, pch = 19, xlab= "TN : TP",
     ylab = "Species Richness")
lines(lowess(metadata$TN.TP.ratio, metadata$S), col = 2, lwd = 2)
text(40,1780,"Pearson's r = -0.13", cex=0.9)
title(main="Species richness vs TN : TP", cex.main = 1)
plot(metadata$TN.TP.ratio, metadata$H, pch = 19, xlab= "TN : TP",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs TN : TP", cex.main = 1)
plot(metadata$TN.TP.ratio, metadata$J, pch = 19, xlab= "TN : TP",
     ylab = "Evenness", col="grey54")
title(main="Species Evenness vs TN : TP", cex.main = 1)
plot(metadata$TN.TP.ratio, metadata$inv.D, pch = 19, xlab= "TN : TP",
     ylab = "inverse Simpson Diversity Index", col="grey54")
title(main="inverse Simpson Diversity vs TN : TP", cex.main = 1)

#Total Phosphate
par(mfrow = c(2,2), mar = c(4, 4, 3, 3)) ## all plots on one graph
plot(metadata$Phosphate.Ortho, metadata$S, pch = 19, xlab= "Total Phosphate (mg/L)",
     ylab = "Species Richness", col="grey54")
title(main="Species richness vs Total Phosphate", cex.main = 1)
plot(metadata$Phosphate.Ortho, metadata$H, pch = 19, xlab= "Total Phosphate (mg/L)",
     ylab = "Shannon Diversity Index", col="grey54")
title(main="Shannon Diversity vs Total Phosphate", cex.main = 1)
plot(metadata$Phosphate.Ortho, metadata$J, pch = 19, xlab= "Total Phosphate (mg/L)",
     ylab = "Evenness")
lines(lowess(metadata$Phosphate.Ortho, metadata$J), col = 2, lwd = 2)
text(0.17,0.82,"Pearson's r = -0.11", cex=0.9)
title(main="Species Evenness vs Total Phosphate", cex.main = 1)
plot(metadata$Phosphate.Ortho, metadata$inv.D, pch = 19, xlab= "Total Phosphate (mg/L)",
     ylab = "inverse Simpson Diversity Index")
lines(lowess(metadata$Phosphate.Ortho, metadata$inv.D), col = 2, lwd = 2)
text(0.17,400,"Pearson's r = -0.12", cex=0.9)
title(main="inverse Simpson Diversity vs Total Phosphate", cex.main = 1)


###### Alpha Diversity vs Environmental Variables - Correlation Heat map ######
library(corrplot)
library(reshape2)

#load in metadata
metadata <- read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)

#Making a dataframe with only env. variables and a-div measures
alphaenv <- metadata[, c(7:24,33,36:38,40:42)] 
#changing some column names
colnames(alphaenv)[8] ="Pheophytin-a"
colnames(alphaenv)[9] ="Chlorophyll-a"
colnames(alphaenv)[12] ="Nitrate + Nitrite"
colnames(alphaenv)[13] ="Total.Phosphate"
colnames(alphaenv)[14] ="Total.Phosphorus"
colnames(alphaenv)[19] ="Microcystin"
colnames(alphaenv)[20] ="Anatoxin-a"

#Before making heatmap, we must first calculate the correlation coefficient 
#between each variable using cor() and then transform the results into a usable 
#format using the melt() function from the reshape2 package

#calculate correlation coefficients, rounded to 2 decimal places
envcor <- round(cor(alphaenv), 2) #this a correlation matrix
testRes <- cor.mtest(alphaenv, conf.level = 0.95) #generates a table of p-values

#creating heatmap
corrplot(envcor, 
         type = "lower",
         method = 'color',
         col = COL2('BrBG', 10),
         p.mat = testRes$p,
         insig = 'label_sig',
         pch.cex = 0.98,
         pch.col = 'grey8',
         sig.level = c(0.001, 0.01, 0.05),
         order = 'original',
         number.cex = 0.8,
         tl.col = 'black',
         cl.ratio = 0.2, 
         tl.srt = 45)


###### Beta Diversity - Creating Bray Curtis matrix #################
## re-creating relative abundance table
set.seed(1998)
dat<-read.csv("feature_Y123_ADJUSTED.csv", header=TRUE, row.names = 1)
dat<-data.matrix(dat)
typeof(dat)
dat <- t(dat)
row.names(dat)
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
typeof(metadata) 
dat <- as.data.frame(dat)
typeof(dat)
common.rownames <- intersect(rownames(dat), rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata))
otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund] 
dat.pa<-decostand(dat.dom, method ="pa")
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))
dat.01per<-dat.dom[,dat.otus.01per]
dat.otus.001per<-which(colSums(dat.pa) > (0.001*nrow(dat.pa)))
dat.001per<-dat.dom[,dat.otus.001per]
dat.ra<-decostand(dat.01per, method = "total") 

#use relative abundance table created
#creating Bray-Curtis dissimilarity distance matrix
ra.bc.dist<-vegdist(dat.ra, method = "bray")

#Separating into the different years
Y1r <- dat.ra[grep("_19$", rownames(dat.ra)),]
Y2r <- dat.ra[grep("_20$", rownames(dat.ra)),]
Y3r <- dat.ra[grep("_21$", rownames(dat.ra)),]
ra.bc.d.Y1<-vegdist(Y1r, method = "bray")
ra.bc.d.Y2<-vegdist(Y2r, method = "bray")
ra.bc.d.Y3<-vegdist(Y3r, method = "bray")
metadata <-read.csv("Metadata-Diversity_BATCH.csv", row.names = 1)



###### Plotting NMDS by Year - 2D ######
nmds2d <- metaMDS(ra.bc.dist,k=2,autotransform = F,trymax=20)
#Dimensions = 2
#Stress = 0.1705273 
stressplot(nmds2d)
#Shepard plot "shows scatter around the regression between the inter-point 
#distances in the final configuration (i.e., the distances between each pair of communities) 
#against their original dissimilarities"


#Fitting environmental vectors to NMDS plot
ef.cca<- envfit(cca.p,metadata[,c(7,8,16)])
ef.cca$vectors$pvals

nmds.plot <- ordiplot(nmds2d,display="sites")
## Adding ellipses to group years
ordihull(nmds.plot,groups=metadata$Year,draw="lines",col=c("tomato3","steelblue3","springgreen3"))
##adjust colors to match each year, pch=20 makes it bullet points 
points(nmds.plot,"sites", pch=20, col= "tomato4", select = metadata$Year == "1")
points(nmds.plot,"sites", pch=20, col= "steelblue4", select = metadata$Year == "2")
points(nmds.plot,"sites", pch=20, col= "springgreen4", select = metadata$Year == "3")
##Add Stress Value
text(1.2,1.5,"2D Stress: 0.17", cex=0.9)
##Adding legend
legend("topleft",legend= c("Year 1","Year 2", "Year 3"), 
       title = "Year",
       col=c("tomato4","steelblue4","springgreen4"), 
       pch=19, cex=1)
##Adding title
title(main="nMDS of Relative Abundances by Year")
#NMDS by Season
nmds.plot <- ordiplot(nmds2d,display="sites")
ordihull(nmds.plot,groups=metadata$Season,draw="lines",col = c("sienna4","royalblue3"))
points(nmds.plot,"sites", pch=20, col= "sienna4", select = metadata$Season == "dry")
points(nmds.plot,"sites", pch=20, col= "royalblue3", select = metadata$Season == "wet")
text(1.2,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("dry","wet"), 
       title = "Season",
       col=c("sienna4","royalblue3"), 
       pch=19, cex=1)
title(main="nMDS of Relative Abundances by Season") 
#NMDS by Zone
nmds.plot <- ordiplot(nmds2d,display="sites")
ordihull(nmds.plot,groups=metadata$Zone,draw="lines",col = c("palegreen3","wheat4","cornflowerblue","violetred2"))
points(nmds.plot,"sites", pch=20, col= "palegreen3", select = metadata$Zone == "Inflow")
points(nmds.plot,"sites", pch=20, col= "wheat4", select = metadata$Zone == "Nearshore")
points(nmds.plot,"sites", pch=20, col= "cornflowerblue", select = metadata$Zone == "Pelagic")
points(nmds.plot,"sites", pch=20, col= "violetred2", select = metadata$Zone == "S79")
text(1.2,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("Inflow","Nearshore","Pelagic", "S79"), 
       title = "Zone",
       col=c("palegreen3","wheat4","cornflowerblue","violetred2"), 
       pch=19, cex=1)
title(main="nMDS of Relative Abundances by Zone")
#NMDS by Month
nmds.plot <- ordiplot(nmds2d,display="sites")
ordihull(nmds.plot,groups=metadata$Month,draw="lines",col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2",
                                                                       "mediumpurple2","hotpink","tan","sienna","purple4"))
points(nmds.plot,"sites", pch=19, col= "firebrick2", select = metadata$Month == "1")
points(nmds.plot,"sites", pch=19, col= "darkorange1", select = metadata$Month == "2")
points(nmds.plot,"sites", pch=19, col= "gray34", select = metadata$Month == "3")
points(nmds.plot,"sites", pch=19, col= "goldenrod2", select = metadata$Month == "4")
points(nmds.plot,"sites", pch=19, col= "green3", select = metadata$Month == "5")
points(nmds.plot,"sites", pch=19, col= "cadetblue2", select = metadata$Month == "6")
points(nmds.plot,"sites", pch=19, col= "dodgerblue2", select = metadata$Month == "7")
points(nmds.plot,"sites", pch=19, col= "mediumpurple2", select = metadata$Month == "8")
points(nmds.plot,"sites", pch=19, col= "hotpink", select = metadata$Month == "9")
points(nmds.plot,"sites", pch=19, col= "tan", select = metadata$Month == "10")
points(nmds.plot,"sites", pch=19, col= "sienna", select = metadata$Month == "11")
points(nmds.plot,"sites", pch=19, col= "purple4", select = metadata$Month == "12")
text(1.8,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("1","2","3","4","5", "6","7","8","9","10","11","12"), title = "Month",
       col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2",
                         "mediumpurple2","hotpink","tan","sienna","purple4"), pch=19,ncol=2, cex=0.88)
title(main="nMDS of Relative Abundances by Month")
#NMDS by Station
nmds.plot <- ordiplot(nmds2d,display="sites")
ordihull(nmds.plot,groups=metadata$Station,draw="lines",col=c("#A6CEE3","#579CC7","#3688AD",
                                                                      "#8BC395","#89CB6C",
                                                                      "#40A635","#919D5F",
                                                                      "#F99392","#EB494A",
                                                                      "#E83C2D","#F79C5D",
                                                                      "#FDA746","#FE8205",
                                                                      "#E39970", "#BFA5CF",
                                                                      "#8861AC","#917099",
                                                                      "#E7E099","#DEB969",
                                                                      "#B15928"))
points(nmds.plot,"sites", pch=19, col= "#A6CEE3", select = metadata$Station == "CLV10A")
points(nmds.plot,"sites", pch=19, col= "#579CC7", select = metadata$Station == "KISSR0.0")
points(nmds.plot,"sites", pch=19, col= "#3688AD", select = metadata$Station == "L001")
points(nmds.plot,"sites", pch=19, col= "#8BC395", select = metadata$Station == "L004")
points(nmds.plot,"sites", pch=19, col= "#89CB6C", select = metadata$Station == "L005")
points(nmds.plot,"sites", pch=19, col= "#40A635", select = metadata$Station == "L006")
points(nmds.plot,"sites", pch=19, col= "#919D5F", select = metadata$Station == "L007")
points(nmds.plot,"sites", pch=19, col= "#F99392", select = metadata$Station == "L008")
points(nmds.plot,"sites", pch=19, col= "#EB494A", select = metadata$Station == "LZ2")
points(nmds.plot,"sites", pch=19, col= "#E83C2D", select = metadata$Station == "LZ25A")
points(nmds.plot,"sites", pch=19, col= "#F79C5D", select = metadata$Station == "LZ30")
points(nmds.plot,"sites", pch=19, col= "#FDA746", select = metadata$Station == "LZ40")
points(nmds.plot,"sites", pch=19, col= "#FE8205", select = metadata$Station == "PALMOUT")
points(nmds.plot,"sites", pch=19, col= "#E39970", select = metadata$Station == "PELBAY3")
points(nmds.plot,"sites", pch=19, col= "#BFA5CF", select = metadata$Station == "POLE3S")
points(nmds.plot,"sites", pch=19, col= "#8861AC", select = metadata$Station == "POLESOUT")
points(nmds.plot,"sites", pch=19, col= "#917099", select = metadata$Station == "RITTAE2")
points(nmds.plot,"sites", pch=19, col= "#E7E099", select = metadata$Station == "S308")
points(nmds.plot,"sites", pch=19, col= "#DEB969", select = metadata$Station == "S77")                                              
points(nmds.plot,"sites", pch=19, col= "#B15928", select = metadata$Station == "S79")
text(1.8,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("CLV10A","KISSR0.0","L001","L004","L005","L006","L007",
                           "L008","LZ2","LZ25A","LZ30","LZ40","PALMOUT","PELBAY3",
                           "POLE3S","POLESOUT","RITTAE2","S308","S77","S79"),title = "Station",
       col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F", 
                      "#F99392","#EB494A","#E83C2D","#F79C5D","#FDA746","#FE8205",
                      "#E39970","#BFA5CF","#8861AC","#917099","#E7E099","#DEB969",
                      "#B15928"),ncol=2,pch=19, cex=0.74)
title(main="nMDS of Relative Abundances by Station")

#Statistics
anosim(ra.bc.dist, metadata$Year, permutations = 999, distance = "bray")
# ANOSIM statistic R: -0.003354 
# Significance: 0.748 -> NOT SIGNIFICANT
anosim(ra.bc.dist, metadata$Season, permutations = 999, distance = "bray")
# ANOSIM statistic R: -0.004122 
# Significance: 0.78  -> NOT SIGNIFICANT
anosim(ra.bc.dist, metadata$Month, permutations = 999, distance = "bray")
# ANOSIM statistic R: -0.00777  
# Significance: 0.913 -> NOT SIGNIFICANT
anosim(ra.bc.dist, metadata$Zone, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.01493
# Significance: 0.191 -> NOT SIGNIFICANT
anosim(ra.bc.dist, metadata$Station, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.1967
# Significance: 0.001

###### Plotting NMDS separated by Year - 2D ONLY ####
###Year 1
nmdsY1 <- metaMDS(ra.bc.d.Y1,k=2,autotransform = F,trymax=20)
# Dimensions: 2 
# Stress: 0.1672539 
stressplot(nmdsY1)
#Base Plot and title
nmds.plot.Y1 <- ordiplot(nmdsY1,display="sites")
title(main="nMDS of Relative Abundances - Year 1")
text(1,1.5,"2D Stress: 0.17", cex=0.9)

#Month
ordihull(nmds.plot.Y1,groups=met1$Month,draw="lines",col=c("gray34","goldenrod2","green3","cadetblue2","dodgerblue2",
                                                                     "mediumpurple2","hotpink","tan","sienna","purple4"))
points(nmds.plot.Y1,"sites", pch=19, col= "gray34", select = met1$Month == "3")
points(nmds.plot.Y1,"sites", pch=19, col= "goldenrod2", select = met1$Month == "4")
points(nmds.plot.Y1,"sites", pch=19, col= "green3", select = met1$Month == "5")
points(nmds.plot.Y1,"sites", pch=19, col= "cadetblue2", select = met1$Month == "6")
points(nmds.plot.Y1,"sites", pch=19, col= "dodgerblue2", select = met1$Month == "7")
points(nmds.plot.Y1,"sites", pch=19, col= "mediumpurple2", select = met1$Month == "8")
points(nmds.plot.Y1,"sites", pch=19, col= "hotpink", select = met1$Month == "9")
points(nmds.plot.Y1,"sites", pch=19, col= "tan", select = met1$Month == "10")
points(nmds.plot.Y1,"sites", pch=19, col= "sienna", select = met1$Month == "11")
points(nmds.plot.Y1,"sites", pch=19, col= "purple4", select = met1$Month == "12")
text(1,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("3","4","5", "6","7","8","9","10","11","12"), 
       title = "Month",ncol=2, col=c("gray34","goldenrod2","green3","cadetblue2",
                                     "dodgerblue2","mediumpurple2","hotpink",
                                     "tan","sienna","purple4"), 
                                                                                          pch=19, cex=0.92)
title(main="nMDS of Relative Abundances by Month - Year 1")
                                                                     
                                                                     
#Season
nmds.plot.Y1 <- ordiplot(nmdsY1,display="sites")
ordihull(nmds.plot.Y1,groups=met1$Season,draw="lines",col = c("sienna4","royalblue3"))
points(nmds.plot.Y1,"sites", pch=19, col= "sienna4", select = met1$Season == "dry")
points(nmds.plot.Y1,"sites", pch=19, col= "royalblue3", select = met1$Season == "wet")
text(1,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("dry","wet"), title = "Season", 
       col=c("sienna4","royalblue3"), pch=19, cex=0.92)
title(main="nMDS of Relative Abundances by Season - Year 1")
                                                                     
#Zone
nmds.plot.Y1 <- ordiplot(nmdsY1,display="sites")
ordihull(nmds.plot.Y1,groups=met1$Zone,draw="lines",col = c("palegreen3","wheat4","cornflowerblue","violetred2"))
points(nmds.plot.Y1,"sites", pch=19, col= "palegreen3", select = met1$Zone == "Inflow")
points(nmds.plot.Y1,"sites", pch=19, col= "wheat4", select = met1$Zone == "Nearshore")
points(nmds.plot.Y1,"sites", pch=19, col= "cornflowerblue", select = met1$Zone == "Pelagic")
points(nmds.plot.Y1,"sites", pch=19, col= "violetred2", select = met1$Zone == "S79")
text(1,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("Inflow","Nearshore","Pelagic", "S79"),
       title = "Zone",col=c("palegreen3","wheat4","mediumblue","violetred2"), 
       pch=19, cex=0.92)
title(main="nMDS of Relative Abundances by Zone - Year 1")
                                                                     
#Station
nmds.plot.Y1 <- ordiplot(nmdsY1,display="sites")
ordihull(nmds.plot.Y1,groups=met1$Station,draw="lines",col=c("#A6CEE3","#579CC7","#3688AD",
                                                                        "#8BC395","#89CB6C",
                                                                        "#40A635","#919D5F",
                                                                        "#F99392","#EB494A",
                                                                        "#E83C2D","#F79C5D",
                                                                        "#FDA746","#FE8205",
                                                                        "#E39970", "#BFA5CF",
                                                                        "#8861AC","#917099",
                                                                        "#E7E099","#DEB969",
                                                                        "#B15928"))
points(nmds.plot.Y1,"sites", pch=19, col= "#A6CEE3", select = met1$Station == "CLV10A")
points(nmds.plot.Y1,"sites", pch=19, col= "#579CC7", select = met1$Station == "KISSR0.0")
points(nmds.plot.Y1,"sites", pch=19, col= "#3688AD", select = met1$Station == "L001")
points(nmds.plot.Y1,"sites", pch=19, col= "#8BC395", select = met1$Station == "L004")
points(nmds.plot.Y1,"sites", pch=19, col= "#89CB6C", select = met1$Station == "L005")
points(nmds.plot.Y1,"sites", pch=19, col= "#40A635", select = met1$Station == "L006")
points(nmds.plot.Y1,"sites", pch=19, col= "#919D5F", select = met1$Station == "L007")
points(nmds.plot.Y1,"sites", pch=19, col= "#F99392", select = met1$Station == "L008")
points(nmds.plot.Y1,"sites", pch=19, col= "#EB494A", select = met1$Station == "LZ2")
points(nmds.plot.Y1,"sites", pch=19, col= "#E83C2D", select = met1$Station == "LZ25A")
points(nmds.plot.Y1,"sites", pch=19, col= "#F79C5D", select = met1$Station == "LZ30")
points(nmds.plot.Y1,"sites", pch=19, col= "#FDA746", select = met1$Station == "LZ40")
points(nmds.plot.Y1,"sites", pch=19, col= "#FE8205", select = met1$Station == "PALMOUT")
points(nmds.plot.Y1,"sites", pch=19, col= "#E39970", select = met1$Station == "PELBAY3")
points(nmds.plot.Y1,"sites", pch=19, col= "#BFA5CF", select = met1$Station == "POLE3S")
points(nmds.plot.Y1,"sites", pch=19, col= "#8861AC", select = met1$Station == "POLESOUT")
points(nmds.plot.Y1,"sites", pch=19, col= "#917099", select = met1$Station == "RITTAE2")
points(nmds.plot.Y1,"sites", pch=19, col= "#E7E099", select = met1$Station == "S308")
points(nmds.plot.Y1,"sites", pch=19, col= "#DEB969", select = met1$Station == "S77")                                              
points(nmds.plot.Y1,"sites", pch=19, col= "#B15928", select = met1$Station == "S79")
text(1,1.5,"2D Stress: 0.17", cex=0.9)
legend("topleft",legend= c("CLV10A","KISSR0.0","L001","L004","L005","L006","L007",
                           "L008","LZ2","LZ25A","LZ30","LZ40","PALMOUT","PELBAY3",
                           "POLE3S","POLESOUT","RITTAE2","S308","S77","S79"),title = "Station",
       col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F", 
                      "#F99392","#EB494A","#E83C2D","#F79C5D","#FDA746","#FE8205",
                      "#E39970","#BFA5CF","#8861AC","#917099","#E7E099","#DEB969",
                      "#B15928"),ncol=2,pch=19, cex=0.72)
title(main="nMDS of Relative Abundances by Station - Year 1")
                                                                                                                                             
### Year 2
nmdsY2 <- metaMDS(ra.bc.d.Y2,k=2,autotransform = F,trymax=20)
# Dimensions: 2 
# Stress: 0.1773041  
stressplot(nmdsY2)
#Base Plot and title
nmds.plot.Y2 <- ordiplot(nmdsY2,display="sites")
title(main="nMDS of Relative Abundances - Year 2")


#Month
nmds.plot.Y2 <- ordiplot(nmdsY2,display="sites")
ordihull(nmds.plot.Y2,groups=met2$Month,draw="lines",col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2",
                                                                         "mediumpurple2","hotpink","tan","sienna","purple4"))
points(nmds.plot.Y2,"sites", pch=19, col= "firebrick2", select = met2$Month == "1")
points(nmds.plot.Y2,"sites", pch=19, col= "darkorange1", select = met2$Month == "2")
points(nmds.plot.Y2,"sites", pch=19, col= "gray34", select = met2$Month == "3")
points(nmds.plot.Y2,"sites", pch=19, col= "goldenrod2", select = met2$Month == "4")
points(nmds.plot.Y2,"sites", pch=19, col= "green3", select = met2$Month == "5")
points(nmds.plot.Y2,"sites", pch=19, col= "cadetblue2", select = met2$Month == "6")
points(nmds.plot.Y2,"sites", pch=19, col= "dodgerblue2", select = met2$Month == "7")
points(nmds.plot.Y2,"sites", pch=19, col= "mediumpurple2", select = met2$Month == "8")
points(nmds.plot.Y2,"sites", pch=19, col= "hotpink", select = met2$Month == "9")
points(nmds.plot.Y2,"sites", pch=19, col= "tan", select = met2$Month == "10")
points(nmds.plot.Y2,"sites", pch=19, col= "sienna", select = met2$Month == "11")
points(nmds.plot.Y2,"sites", pch=19, col= "purple4", select = met2$Month == "12")
text(1.8,1.4,"2D Stress: 0.18", cex=0.9)
legend("topleft",legend= c("1","2","3","4","5", "6","7","8","9","10","11","12"), title = "Month",
       col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2",
                         "mediumpurple2","hotpink","tan","sienna","purple4"), pch=19,ncol=2, cex=0.88)
title(main="nMDS of Relative Abundances by Month - Year 2")
                                                                         
#Season
nmds.plot.Y2 <- ordiplot(nmdsY2,display="sites")
ordihull(nmds.plot.Y2,groups=met2$Season,draw="lines",col = c("sienna4","royalblue3"))
points(nmds.plot.Y2,"sites", pch=19, col= "sienna4", select = met2$Season == "dry")
points(nmds.plot.Y2,"sites", pch=19, col= "royalblue3", select = met2$Season == "wet")
text(1.8,1.4,"2D Stress: 0.18", cex=0.9)
legend("topleft",legend= c("dry","wet"),  title = "Season",col=c("sienna4","royalblue3"),  pch=19, cex=0.92)
title(main="nMDS of Relative Abundances by Season - Year 2")
                                                                         
#Zone
nmds.plot.Y2 <- ordiplot(nmdsY2,display="sites")
ordihull(nmds.plot.Y2,groups=met2$Zone,draw="lines",col = c("palegreen3","wheat4","cornflowerblue","violetred2"))
points(nmds.plot.Y2,"sites", pch=19, col= "palegreen3", select = met2$Zone == "Inflow")
points(nmds.plot.Y2,"sites", pch=19, col= "wheat4", select = met2$Zone == "Nearshore")
points(nmds.plot.Y2,"sites", pch=19, col= "cornflowerblue", select = met2$Zone == "Pelagic")
points(nmds.plot.Y2,"sites", pch=19, col= "violetred2", select = met2$Zone == "S79")
text(1.8,1.4,"2D Stress: 0.18", cex=0.9)
legend("topleft",legend= c("Inflow","Nearshore","Pelagic", "S79"), title = "Zone",
       col=c("palegreen3","wheat4","mediumblue","violetred2"), pch=19, cex=0.92)
title(main="nMDS of Relative Abundances by Zone - Year 2")
                                                                         
#Station
nmds.plot.Y2 <- ordiplot(nmdsY2,display="sites")
ordihull(nmds.plot.Y2,groups=met2$Station,draw="lines",col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F","#F99392","#EB494A","#E83C2D","#F79C5D","#FDA746","#FE8205","#E39970","#BFA5CF","#8861AC","#917099","#E7E099","#DEB969","#B15928"))
points(nmds.plot.Y2,"sites", pch=19, col= "#A6CEE3", select = met2$Station == "CLV10A")
points(nmds.plot.Y2,"sites", pch=19, col= "#579CC7", select = met2$Station == "KISSR0.0")
points(nmds.plot.Y2,"sites", pch=19, col= "#3688AD", select = met2$Station == "L001")
points(nmds.plot.Y2,"sites", pch=19, col= "#8BC395", select = met2$Station == "L004")
points(nmds.plot.Y2,"sites", pch=19, col= "#89CB6C", select = met2$Station == "L005")
points(nmds.plot.Y2,"sites", pch=19, col= "#40A635", select = met2$Station == "L006")
points(nmds.plot.Y2,"sites", pch=19, col= "#919D5F", select = met2$Station == "L007")
points(nmds.plot.Y2,"sites", pch=19, col= "#F99392", select = met2$Station == "L008")
points(nmds.plot.Y2,"sites", pch=19, col= "#EB494A", select = met2$Station == "LZ2")
points(nmds.plot.Y2,"sites", pch=19, col= "#E83C2D", select = met2$Station == "LZ25A")
points(nmds.plot.Y2,"sites", pch=19, col= "#F79C5D", select = met2$Station == "LZ30")
points(nmds.plot.Y2,"sites", pch=19, col= "#FDA746", select = met2$Station == "LZ40")
points(nmds.plot.Y2,"sites", pch=19, col= "#FE8205", select = met2$Station == "PALMOUT")
points(nmds.plot.Y2,"sites", pch=19, col= "#E39970", select = met2$Station == "PELBAY3")
points(nmds.plot.Y2,"sites", pch=19, col= "#BFA5CF", select = met2$Station == "POLE3S")
points(nmds.plot.Y2,"sites", pch=19, col= "#8861AC", select = met2$Station == "POLESOUT")
points(nmds.plot.Y2,"sites", pch=19, col= "#917099", select = met2$Station == "RITTAE2")
points(nmds.plot.Y2,"sites", pch=19, col= "#E7E099", select = met2$Station == "S308")
points(nmds.plot.Y2,"sites", pch=19, col= "#DEB969", select = met2$Station == "S77")
points(nmds.plot.Y2,"sites", pch=19, col= "#B15928", select = met2$Station == "S79")
text(1.8,1.4,"2D Stress: 0.18", cex=0.9)
legend("topleft",legend= c("CLV10A","KISSR0.0","L001","L004","L005",
                               "L006","L007","L008","LZ2","LZ25A","LZ30","LZ40",
                               "PALMOUT","PELBAY3","POLE3S","POLESOUT","RITTAE2",
                               "S308","S77","S79"),title = "Station", 
       col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F",
                      "#F99392","#EB494A","#E83C2D","#F79C5D","#FDA746","#FE8205",
                      "#E39970", "#BFA5CF","#8861AC","#917099","#E7E099","#DEB969",
                      "#B15928"),pch=19, ncol=2,cex=0.64)
title(main="nMDS of Relative Abundances by Station - Year 2")
                                                                                                                                                 
### Year 3
nmdsY3 <- metaMDS(ra.bc.d.Y3,k=2,autotransform = F,trymax=20)
# Dimensions: 2 
# Stress: 0.1471427 
stressplot(nmdsY3)
#Base Plot and title
nmds.plot.Y3 <- ordiplot(nmdsY3,display="sites")
title(main="nMDS of Relative Abundances - Year 3")



#Month
nmds.plot.Y3 <- ordiplot(nmdsY3,display="sites")
ordihull(nmds.plot.Y3,groups=met3$Month,draw="lines",col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2","mediumpurple2","hotpink","tan"))
points(nmds.plot.Y3,"sites", pch=19, col= "firebrick2", select = met3$Month == "1")
points(nmds.plot.Y3,"sites", pch=19, col= "darkorange1", select = met3$Month == "2")
points(nmds.plot.Y3,"sites", pch=19, col= "gray34", select = met3$Month == "3")
points(nmds.plot.Y3,"sites", pch=19, col= "goldenrod2", select = met3$Month == "4")
points(nmds.plot.Y3,"sites", pch=19, col= "green3", select = met3$Month == "5")
points(nmds.plot.Y3,"sites", pch=19, col= "cadetblue2", select = met3$Month == "6")
points(nmds.plot.Y3,"sites", pch=19, col= "dodgerblue2", select = met3$Month == "7")
points(nmds.plot.Y3,"sites", pch=19, col= "mediumpurple2", select = met3$Month == "8")
points(nmds.plot.Y3,"sites", pch=19, col= "hotpink", select = met3$Month == "9")
points(nmds.plot.Y3,"sites", pch=19, col= "tan", select = met3$Month == "10")
text(-0.85,1.3,"2D Stress: 0.15", cex=0.9)
legend("topright",legend= c("1","2","3","4","5", "6","7","8","9","10"),
       title = "Month",
       col=c("firebrick2","darkorange1","gray34","goldenrod2","green3","cadetblue2","dodgerblue2","mediumpurple2","hotpink","tan"), 
       pch=19, ncol=2,cex=1)
title(main="nMDS of Relative Abundances by Month - Year 3")
                                                                         
#Season
nmds.plot.Y3 <- ordiplot(nmdsY3,display="sites")
ordihull(nmds.plot.Y3,groups=met3$Season,draw="lines",col = c("sienna4","royalblue3"))
points(nmds.plot.Y3,"sites", pch=19, col= "sienna4", select = met3$Season == "dry")
points(nmds.plot.Y3,"sites", pch=19, col= "royalblue3", select = met3$Season == "wet")
text(-0.85,1.3,"2D Stress: 0.15", cex=0.9)
legend("topright",legend= c("dry","wet"), title = "Season",col=c("sienna4","royalblue3"),pch=19, cex=1.4)
title(main="nMDS of Relative Abundances by Season - Year 3")
                                                                         
#Zone
nmds.plot.Y3 <- ordiplot(nmdsY3,display="sites")
ordihull(nmds.plot.Y3,groups=met3$Zone,draw="lines",col = c("palegreen3","wheat4","cornflowerblue","violetred2"))
points(nmds.plot.Y3,"sites", pch=19, col= "palegreen3", select = met3$Zone == "Inflow")
points(nmds.plot.Y3,"sites", pch=19, col= "wheat4", select = met3$Zone == "Nearshore")
points(nmds.plot.Y3,"sites", pch=19, col= "cornflowerblue", select = met3$Zone == "Pelagic")
points(nmds.plot.Y3,"sites", pch=19, col= "violetred2", select = met3$Zone == "S79")
text(-0.85,1.3,"2D Stress: 0.15", cex=0.9)
legend("topright",legend= c("Inflow","Nearshore","Pelagic", "S79"),title = "Zone",
       col=c("palegreen3","wheat4","mediumblue","violetred2"),pch=19, cex=0.9)
title(main="nMDS of Relative Abundances by Zone - Year 3")
                                                                         
#Station
nmds.plot.Y3 <- ordiplot(nmdsY3,display="sites")
ordihull(nmds.plot.Y3,groups=met3$Station,draw="lines",col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F","#F99392","#EB494A","#E83C2D","#F79C5D","#FDA746","#FE8205","#E39970", "#BFA5CF","#8861AC","#917099","#E7E099","#DEB969","#B15928"))
points(nmds.plot.Y3,"sites", pch=19, col= "#A6CEE3", select = met3$Station == "CLV10A")
points(nmds.plot.Y3,"sites", pch=19, col= "#579CC7", select = met3$Station == "KISSR0.0")
points(nmds.plot.Y3,"sites", pch=19, col= "#3688AD", select = met3$Station == "L001")
points(nmds.plot.Y3,"sites", pch=19, col= "#8BC395", select = met3$Station == "L004")
points(nmds.plot.Y3,"sites", pch=19, col= "#89CB6C", select = met3$Station == "L005")
points(nmds.plot.Y3,"sites", pch=19, col= "#40A635", select = met3$Station == "L006")
points(nmds.plot.Y3,"sites", pch=19, col= "#919D5F", select = met3$Station == "L007")
points(nmds.plot.Y3,"sites", pch=19, col= "#F99392", select = met3$Station == "L008")
points(nmds.plot.Y3,"sites", pch=19, col= "#EB494A", select = met3$Station == "LZ2")
points(nmds.plot.Y3,"sites", pch=19, col= "#E83C2D", select = met3$Station == "LZ25A")
points(nmds.plot.Y3,"sites", pch=19, col= "#F79C5D", select = met3$Station == "LZ30")
points(nmds.plot.Y3,"sites", pch=19, col= "#FDA746", select = met3$Station == "LZ40")
points(nmds.plot.Y3,"sites", pch=19, col= "#FE8205", select = met3$Station == "PALMOUT")
points(nmds.plot.Y3,"sites", pch=19, col= "#E39970", select = met3$Station == "PELBAY3")
points(nmds.plot.Y3,"sites", pch=19, col= "#BFA5CF", select = met3$Station == "POLE3S")
points(nmds.plot.Y3,"sites", pch=19, col= "#8861AC", select = met3$Station == "POLESOUT")
points(nmds.plot.Y3,"sites", pch=19, col= "#917099", select = met3$Station == "RITTAE2")
points(nmds.plot.Y3,"sites", pch=19, col= "#E7E099", select = met3$Station == "S308")
points(nmds.plot.Y3,"sites", pch=19, col= "#DEB969", select = met3$Station == "S77")
points(nmds.plot.Y3,"sites", pch=19, col= "#B15928", select = met3$Station == "S79")
text(-0.85,1.3,"2D Stress: 0.15", cex=0.9)
legend("topright",legend= c("CLV10A","KISSR0.0","L001","L004","L005","L006",
                           "L007","L008","LZ2","LZ25A","LZ30","LZ40","PALMOUT",
                           "PELBAY3","POLE3S","POLESOUT","RITTAE2","S308","S77","S79"),
       title = "Station",col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C",
                                        "#40A635","#919D5F","#F99392","#EB494A",
                                        "#E83C2D","#F79C5D","#FDA746","#FE8205",
                                        "#E39970","#BFA5CF","#8861AC","#917099",
                                        "#E7E099","#DEB969","#B15928"), 
                                        ncol=2,pch=19, cex=0.8)
title(main="nMDS of Relative Abundances by Station - Year 3")                                                                                                                                                 

###### Beta Diversity Stat. Analyses for each year ######
##betadisper calculates dispersion (variances) within each group 

#Loading in metadata
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
#Subsetting metadata table by year
met1 <- metadata[grep("_19$", rownames(metadata)),]
met2 <- metadata[grep("_20$", rownames(metadata)),]
met3 <- metadata[grep("_21$", rownames(metadata)),]

#Year 1
dis.Z1 <-betadisper(ra.bc.d.Y1,met1$Zone)
dis.S1 <-betadisper(ra.bc.d.Y1,met1$Season)
dis.St1 <-betadisper(ra.bc.d.Y1,met1$Station)
dis.M1 <-betadisper(ra.bc.d.Y1,met1$Month)
#Year 2
dis.Z2 <-betadisper(ra.bc.d.Y2,met2$Zone)
dis.S2 <-betadisper(ra.bc.d.Y2,met2$Season)
dis.St2 <-betadisper(ra.bc.d.Y2,met2$Station)
dis.M2 <-betadisper(ra.bc.d.Y2,met2$Month)
#Year 3
dis.Z3 <-betadisper(ra.bc.d.Y3,met3$Zone)
dis.S3 <-betadisper(ra.bc.d.Y3,met3$Season)
dis.St3 <-betadisper(ra.bc.d.Y3,met3$Station)
dis.M3 <-betadisper(ra.bc.d.Y3,met3$Month)

##permutest determines if the variances differ by groups (If differences are SIGNIFICANT - use ANOSIM
##                                                        if not use PERMANOVA (adonis))
#Year 1
permutest(dis.Z1, pairwise=TRUE, permutations=999)
#            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups      3 0.0448 0.014934 1.4207    999  0.238 -> NOT SIGNIFICANT
# Residuals 153 1.6082 0.010511  
# ---

permutest(dis.S1, pairwise=TRUE, permutations=999)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      1 0.00001 0.0000127 0.0013    999  0.968 -> NOT SIGNIFICANT
# Residuals 155 1.45375 0.0093790                        
# ---
  
permutest(dis.M1, pairwise=TRUE, permutations=999)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      9 0.07056 0.0078398 0.7765    999  0.651 -> NOT SIGNIFICANT
# Residuals 147 1.48410 0.0100959                          
# ---
 
permutest(dis.St1, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     19 0.28943 0.015233 1.1881    999  0.279 -> NOT SIGNIFICANT
# Residuals 137 1.75652 0.012821


## USE PERMANOVA/adonis!!


##PERMANOVA - determining if the differences between two or more groups are significant
adonis2(ra.bc.d.Y1~met1$Station, permutations = 999)
#               Df SumOfSqs      R2      F Pr(>F)    
# met1$Station  19   10.764 0.23512 2.2165  0.001 ***
# Residual     137   35.016 0.76488                  
# Total        156   45.779 1.00000
# ___
#Pairwise perMANOVA to see what sites have the differences
Y1Stat <- pairwise.perm.manova(ra.bc.d.Y1, met1$Station,nperm = 999,p.method = "fdr")
# Get p-values in a dataframe
Y1Stp <- Y1Stat$p.value
# Convert the data to a table
m <- as.data.frame(Y1Stp)
# Plot p-values
library(gplots)
ggballoonplot(m, 
            main ="p.values", 
            xlab ="", 
            ylab="",
            label = T, label.size=0.6, #adds the p value number to the plot
            show.margins = F)
ggballoonplot(
  m, main = "Year 1 by Station - p-value comparison",
  size = "value",
  size.range = c(1, 10),
  shape = 21,
  color = "black",
  fill = "value",
  show.label = F, legend = ggplot2::lims(0.05,0.8),
  font.label = list(size = 6, color = "black"),
  rotate.x.text = TRUE,
  ggtheme = theme_minimal())
#________________________________________ 

adonis2(ra.bc.d.Y1~met1$Season, permutations = 999)
#              Df SumOfSqs      R2      F Pr(>F)
# met1$Season   1    0.244 0.00533 0.8308  0.672 -> NOT SIGNIFICANT
# Residual    155   45.535 0.99467              
# Total       156   45.779 1.00000 

adonis2(ra.bc.d.Y1~met1$Zone, permutations = 999)
#            Df SumOfSqs      R2      F Pr(>F)    
# met1$Zone   3    1.791 0.03911 2.0759  0.001 ***
# Residual  153   43.989 0.96089                  
# Total     156   45.779 1.00000 
# ___
#PerMANOVA to see what sites have the differences
Y1Zone <- pairwise.perm.manova(ra.bc.d.Y1, met1$Zone,nperm = 999,p.method = "fdr")
# Significant differences found between all zones

adonis2(ra.bc.d.Y1~met1$Month, permutations = 999)
#             Df SumOfSqs      R2      F Pr(>F)
# met1$Month   1    0.157 0.00342 0.5322  0.994 -> NOT SIGNIFICANT
# Residual   155   45.622 0.99658              
# Total      156   45.779 1.00000   


#Year 2
permutest(dis.Z2, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)   
# Groups      3 0.17468 0.058226 6.558    999  0.002 **
# Residuals 206 1.82900 0.008879                       
# ---
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#              Inflow  Nearshore    Pelagic   S79
# Inflow               2.2100e-01 3.1000e-02 0.018
# Nearshore 2.2085e-01            3.6200e-01 0.002
# Pelagic   1.9483e-02 3.3873e-01            0.001
# S79       2.3672e-02 8.2715e-04 3.5696e-05    
permutest(dis.S2, pairwise=TRUE, permutations=999)
#            Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups      1 0.00219 0.0021948 0.258    999  0.614
# Residuals 208 1.76932 0.0085063                      
# ---
permutest(dis.M2, pairwise=TRUE, permutations=999)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     11 0.05232 0.0047561 0.5497    999  0.858
# Residuals 198 1.71297 0.0086514
permutest(dis.St2, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     19 0.67528 0.035541 2.8946    999  0.001 ***
# Residuals 190 2.33290 0.012278   

## USE ANOSIM FOR ZONE AND STATION, USE PERMANOVA FOR SEASON AND MONTH!!


##ANOSIM - determining if the differences between two or more groups are significant
anosim(ra.bc.d.Y2,met2$Zone, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.01148 
# Significance: 0.314 -> NOT SIGINFICANT

anosim(ra.bc.d.Y2,met2$Station, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.2535 
# Significance: 0.001
Y2Stat <- pairwise.perm.manova(ra.bc.d.Y2, met2$Station,nperm = 999,p.method = "fdr")

##PERMANOVA
adonis2(ra.bc.d.Y2~met2$Month, permutations = 999)
#             Df SumOfSqs    R2      F Pr(>F)
# met2$Month   1    0.184 0.003 0.6265  0.945 -> NOT SIGINFICANT
# Residual   208   61.122 0.997              
# Total      209   61.306 1.000 
adonis2(ra.bc.d.Y2~met2$Season, permutations = 999)
#              Df SumOfSqs      R2      F Pr(>F)
# met2$Season   1    0.172 0.00281 0.5857  0.977 -> NOT SIGINFICANT
# Residual    208   61.134 0.99719              
# Total       209   61.306 1.00000   

#Year 3
permutest(dis.Z3, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
# Groups      3 0.18912 0.063039 5.1907    999  0.007 **
# Residuals 170 2.06459 0.012145                      
# ---
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#              Inflow  Nearshore    Pelagic   S79
# Inflow               0.01000000 0.16800000 0.463
# Nearshore 0.01207560            0.00100000 0.068
# Pelagic   0.15407191 0.00012197            0.975
# S79       0.46792457 0.05697194 0.96831209  
permutest(dis.S3, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
# Groups      1 0.03421 0.034209 3.3793    999  0.074 . -> NOT SIGNIFICANT
# Residuals 172 1.74117 0.010123                        
# ---
permutest(dis.M3, pairwise=TRUE, permutations=999)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      9 0.06587 0.0073193 0.7267    999  0.721 -> NOT SIGNIFICANT 
# Residuals 164 1.65174 0.0100716 
permutest(dis.St3, pairwise=TRUE, permutations=999)
#            Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)   
# Groups     19 0.70017 0.036851 3.009    999  0.003 **
# Residuals 154 1.88604 0.012247

## USE ANOSIM FOR ZONE AND STATION, USE PERMANOVA FOR SEASON AND MONTH!!


##ANOSIM - determining if the differences between two or more groups are significant
anosim(ra.bc.d.Y3,met3$Zone, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.4239
# Significance: 0.001
Y3Zone <- pairwise.perm.manova(ra.bc.d.Y3, met3$Zone,nperm = 999,p.method = "fdr")
# Significant differences found between all zones

anosim(ra.bc.d.Y3,met3$Station, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.2877
# Significance: 0.001
Y3Stat <- pairwise.perm.manova(ra.bc.d.Y3, met3$Station,nperm = 999,p.method = "fdr")

##PERMANOVA
adonis2(ra.bc.d.Y3~met3$Season, permutations = 999)
#              Df SumOfSqs      R2      F Pr(>F)
# met3$Season   1    0.265 0.00598 1.0348   0.33 -> NOT SIGNIFICANT
# Residual    172   44.122 0.99402              
# Total       173   44.387 1.00000  
adonis2(ra.bc.d.Y3~met3$Month, permutations = 999)
#             Df SumOfSqs      R2      F Pr(>F)
# met3$Month   1    0.193 0.00434 0.7504  0.735 -> NOT SIGNIFICANT
# Residual   172   44.195 0.99566              
# Total      173   44.387 1.00000 


###### Beta Diversity - Stat. Analyses - ALL YEARS TOGETHER ######
set.seed(1998)

##betadisper calculates dispersion (variances) within each group 
#values should be non-significant in order to use PERMANOVA
dis.Zone <-betadisper(ra.bc.dist,metadata$Zone)
dis.Season <-betadisper(ra.bc.dist,metadata$Season)
dis.Year <-betadisper(ra.bc.dist,metadata$Year)
dis.Station <-betadisper(ra.bc.dist,metadata$Station)
dis.Month <-betadisper(ra.bc.dist,metadata$Month)

##permutest determines if the variances differ by groups (If differences are SIGNIFICANT - use ANOSIM
##                                                        if not use PERMANOVA (adonis))
permutest(dis.Zone, pairwise=TRUE, permutations=999)
#            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups      3 0.1605 0.053487 5.3955    999  0.001 ***
# Residuals 537 5.3235 0.009913 
# ---
# Pairwise comparisons:
# (Observed p-value below diagonal, permuted p-value above diagonal)
#              Inflow Nearshore   Pelagic   S79
# Inflow              0.0030000 0.5910000 0.051
# Nearshore 0.0025931           0.0010000 0.713
# Pelagic   0.5842551 0.0011149           0.057
# S79       0.0309406 0.7291803 0.0427081 
permutest(dis.Season, pairwise=TRUE, permutations=999)
#            Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      1 0.0038 0.0037558 0.4045    999  0.532 -> NOT SIGNIFICANT
# Residuals 539 5.0041 0.0092840                          
# ---
permutest(dis.Year, pairwise=TRUE, permutations=999)
#            Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      2 0.0042 0.0021079 0.2258    999  0.809 -> NOT SIGNIFICANT
# Residuals 538 5.0226 0.0093358    
# ---
permutest(dis.Station, pairwise=TRUE, permutations=999) #look at pairwise in R (very large)
#            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups     19 1.0197 0.053670 5.1682    999  0.001 ***
# Residuals 521 5.4105 0.010385     
permutest(dis.Month, pairwise=TRUE, permutations=999)
#            Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     11 0.0580 0.0052772 0.5639    999  0.851 -> NOT SIGNIFICANT
# Residuals 529 4.9508 0.0093589                     
# ---


## USE ANOSIM FOR ZONE AND STATION AND USE PERMANOVA FOR SEASON, YEAR, AND MONTH


##ANOSIM - determining if the differences between two or more groups are significant. 
## The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to
## the mean of ranked dissimilarities within groups. An R value close to “1" suggests 
## dissimilarity between groups while an R value close to “0” suggests an even distribution of
## high and low ranks within and between groups”
## the higher the R value, the more dissimilar your groups are in terms of microbial community composition.

anosim(ra.bc.dist, metadata$Zone, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.01493 
# Significance: 0.205 -> NOT SIGNIFICANT
anosim(ra.bc.dist, metadata$Station, permutations = 999, distance = "bray")
# ANOSIM statistic R: 0.1967 
# Significance: 0.001 

##PERMANOVA
adonis2(ra.bc.dist~metadata$Month, permutations = 999)
#                 Df SumOfSqs      R2     F Pr(>F)
# metadata$Month   1    0.195 0.00127 0.683  0.909 -> NOT SIGNIFICANT
# Residual       539  154.113 0.99873             
# Total          540  154.309 1.00000  
adonis2(ra.bc.dist~metadata$Year, permutations = 999)
#                Df SumOfSqs      R2      F Pr(>F)
# metadata$Year   1    0.171 0.00111 0.5987  0.974 -> NOT SIGNIFICANT
# Residual      539  154.137 0.99889              
# Total         540  154.309 1.00000 
adonis2(ra.bc.dist~metadata$Season, permutations = 999)
#                  Df SumOfSqs      R2      F Pr(>F)
# metadata$Season   1    0.204 0.00132 0.7127  0.881 -> NOT SIGNIFICANT
# Residual        539  154.105 0.99868              
# Total           540  154.309 1.00000    

## USE MANTEL TEST FOR CONTINUOUS VARIABLES
##Mantel tests are correlation tests that determine the correlation between two 
##matrices (rather than two variables). A significant Mantel test will tell you
##that the distances between samples in one matrix are correlated with the distances 
##between samples in the other matrix. Therefore, as the distance between samples 
##increases with respect to one matrix, the distances between the same samples also 
##increases in the other matrix

#abundance dissim. matrix
dist.abund <- ra.bc.dist
#Microcystis/Bloom distance using euclidean
MA <- metadata$Microcystis.Abundance
CHL <- metadata$Chlorophyll.a
dist.MA <- dist(MA, method = "euclidean")
dist.CHL <- dist(CHL, method = "euclidean")

#Mantel test - Microcystis
mantel(dist.abund, dist.MA, method = "spearman", permutations = 999)
# Mantel statistic r: 0.008024
# Significance: 0.4 -> NOT SIGINIFCANT

#Mantel test -  Chlorophyll a
mantel(dist.abund, dist.CHL, method = "spearman", permutations = 999)
# Mantel statistic r: 0.01756 
# Significance: 0.225 -> NOT SIGNIFICANT

##Plotting beta diversity against significant variables 
#create vectors of matrices
cc <- as.vector(dist.CHL)
mm <- as.vector(dist.MA)
aa <- as.vector(dist.abund)
#new data frame with vectorized distance matrices
mat <- data.frame(cc,aa,mm)
#PLOT - Chlorophyll a
ggplot(mat, aes(y = aa, x = cc)) + 
  geom_point(size = 2, alpha = 0.75, colour = "black",shape = 21) + 
  labs(x = "Chlorophyll a (ug/L)", y = "Bray-Curtis Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
#PLOT - Microcystis
ggplot(mat, aes(y = aa, x = mm)) + 
  geom_point(size = 2, alpha = 0.75, colour = "black",shape = 21) +
  labs(x = "Microcystis Relative Abundance", y = "Bray-Curtis Dissimilarity") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))



###### Venn Diagram of ASVs (Year, Zone, Season) ######
##Packages 
library(eulerr)
library(microbiome)
library(microbiomeutilities)
#library(devtools) ##used to install microbiome utilities package
#devtools::install_github('microsud/microbiomeutilities') ## only run if need to install package



## Making phyloseq objects (WHOLE DATA SET)
asvdat <- as.data.frame(t(dat.01per)) #species has to be rows so the df was transformed
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1)
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the data will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
pseq <- phyloseq(ASV,TAX,META)


# simple way to count number of samples in each group
table(meta(pseq)$Year, useNA = "always")
## 
##   1       2        3    <NA> 
##  157     210     174      0
table(meta(pseq)$Zone, useNA = "always")
## 
# Inflow  Nearshore   Pelagic     S79       <NA> 
#   107       131       281        22         0 
table(meta(pseq)$Season, useNA = "always")
## 
##  dry     wet   <NA> 
##  247     294     0

#convert to relative abundance
transform <- microbiome::transform
pseq_rel <- transform(pseq, "compositional")

#Make a list of Years
years <- unique(as.character(meta(pseq_rel)$Year))
print(years)
# [1] "1" "2" "3"

#Make a list of Zones
zones <- unique(as.character(meta(pseq_rel)$Zone))
print(zones)
# [1] "Inflow"    "Pelagic"   "Nearshore" "S79"  

#Make a list of Seasons
seasons <- unique(as.character(meta(pseq_rel)$Season))
print(seasons)
# [1] "dry" "wet"

#### YEAR
#Write a for loop to go through each of the years 
#one by one and combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in years){ # for each variable n in Year
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq_rel, Year == n) # Choose sample from Year by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each year.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
# [1] "No. of core taxa in 1 : 14"
# [1] "No. of core taxa in 2 : 16"        WHOLE DATASET
# [1] "No. of core taxa in 3 : 32"


##Adding taxa information
print(list_core) # can see that its the ASV id w/ NO taxa info

taxa_names(pseq_rel)[1:5] #shows ASV id
# [1] "0885965c051f3034c0e28043193bc5d2" "51e00e866016fba8a19581249b811ec4"
# [3] "dfd3874c0e70ae177e8cdc4fb6961e7d" "ac879ef0bc703ee2637bc55f0ef97afc"
# [5] "41714fa1a258e8098d51d03a1e1b3304"

#format names and checking
pseq_rel_f <- format_to_besthit(pseq_rel)
taxa_names(pseq_rel_f)[1:5]

#rerun 'for' loop with better taxa information
for (n in years){ 
  ps.sub <- subset_samples(pseq_rel_f, Year == n)
  core_m <- core_members(ps.sub, 
                         detection = 0.001,
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  list_core[[n]] <- core_m 
}
print(list_core) #shows ASV id with taxa information
#converting lists to dfs and saving as CSVs
Year1VennTaxa <- as.data.frame(list_core[["1"]])
Year2VennTaxa <- as.data.frame(list_core[["2"]])
Year3VennTaxa <- as.data.frame(list_core[["3"]])
write.csv(Year1VennTaxa, "CoreTaxaYear1-Venn.csv")
write.csv(Year2VennTaxa, "CoreTaxaYear2-Venn.csv")
write.csv(Year3VennTaxa, "CoreTaxaYear3-Venn.csv")

###Comparing venn diagram packages to see which to use (1.31.23)
##Plotting venn diagram using eulerr 
plot(venn(list_core),fills = c("tomato3", "steelblue3", "springgreen3"))

#### ZONE
list_core <- c()
for (n in zones){
  ps.sub <- subset_samples(pseq_rel_f, Zone == n)
  core_m <- core_members(ps.sub,
                         detection = 0.001,
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
# [1] "No. of core taxa in Inflow : 15"
# [1] "No. of core taxa in Pelagic : 45"
# [1] "No. of core taxa in Nearshore : 31"
# [1] "No. of core taxa in S79 : 33"

print(list_core) #shows ASV id with taxa information

#converting lists to dfs and saving as CSVs
InflowVennTaxa <- as.data.frame(list_core[["Inflow"]])
NearVennTaxa <- as.data.frame(list_core[["Nearshore"]])
PelVennTaxa <- as.data.frame(list_core[["Pelagic"]])
S79VennTaxa <- as.data.frame(list_core[["S79"]])
write.csv(InflowVennTaxa, "CoreTaxaInflow-Venn.csv")
write.csv(NearVennTaxa, "CoreTaxaNear-Venn.csv")
write.csv(PelVennTaxa, "CoreTaxaPelagic-Venn.csv")
write.csv(S79VennTaxa, "CoreTaxaS79-Venn.csv")

##Plotting venn diagram
plot(venn(list_core),fills = c("palegreen3","cornflowerblue","wheat4","violetred2"))

##Plotting venn diagram using VennDiagram
#downfall - creates a png file for the venn diagram BUT there is a workaround to view it in R
#         - does not allow for less than 4 variables
install.packages("VennDiagram")
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(
  list_core,
  category.names = c("Inflow" , "Pelagic" , "Nearshore", "S79"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("palegreen3","cornflowerblue","wheat4","violetred2"),
  # Numbers
  cex = 1,
  # Set names
  cat.cex = 1.26,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1)
)

#### SEASON
list_core <- c()
for (n in seasons){
  ps.sub <- subset_samples(pseq_rel_f, Season == n)
  
  core_m <- core_members(ps.sub,
                         detection = 0.001,
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
}
# [1] "No. of core taxa in dry : 29"
# [1] "No. of core taxa in wet : 17"

print(list_core) #shows ASV id with taxa information
#converting lists to dfs and saving as CSVs
DryVennTaxa <- as.data.frame(list_core[["dry"]])
WetVennTaxa <- as.data.frame(list_core[["wet"]])
write.csv(DryVennTaxa, "CoreTaxaDry-Venn.csv")
write.csv(WetVennTaxa, "CoreTaxaWet-Venn.csv")

##Plotting venn diagram
plot(venn(list_core),fills = c("lemonchiffon2","royalblue1"))

##Core line plots
# Determine core microbiota across various abundance/prevalence thresholds with 
# the blanket analysis (Salonen et al. CMI, 2012) based on various signal and 
# prevalences.

# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(pseq_rel_f, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

##Core heatmaps
# This visualization method has been used for instance in Intestinal microbiome 
# landscaping: Insight in community assemblage and implications for microbial 
# modulation strategies. Shetty et al. FEMS Microbiology Reviews fuw045, 2017.

#Note that you can order the taxa on the heatmap with the order.taxa argument.

# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

#Deletes "ASV" from taxa_names, e.g. ASV1 --> 1
#taxa_names(ps.m3.rel) = taxa_names(ps.m3.rel) %>% str_replace("ASV", "")
# Also define gray color palette
gray <- gray(seq(0,1,length=5))

p1 <- plot_core(pseq_rel_f,
                plot.type = "heatmap",
                colours = gray,
                prevalences = prevalences,
                detections = detections, min.prevalence = .05) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1


###### CCA Analysis - Overall and Year-to-Year ######
set.seed(1998)
#ALL YEARS
ccamodel <- cca(dat.ra~., metadata[,c(7:37)]) #run 1
# If VIF>10, the variable presents colinearity with another or other variables. 
# In that case, delete the variable from initial dataset and redo the analysis.
# VIF = 1 for completely independent variables,and values above 10 or 20 
# (depending on your taste) are regarded as highly multicollinear (dependent on others).

ccamodel <- cca(dat.ra~., metadata[,c(7:19,21:24,31,33)]) #run 2
anova.cca(finalmodel, by="terms")
#                          Df ChiSquare       F Pr(>F)    
#   SecchiDiskDepth         1    0.1574  9.9667  0.001 ***
#   Silica                  1    0.0667  4.2218  0.001 ***
#   Sulfate                 1    0.0552  3.4962  0.001 ***
#   Temperature             1    0.1163  7.3647  0.001 ***
#   Turbidity               1    0.1578  9.9912  0.001 ***
#   Alkalinity              1    0.1466  9.2843  0.001 ***
#   Ammonia                 1    0.1299  8.2251  0.006 ** 
#   Pheophytin.a            1    0.0678  4.2934  0.001 ***
#   Chlorophyll.a           1    0.1273  8.0613  0.001 ***
#   TotalDepth              1    0.0952  6.0274  0.001 ***
#   DissolvedOxygen         1    0.0584  3.6952  0.004 ** 
#   Nitrate.Nitrite         1    0.0654  4.1389  0.001 ***
#   Phosphate.Ortho         1    0.0530  3.3573  0.001 ***
#   pH                      1    0.0321  2.0298  0.014 *  
#   Total.Nitrogen          1    0.0360  2.2828  0.004 ** 
#   TN.TP.ratio             1    0.0677  4.2882  0.001 ***
#   Microcystis.Abundance   1    0.1738 11.0048  0.001 ***
#   Microcystin.LA          1    0.0144  0.9097  0.383   -> REMOVE  
#   Microcystin.LR          1    0.0243  1.5367  0.038 *  
#   Residual              521    8.2273                 
# ---    

ccamodel <- cca(dat.ra~., metadata[,c(7:19,21:24,33)]) #run 3 
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel) ## everything is under 10
finalmodel ## Note that "Total Inertia" is the total variance in species (observations matrix) distributions. 
## "Constrained Inertia" is the variance explained by the environmental variables (gradients matrix). 
## The "Proportion" values represent the percentages of variance of species distributions explained  
## by Constrained (environmental) and Unconstrained variables. Eigenvalues of constrained and 
## unconstrained axes represent the amount of variance explained by each CCA axis (graphs usually 
## present the first two constrained axes, so take a look at their values).
#Total Inertia = total variance in species (observdistributions
#Unconstrained Inertia = the variance explained by the environmental variables

#               Inertia Proportion Rank
# Total           9.872      1.000     
# Constrained     1.629      0.165   18
# Unconstrained   8.243      0.835  522
# Inertia is scaled Chi-square

R2.adj.cca <- RsquareAdj(finalmodel) 
# adjusting the R-squared value: The adjusted R2 tells you the percentage of 
# variation explained by only the independent variables that actually affect 
# the dependent variable
#indicates how well terms fit a curve or line, but adjusts for the number of terms in a model
R2.adj.cca
# r.squared: 0.173352
# adj.r.squared: 0.1446893

# Testing the significance of the CCA model
anova.cca(finalmodel) #should be significant
#           Df ChiSquare      F Pr(>F)    
# Model     18    1.6290 5.7307  0.001 ***
# Residual 522    8.2434                  
# ---    

# Testing the significance of terms (environmental variables)
anova.cca(finalmodel, by="terms")
#                          Df ChiSquare       F Pr(>F)    
#   SecchiDiskDepth         1    0.1574  9.9663  0.001 ***
#   Silica                  1    0.0667  4.2216  0.001 ***
#   Sulfate                 1    0.0552  3.4961  0.001 ***
#   Temperature             1    0.1163  7.3644  0.001 ***
#   Turbidity               1    0.1578  9.9908  0.001 ***
#   Alkalinity              1    0.1466  9.2839  0.001 ***
#   Ammonia                 1    0.1299  8.2248  0.003 ** 
#   Pheophytin.a            1    0.0678  4.2932  0.001 ***
#   Chlorophyll.a           1    0.1273  8.0610  0.001 ***
#   TotalDepth              1    0.0952  6.0272  0.001 ***
#   DissolvedOxygen         1    0.0584  3.6951  0.002 ** 
#   Nitrate.Nitrite         1    0.0654  4.1387  0.001 ***
#   Phosphate.Ortho         1    0.0530  3.3572  0.002 ** 
#   pH                      1    0.0321  2.0297  0.008 ** 
#   Total.Nitrogen          1    0.0360  2.2827  0.003 ** 
#   TN.TP.ratio             1    0.0677  4.2880  0.001 ***
#   Microcystis.Abundance   1    0.1738 11.0044  0.001 ***
#   Microcystin.LR          1    0.0225  1.4273  0.064 .   -> Make sure to specify that it had a p-value of 0.06
#   Residual              522    8.2434                   
# ---

summary(finalmodel)

## Correlation between the significant environmental variables
cor(metadata[,c(7:19,21:24,33)], method ="pearson")
#create pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(metadata[,c(7:19,21:24,33)])



#Year-by-year
#Year 1
ccamodel <- cca(Y1r~., met1[,c(7:37)]) #run1
ccamodel <- cca(Y1r~., met1[,c(7:18,21,23,24,28)]) #run2
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel)
finalmodel
#               Inertia Proportion Rank
# Total          8.5646     1.0000     
# Constrained    2.0371     0.2379   14
# Unconstrained  6.5275     0.7621  142
# Inertia is scaled Chi-square 
# 588 species (variables) deleted due to missingness

R2.adj.cca <- RsquareAdj(finalmodel) 
R2.adj.cca
# r.squared:0.2591125
# adj.r.squared: 0.1743835

# Testing the significance of the CCA model
anova.cca(finalmodel)
#           Df ChiSquare      F Pr(>F)    
# Model     16    2.2068 3.0372  0.001 ***
# Residual 140    6.3578                  
# ---

# Testing the significance of terms (environmental variables)
anova.cca(finalmodel, by="terms")
# Microcystin             1    0.0339 0.7469  0.746  -> NOT SIG.  

#create pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(met1[,c(7:18,21,23,24)])


#Year 2
ccamodel <- cca(Y2r~., met2[,c(7:37)]) #run1
ccamodel <- cca(Y2r~., met2[,c(7:19,21:24,28,31,33,36,37)]) #run2
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel)
finalmodel
#               Inertia Proportion Rank
# Total          9.3746     1.0000     
# Constrained    2.3486     0.2505   22
# Unconstrained  7.0260     0.7495  187
# Inertia is scaled Chi-square 


R2.adj.cca <- RsquareAdj(finalmodel) 
R2.adj.cca
# r.squared:0.2593453
# adj.r.squared:0.172592


anova.cca(finalmodel)
#           Df ChiSquare      F Pr(>F)    
# Model     22    2.3486 2.8413  0.001 ***
# Residual 187    7.0260                  
# ---



#create pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(met2[,c(7:19,21:24,28,31,33,36,37)])



#Year 3
ccamodel <- cca(Y3r~., met3[,c(7:37)]) #run1
ccamodel <- cca(Y3r~., met3[,c(7:10,12:19,21,23,24,31,33)]) #run2
finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))
vif.cca(finalmodel)
finalmodel
#               Inertia Proportion Rank
# Total          6.9434     1.0000     
# Constrained    1.9044     0.2743   15
# Unconstrained  5.0390     0.7257  158
# Inertia is scaled Chi-square 
# 669 species (variables) deleted due to missingness


R2.adj.cca <- RsquareAdj(finalmodel) 
R2.adj.cca
# r.squared: 0.2852408
# adj.r.squared: 0.2068729


anova.cca(finalmodel)
#           Df ChiSquare      F Pr(>F)    
# Model     17    1.9617 3.6136  0.001 ***
# Residual 156    4.9817                 
# ---


#create pairs plot to see the correlation statistics between each variable
library(psych)
pairs.panels(met3[,c(7:10,12:19,21,23,24,31,33)])


###### Plotting CCAs ######
cca.p <- plot(finalmodel,type = "none")

#Fitting of the environmental variables to the CCA plot
ef.cca<- envfit(cca.p,met3[,c(7:10,12:19,21,23,24,31,33)])
#Creating R2 threshold for vectors (found function code on research gate)
#Function: select.envfit - Setting r2 cutoff values to display in an 
#                          ordination.r.select<-0.3 # correlation threshold, 
#                          see function below
#__FUNCTION: select.envfit__#
# function (select.envfit) filters the resulting list of function (envfit) based on their p values. This allows to display only significant values in the final plot.
# just run this
select.envfit<-function(fit, r.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$r)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$r[i]<r.select) { #Check wether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    } #close if-loop
  } #close for-loop
  return(fit) #return fit as the result of the function
} #close the function

#Running select function on actual data
ef.cca<- select.envfit(ef.cca, 0.3) #selecting from a weak positive correlation and stronger


## R2 VALUES
#All years 
# SecchiDiskDepth           Silica               Sulfate           Temperature             Turbidity 
# 0.27094972            0.07374569            0.05479263            0.12495445            0.42287517 
# Alkalinity               Ammonia          Pheophytin.a         Chlorophyll.a            TotalDepth 
# 0.25428886            0.33806540            0.05730124            0.23852181            0.21004233 
# DissolvedOxygen       Nitrate.Nitrite       Phosphate.Ortho           pH              Total.Nitrogen 
# 0.42767606            0.54789964            0.47798414            0.34217550            0.05233525 
# TN.TP.ratio   Microcystis.Abundance        Microcystin.LR 
# 0.57227444            0.03451549            0.03085789 

#Year 1
# SecchiDiskDepth           Silica               Sulfate           Temperature             Turbidity 
# 0.304765766           0.059737589           0.006162602           0.025615940           0.314931560 
# Alkalinity               Ammonia          Pheophytin.a         Chlorophyll.a            TotalDepth 
# 0.210544801           0.597196549           0.054743998           0.175220168           0.220000596 
# DissolvedOxygen       Nitrate.Nitrite             pH             TN.TP.ratio    Microcystis.Abundance 
# 0.485019703           0.462306509           0.514576526           0.652323571           0.004472664 
# Microcystin 
# 0.006837999 

#Year 2
# SecchiDiskDepth          Silica               Sulfate           Temperature             Turbidity 
# 0.18704153            0.07288965            0.14544517            0.14922633            0.52276802 
# Alkalinity               Ammonia          Pheophytin.a         Chlorophyll.a            TotalDepth 
# 0.25794197            0.35220927            0.08725580            0.35052294            0.18683669 
# DissolvedOxygen       Nitrate.Nitrite       Phosphate.Ortho           pH              Total.Nitrogen 
# 0.51390253            0.54838242            0.34838408            0.68891135            0.01746031 
# TN.TP.ratio   Microcystis.Abundance        Microcystin        Microcystin.LA        Microcystin.LR 
# 0.62322767            0.01788581            0.00223627            0.02135175            0.03884635 
# Anatoxin.a    Cylindrospermopsin 
# 0.04925972            0.03583364 

#Year 3
# SecchiDiskDepth           Silica               Sulfate           Temperature            Alkalinity 
# 0.12798686            0.14790446            0.16111518            0.36344282            0.30968020 
# Ammonia             Pheophytin.a         Chlorophyll.a            TotalDepth       DissolvedOxygen 
# 0.18427317            0.09774378            0.38622539            0.20853791            0.30090864 
# Nitrate.Nitrite       Phosphate.Ortho           pH                TN.TP.ratio   Microcystis.Abundance 
# 0.67163554            0.44076153            0.11917088            0.36285678            0.55155892 
# Microcystin.LA        Microcystin.LR 
# 0.03009204            0.38899517 

#Microcystin LR strongly correlated to Microcystis abundance so removing that vector
ef.cca$vectors$arrows["Microcystin.LR",]=NA


#Setting up base plot
#ALL Years
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()
#Year 1
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()
#Year 2
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()
#Year 3
par(mar=c(5.1, 6.1, 3.1, 4.1))
plot(finalmodel,type = "none")
abline(h = 0, v = 0, col = "white", lwd = 2)
box()


#Adding the points 

#Year
#Adding the points 
points(cca.p,"sites", pch=19, col= "goldenrod3", select = metadata$Year == "1")
points(cca.p,"sites", pch=19, col= "mediumpurple2", select = metadata$Year == "2")
points(cca.p,"sites", pch=19, col= "springgreen4", select = metadata$Year == "3")
#Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
#Add legend (click to place legend on the outside of the plot) & Title
legend(locator(1),legend=c("1","2", "3"), 
       col=c("goldenrod3","mediumpurple2", "springgreen4"), pch=19, cex=1.2, 
       title = "Year")
title(main="Years 1 - 3 (2019 - 2021)")

#Zone
points(cca.p,"sites", pch=19, col= "palegreen3", select = met3$Zone == "Inflow")
points(cca.p,"sites", pch=19, col= "cornflowerblue", select = met3$Zone == "Pelagic")
points(cca.p,"sites", pch=19, col= "wheat4", select = met3$Zone == "Nearshore")
points(cca.p,"sites", pch=19, col= "violetred2", select = met3$Zone == "S79")
#Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
#Add legend (click to place legend on the outside of the plot) & Title
legend(locator(1),legend=c("Inflow","Nearshore","Pelagic","S79"), 
       col=c("palegreen3","wheat4","cornflowerblue","violetred2"), pch=19, cex=1.2, 
       title = "Ecological Zone")
title(main="Years 1 - 3 (2019 - 2021)")
title(main="Year 1 - 2019")
title(main="Year 2 - 2020")
title(main="Year 3 - 2021")
#Season
#Adding the points 
points(cca.p,"sites", pch=19, col= "lemonchiffon3", select = met3$Season == "dry")
points(cca.p,"sites", pch=19, col= "royalblue1", select = met3$Season == "wet")
#Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
#Add legend (click to place legend on the outside of the plot) & Title
legend(locator(1),legend=c("Dry","Wet"), 
       col=c("lemonchiffon3","royalblue1"), pch=19, cex=1.2, title = "Season")
title(main="Years 1 - 3 (2019 - 2021)")
title(main="Year 1 - 2019")
title(main="Year 2 - 2020")
title(main="Year 3 - 2021")

#Month
#Adding the points 
points(cca.p,"sites", pch=19, col= "firebrick2", select = met3$Month == "1")
points(cca.p,"sites", pch=19, col= "darkorange1", select = met3$Month == "2")
points(cca.p,"sites", pch=19, col= "gray38", select = met3$Month == "3")
points(cca.p,"sites", pch=19, col= "goldenrod1", select = met3$Month == "4")
points(cca.p,"sites", pch=19, col= "green4", select = met3$Month == "5")
points(cca.p,"sites", pch=19, col= "cadetblue2", select = met3$Month == "6")
points(cca.p,"sites", pch=19, col= "dodgerblue2", select = met3$Month == "7")
points(cca.p,"sites", pch=19, col= "mediumpurple2", select = met3$Month == "8")
points(cca.p,"sites", pch=19, col= "hotpink", select = met3$Month == "9")
points(cca.p,"sites", pch=19, col= "tan", select = met3$Month == "10")
points(cca.p,"sites", pch=19, col= "saddlebrown", select = met3$Month == "11")
points(cca.p,"sites", pch=19, col= "purple4", select = met3$Month == "12")
#Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
#Add legend (click to place legend on the outside of the plot) & Title
legend(locator(1),legend= c("3","4","5","6","7","8","9","10","11","12"), 
       title = "Month",ncol = 2,
       col=c("gray34","goldenrod2","green3",
                         "cadetblue2","dodgerblue2","mediumpurple2","hotpink","tan","saddlebrown","purple4"), 
                         pch=19, cex=1.2)
legend(locator(1),legend= c("1","2","3","4","5","6","7","8","9","10","11","12"), 
       title = "Month",ncol = 2,
       col=c("firebrick2","darkorange1","gray34","goldenrod2","green3",
                         "cadetblue2","dodgerblue2","mediumpurple2","hotpink","tan","saddlebrown","purple4"), 
                         pch=19, cex=1.2)
legend(locator(1),legend= c("1","2","3","4","5","6","7","8","9","10"), 
       title = "Month",ncol = 2,
       col=c("firebrick2","darkorange1","gray34","goldenrod2","green3",
                         "cadetblue2","dodgerblue2","mediumpurple2","hotpink","tan"), 
                         pch=19, cex=1.2)
title(main="Years 1 - 3 (2019 - 2021)")
title(main="Year 1 - 2019")
title(main="Year 2 - 2020")
title(main="Year 3 - 2021")
#Station
#Adding the points 
points(cca.p,"sites", pch=19, col= "#A6CEE3", select = met3$Station == "CLV10A")
points(cca.p,"sites", pch=19, col= "#579CC7", select = met3$Station == "KISSR0.0")
points(cca.p,"sites", pch=19, col= "#3688AD", select = met3$Station == "L001")
points(cca.p,"sites", pch=19, col= "#8BC395", select = met3$Station == "L004")
points(cca.p,"sites", pch=19, col= "#89CB6C", select = met3$Station == "L005")
points(cca.p,"sites", pch=19, col= "#40A635", select = met3$Station == "L006")
points(cca.p,"sites", pch=19, col= "#919D5F", select = met3$Station == "L007")
points(cca.p,"sites", pch=19, col= "#F99392", select = met3$Station == "L008")
points(cca.p,"sites", pch=19, col= "#EB444A", select = met3$Station == "LZ2")
points(cca.p,"sites", pch=19, col= "red", select = met3$Station == "LZ25A")
points(cca.p,"sites", pch=19, col= "#F79C5D", select = met3$Station == "LZ30")
points(cca.p,"sites", pch=19, col= "#FDA746", select = met3$Station == "LZ40")
points(cca.p,"sites", pch=19, col= "#FE8205", select = met3$Station == "PALMOUT")
points(cca.p,"sites", pch=19, col= "#E39970", select = met3$Station == "PELBAY3")
points(cca.p,"sites", pch=19, col= "#BFA5CF", select = met3$Station == "POLE3S")
points(cca.p,"sites", pch=19, col= "#8861AC", select = met3$Station == "POLESOUT")
points(cca.p,"sites", pch=19, col= "violet", select = met3$Station == "RITTAE2")
points(cca.p,"sites", pch=19, col= "#E7E099", select = met3$Station == "S308")
points(cca.p,"sites", pch=19, col= "#DEB969", select = met3$Station == "S77")
points(cca.p,"sites", pch=19, col= "#B15928", select = met3$Station == "S79")
#Plotting envfit vectors
plot(ef.cca, col = "black", p.max=0.05)
#Add legend (click to place legend on the outside of the plot) & Title
legend(locator(1),legend= c("CLV10A","KISSR0.0","L001","L004","L005","L006","L007",
                            "L008","LZ2","LZ25A","LZ30","LZ40","PALMOUT","PELBAY3",
                            "POLE3S","POLESOUT","RITTAE2","S308","S77","S79"),
       title = "Station",ncol=2,
       col=c("#A6CEE3","#579CC7","#3688AD","#8BC395","#89CB6C","#40A635","#919D5F",
                      "#F99392","#EB444A","red","#F79C5D","#FDA746","#FE8205","#E39970",
                      "#BFA5CF","#8861AC","violet","#E7E099","#DEB969","#B15928"), 
                      pch=19, cex=0.9)
title(main="Years 1 - 3 (2019 - 2021)")
title(main="Year 1 - 2019")
title(main="Year 2 - 2020")
title(main="Year 3 - 2021")


###### Differential Abundance Analysis - DESEQ2 ######
## USING DESEQ2 (following lashlock github tutorial)
library(DESeq2)

##Differences between years
#load in data WITHOUT rownames
years <- read.csv("feature_Y123_0.01per.csv")
met <- read.csv("Metadata-Diversity_BATCH.csv")
#turning Year into a factor (since it may be read as a number)
met$Year <- as.factor(met$Year)

##Constructing Deseq2 object from data frame
dds <- DESeqDataSetFromMatrix(countData=years, 
                              colData=met, 
                              design=~Year, tidy = TRUE)
#Design specifies how the counts from each gene depend on our variables in the metadata
#For this dataset the factor we care about is the Zone
#tidy=TRUE argument = tells DESeq2 to output the results table with row names as a first #column called 'row.


#let's see what this object looks like
dds
# class: DESeqDataSet 
# dim: 8340 541 
# metadata(1): version
# assays(1): counts
# rownames(8340): 0885965c051f3034c0e28043193bc5d2 51e00e866016fba8a19581249b811ec4 ...
# f9fe4768ad3ef514b97950516e4af5b2 fe2896a859ec05fd0b600b2f633a3bc7
# rowData names(0):
#   colnames(541): KISSR0.0_3_19 L001_3_19 ... S77_10_21 S79_10_21
# colData names(43): Sample Month ... J inv.D


##Running the DESeq function
dds <- DESeq(dds)
#Error in estimateSizeFactorsForMatrix(counts(object),locfunc = 
#locfunc,: every gene contains at least one zero, cannot compute log geometric
#means -> got this error so going to add a pseudocount of 1 to eliminate zeroes
#         (may add bias to the data according to vegan HELP)

##Adding pseudocount of 1 to feature table
#looking at the structure of the data frame
str(years)
#first column is a character so don't include in the transformation

#Adding 1 excluding the first column (ASV column)
years[-1] <- years[-1] + 1


##Retrying the constructing DESeq object and running the DESeq function
dds <- DESeqDataSetFromMatrix(countData=years, 
                              colData=met, 
                              design=~Year, tidy = TRUE)
dds <- DESeq(dds)

##What just happen?
#estimateSizeFactors
#This calculates the relative library depth of each sample 

#estimateDispersions
#estimates the dispersion of counts for each gene 

#nbinomWaldTest
#calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs


##Looking at the results table
res31 <- results(dds)
res31 #looking at the results table
# log2 fold change (MLE): Year 3 vs 1 
# Wald test p-value: Year 3 vs 1 
# DataFrame with 8340 rows and 6 columns
#                                   baseMean log2FoldChange     lfcSE      stat     pvalue      padj
#                                    <numeric>      <numeric> <numeric> <numeric>  <numeric> <numeric>
# 0885965c051f3034c0e28043193bc5d2   1.17377      0.2149733  0.152042  1.413911  0.1573881        NA
# 51e00e866016fba8a19581249b811ec4   1.14815      0.0699605  0.157170  0.445127  0.6562281        NA
# dfd3874c0e70ae177e8cdc4fb6961e7d   1.22257      0.0762662  0.152455  0.500255  0.6168956 0.7984313
# ac879ef0bc703ee2637bc55f0ef97afc   1.24454      0.3709705  0.155215  2.390037  0.0168467 0.0805026
# 41714fa1a258e8098d51d03a1e1b3304   1.20327     -0.3581498  0.149734 -2.391912  0.0167608 0.0802964


##NOTE: If there are more than 2 levels for the variable – as is the case
##for Year w/ 3 levels – results will extract the results table for a comparison 
##of the last level over the first level (so year 3 vs year 1)


##Other comparisons
res23 <- results(dds, contrast = c("Year", "3", "2") )
res23
# log2 fold change (MLE): Year 3 vs 2 
# Wald test p-value: Year 3 vs 2 
# DataFrame with 8340 rows and 6 columns
#                                   baseMean log2FoldChange     lfcSE      stat      pvalue        padj
#                                   <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# 0885965c051f3034c0e28043193bc5d2   1.17377      0.2169855  0.140908  1.539912   0.1235819          NA
# 51e00e866016fba8a19581249b811ec4   1.14815     -0.0866411  0.142640 -0.607409   0.5435795          NA
# dfd3874c0e70ae177e8cdc4fb6961e7d   1.22257     -0.0207447  0.139635 -0.148564   0.8818978   0.9376288
# ac879ef0bc703ee2637bc55f0ef97afc   1.24454      0.3451806  0.142798  2.417259   0.0156379   0.0690231
# 41714fa1a258e8098d51d03a1e1b3304   1.20327     -0.0463689  0.146971 -0.315496   0.7523851          NA

res12 <- results(dds, contrast = c("Year", "1", "2") )
res12 
# log2 fold change (MLE): Year 1 vs 2 
# Wald test p-value: Year 1 vs 2 
# DataFrame with 8340 rows and 6 columns
#                                   baseMean log2FoldChange     lfcSE       stat    pvalue      padj
#                                   <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
# 0885965c051f3034c0e28043193bc5d2   1.17377     0.00201214  0.150700  0.0133519 0.9893470        NA
# 51e00e866016fba8a19581249b811ec4   1.14815    -0.15660154  0.149012 -1.0509325 0.2932896        NA
# dfd3874c0e70ae177e8cdc4fb6961e7d   1.22257    -0.09701090  0.145892 -0.6649517 0.5060814  0.726901
# ac879ef0bc703ee2637bc55f0ef97afc   1.24454    -0.02578994  0.155971 -0.1653507 0.8686680  0.941353
# 41714fa1a258e8098d51d03a1e1b3304   1.20327     0.31178087  0.141655  2.2009921 0.0277366  0.116155

##Saving all comparisons as CSVs
write.csv(res31, "DESEQ-Y13_results.csv")
write.csv(res23, "DESEQ-Y23_results.csv")
write.csv(res12, "DESEQ-Y12_results.csv")



#Visualizing using Volcano plots
##Volcano Plot
par(mfrow=c(1,3))
#Year 3 vs Year 1
# Make a basic volcano plot
with(res31, plot(log2FoldChange, -log10(pvalue), pch=20, main="Year 3 vs. Year 1", xlim=c(-2,2)))
# Add colored points: red = padj<0.05 AND log2FC >1, black = pdj>0.05
with(subset(res31, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res31, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Year 3 vs Year 2
# Make a basic volcano plot
with(res23, plot(log2FoldChange, -log10(pvalue), pch=20, main="Year 3 vs. Year 2", xlim=c(-3,3)))
# Add colored points: red = padj<0.05 AND log2FC >1, black = pdj>0.05
with(subset(res23, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res23, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Year 1 vs Year 2
# Make a basic volcano plot
with(res12, plot(log2FoldChange, -log10(pvalue), pch=20, main="Year 1 vs. Year 2", xlim=c(-3,3)))
# Add colored points: red = padj<0.05 AND log2FC >1, black = pdj>0.05
with(subset(res12, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res12, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

##PCA
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
par(mfrow=c(1,1))
vsdata <- vst(dds, blind=FALSE) #using the DESEQ2 plotPCA function we can 
#look at how our samples group by treatment
plotPCA(vsdata, intgroup="Year")+
  labs(title = "Years 1-3 (2019-2021)")+ 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) 



#### Differences in Zone for EACH YEAR
#loading in data
Y1 <- dat.01per[grep("_19$", rownames(dat.01per)),]
Y2 <- dat.01per[grep("_20$", rownames(dat.01per)),]
Y3 <- dat.01per[grep("_21$", rownames(dat.01per)),]
write.csv(t(Y1), "feature_Y1_0.01per.csv")
write.csv(t(Y2), "feature_Y2_0.01per.csv")
write.csv(t(Y3), "feature_Y3_0.01per.csv")


##Differences found in Zone of Year 1
Y1 <- read.csv("feature_Y1_0.01per.csv")
met1 <- read.csv("Metadata_BATCH_Y1.csv")

##Adding pseudocount of 1
Y1[-1] <- Y1[-1] + 1

##Constructing Deseq2 object and running DESeq function
dds <- DESeqDataSetFromMatrix(countData=Y1, 
                              colData=met1, 
                              design=~Zone, tidy = TRUE)
dds <- DESeq(dds)

##Retrieving results tables for each comparison
resIP <- results(dds, contrast = c("Zone", "Inflow", "Pelagic") )
resIN <- results(dds, contrast = c("Zone", "Inflow", "Nearshore") )
resNP <- results(dds, contrast = c("Zone", "Nearshore", "Pelagic") )
resNS <- results(dds, contrast = c("Zone", "Nearshore", "S79") )
resPS <- results(dds, contrast = c("Zone", "Pelagic", "S79") )
resSI <- results(dds)

##Saving all comparisons as CSVs
write.csv(resIP, "DESEQ-Y1IP_results.csv")
write.csv(resIN, "DESEQ-Y1IN_results.csv")
write.csv(resNP, "DESEQ-Y1NP_results.csv")
write.csv(resNS, "DESEQ-Y1NS_results.csv")
write.csv(resPS, "DESEQ-Y1PS_results.csv")
write.csv(resSI, "DESEQ-Y1SI_results.csv")

##Volcano Plots 
par(mfrow=c(2,3))
#Inflow vs Pelagic
with(resIP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Pelagic", xlim=c(-6,6)))
with(subset(resIP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Inflow vs Nearshore
with(resIN, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Nearshore", xlim=c(-6,6)))
with(subset(resIN, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIN, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. Pelagic
with(resNP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. Pelagic", xlim=c(-4,4)))
with(subset(resNP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. S79
with(resNS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. S79", xlim=c(-8,8)))
with(subset(resNS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Pelagic vs. S79
with(resPS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Pelagic vs. S79", xlim=c(-8,8)))
with(subset(resPS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resPS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#S79 vs Inflow
with(resSI, plot(log2FoldChange, -log10(pvalue), pch=20, main="S79 vs. Inflow", xlim=c(-7,7)))
with(subset(resSI, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resSI, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

##PCA
par(mfrow=c(1,1))
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Zone")+
  labs(title = "Year 1 - Ecological zones")+ 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) 

##Year 2 Zone (No significant differences found but doing it anyway)
Y2 <- read.csv("feature_Y2_0.01per.csv")
met2 <- read.csv("Metadata_BATCH_Y2.csv")

##Adding pseudocount of 1
Y2[-1] <- Y2[-1] + 1

##Constructing Deseq2 object and running DESeq function
dds <- DESeqDataSetFromMatrix(countData=Y2, 
                              colData=met2, 
                              design=~Zone, tidy = TRUE)
dds <- DESeq(dds)

##Retrieving results tables for each comparison
resIP <- results(dds, contrast = c("Zone", "Inflow", "Pelagic") )
resIN <- results(dds, contrast = c("Zone", "Inflow", "Nearshore") )
resNP <- results(dds, contrast = c("Zone", "Nearshore", "Pelagic") )
resNS <- results(dds, contrast = c("Zone", "Nearshore", "S79") )
resPS <- results(dds, contrast = c("Zone", "Pelagic", "S79") )
resSI <- results(dds)

##Saving all comparisons as CSVs
write.csv(resIP, "DESEQ-Y2IP_results.csv")
write.csv(resIN, "DESEQ-Y2IN_results.csv")
write.csv(resNP, "DESEQ-Y2NP_results.csv")
write.csv(resNS, "DESEQ-Y2NS_results.csv")
write.csv(resPS, "DESEQ-Y2PS_results.csv")
write.csv(resSI, "DESEQ-Y2SI_results.csv")

##Volcano Plots 
par(mfrow=c(2,3))
#Inflow vs Pelagic
with(resIP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Pelagic", xlim=c(-5,5)))
with(subset(resIP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Inflow vs Nearshore
with(resIN, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Nearshore", xlim=c(-6,6)))
with(subset(resIN, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIN, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. Pelagic
with(resNP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. Pelagic", xlim=c(-6,6)))
with(subset(resNP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. S79
with(resNS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. S79", xlim=c(-7,7)))
with(subset(resNS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Pelagic vs. S79
with(resPS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Pelagic vs. S79", xlim=c(-7,7)))
with(subset(resPS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resPS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#S79 vs Inflow
with(resSI, plot(log2FoldChange, -log10(pvalue), pch=20, main="S79 vs. Inflow", xlim=c(-7,7)))
with(subset(resSI, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resSI, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

##PCA
par(mfrow=c(1,1))
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Zone")+
  labs(title = "Year 2 - Ecological zones")+ 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) 

##Differences found in Zone of Year 3
Y3 <- read.csv("feature_Y3_0.01per.csv")
met3 <- read.csv("Metadata_BATCH_Y3.csv")

##Adding pseudocount of 1
Y3[-1] <- Y3[-1] + 1

##Constructing Deseq2 object and running DESeq function
dds <- DESeqDataSetFromMatrix(countData=Y3, 
                              colData=met3, 
                              design=~Zone, tidy = TRUE)
dds <- DESeq(dds)

##Retrieving results tables for each comparison
resIP <- results(dds, contrast = c("Zone", "Inflow", "Pelagic") )
resIN <- results(dds, contrast = c("Zone", "Inflow", "Nearshore") )
resNP <- results(dds, contrast = c("Zone", "Nearshore", "Pelagic") )
resNS <- results(dds, contrast = c("Zone", "Nearshore", "S79") )
resPS <- results(dds, contrast = c("Zone", "Pelagic", "S79") )
resSI <- results(dds)

##Saving all comparisons as CSVs
write.csv(resIP, "DESEQ-Y3IP_results.csv")
write.csv(resIN, "DESEQ-Y3IN_results.csv")
write.csv(resNP, "DESEQ-Y3NP_results.csv")
write.csv(resNS, "DESEQ-Y3NS_results.csv")
write.csv(resPS, "DESEQ-Y3PS_results.csv")
write.csv(resSI, "DESEQ-Y3SI_results.csv")

##Volcano Plots 
par(mfrow=c(2,3))
#Inflow vs Pelagic
with(resIP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Pelagic", xlim=c(-5,5)))
with(subset(resIP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Inflow vs Nearshore
with(resIN, plot(log2FoldChange, -log10(pvalue), pch=20, main="Inflow vs. Nearshore", xlim=c(-6,6)))
with(subset(resIN, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resIN, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. Pelagic
with(resNP, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. Pelagic", xlim=c(-7,7)))
with(subset(resNP, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNP, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Nearshore vs. S79
with(resNS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Nearshore vs. S79", xlim=c(-6,6)))
with(subset(resNS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resNS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#Pelagic vs. S79
with(resPS, plot(log2FoldChange, -log10(pvalue), pch=20, main="Pelagic vs. S79", xlim=c(-7,7)))
with(subset(resPS, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resPS, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

#S79 vs Inflow
with(resSI, plot(log2FoldChange, -log10(pvalue), pch=20, main="S79 vs. Inflow", xlim=c(-7,7)))
with(subset(resSI, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resSI, padj>.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))

##PCA
par(mfrow=c(1,1))
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Zone")+
  labs(title = "Year 3 - Ecological zones")+ 
  theme(plot.title.position = "panel")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5)) 


###### Species Cooccurrence (Correlations) ######
library(Hmisc) 

#All Years
x<-read.csv("feature_Y123_0.01per.csv", header=TRUE, row.names=1)
x<-t(x)
y<-rcorr(as.matrix(x, type = c("pearson")))  ## or spearman (pearson may be best here)
yR<-y$r 
yP<-y$P 

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    corr =(cormat)[ut],
    p = pmat[ut]
  )
}
corr_data<-flattenCorrMatrix(y$r, y$P)  


#Note:
# Sort in R or in excel... may want to only keep significant correlations that are 
# to Microcystis specifically to keep it simple. then retain R2 values that are the 
# highest (>0.9 or <-0.9 -- you can change that if you want.) <- cut off will have 
# to be 0.3 since that's the highest 

# Use these values to create network in Cytoscape to visualize the correlations of taxa 
# to Microcystis. 

#Excluding any non-significant correlations (including zeros) and exporting
corr_data <- corr_data[order(corr_data$p),] #sort from smallest to largest
corr_sig <- corr_data[corr_data$p < 0.05, ] #Subsetting data to ONLY include significant correlations
write.csv(corr_sig, "LakeOCorrelationsSigONLY.csv")

#Created network in Cytoscape, merging nodes with taxonomy
node <- read.csv("LakeOCorrelations_Nodes.csv")
tax <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE)
merged <- merge(node,tax, by="FeatureID")
write.csv(merged,"LakeOCorrelations_NodeTaxa.csv")
#Microcystis with corr = 0.7 and up, merging with taxonomy
node <- read.csv("Microcystis Network-0.7+_Node.csv")
tax <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE)
merged <- merge(node,tax, by="FeatureID")
write.csv(merged,"LakeOCorrelations_Microcystis0.7NodeTaxa.csv")


#### ####


