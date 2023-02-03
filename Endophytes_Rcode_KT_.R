#####################################################################################
#Diversity calculation CALCULATION Honours project Kyl Tan 2022/23 
#####################################################################################
#############################################################
# Clean-up the memory and start a new session
#############################################################
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages 
library("ggplot2")
library("vegan")
library("phyloseq") 
library("tidyverse") 
library("RColorBrewer") 
library("ggsci")

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of the code
sessionInfo()

#########################
##Phyloseq
########################

##Manually created OTU table

OTU_table<- read.delim("OTU_table_kt.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#design file
design <- read.delim("map_kt.txt", sep = "\t", header=TRUE, row.names=1)

design <- design[1:12,]
design

##generarte the count file
dat_count <- OTU_table[ ,1:12]
rownames(dat_count) <- rownames(OTU_table)
dat_count[1:5, ]
OTU_table[1:5, ]

#The taxonomy information
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and a new taxa table

Fungal_dat_tax <- read.delim ("Taxonomy_kt.txt", sep = "\t", header=T,row.names=1, blank.lines.skip = FALSE)
Fungal_taxa <- tax_table(as.matrix(Fungal_dat_tax))
dim(Fungal_taxa)


####################################################################################################
#Construction of the phyloseq object

####################################################################################################
#We need design, the data_count and the data_taxonomy file) Optional: data_tree (phylogenentic tree)

#The phyloseq package recapitulates 
#a) The OTU Table counts. This is a n rows x m columns dataframe, where n represents the individual isolates, with their associated abundance counts based on morphotypes counts
#b) The taxonomy information. For each of the n OTUs, a representative sequence (based in morphotypes was Sanger seqeunced) 
#c) The mapping file (often reported as a design file or metadata): it says what the samples are, from where they come from and also it includes additional attributes that can be used in the analysis 


#The OTU Table counts

Fungal_OTU <- otu_table(dat_count, taxa_are_rows=TRUE)

#The taxonomy information
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
dim(Fungal_taxa)

#The mapping file 
Fungal_map <- sample_data(design)

#merge the files and create the phyloseq object
Fungal_data_phyloseq <- merge_phyloseq(Fungal_OTU, Fungal_taxa, Fungal_map)
Fungal_data_phyloseq

#inspect the generated data
Fungal_data_phyloseq
sum(colSums(otu_table(Fungal_data_phyloseq)))
dim(dat_count)
sum(colSums(dat_count))
###################################################
##Figure 5: Species Accumulation Curves  
####################################################
#https://search.r-project.org/CRAN/refmans/vegan/html/specaccum.html

##gives the expected number of observed species or distinct classes as a function of sampling effort
#Number of fungal species (y-axis) vs Number of plants (x-axis)


##Data has species as rows and genotypes as columns, but we need the opposite. Thus, transpose the original data

dat_count_t<-t(dat_count)

#Species accumulation method (partial match). Method "collector" adds sites in the order they happen to be in the data, 
#"random" adds sites in random order, "exact" finds the expected (mean) species richness, 
#"coleman" finds the expected richness following Coleman et al. 1982, and 
#"rarefaction" finds the mean when accumulating individuals instead of sites. 
dat_count_GP<- dat_count[c(1:6)]
dat_count_GP_t<-t(dat_count_GP)
dat_count_RGT<- dat_count[c(7:12)]
dat_count_RGT_t<-t(dat_count_RGT)


rare_curv_RGT<- specaccum(dat_count_RGT_t)
mod1 <- fitspecaccum(rare_curv_RGT, "lomolino")
coef(mod1)
fitted(mod1)
plot(rare_curv_RGT)

specpool(dat_count_RGT_t)

rare_curv_GP<- specaccum(dat_count_GP_t)
mod2 <- fitspecaccum(rare_curv_GP, "lomolino")
coef(mod2)
fitted(mod2)
plot(rare_curv_GP)
specpool(dat_count_GP_t)

Two_genotypes<-c(rare_curv_RGT,rare_curv_GP)

plot(rare_curv_RGT,xlab="Number of Plants", ylab="Number of taxa", col = "#CC79A7",axes = FALSE, lwd = 2)
plot(rare_curv_GP, add = TRUE, col = "#009E73", axes= FALSE, lwd = 2)
plot(rare_curv_124_17, add = TRUE, col = "black", axes= FALSE, lwd = 2)
axis(side=1, at=c(1:8))
axis(side=2, at=c(1:20))
legend(1,16.5, legend= c("RGT Planet", "Golden Promise"),
       col = c("#CC79A7","#009E73"), lty = 1, lwd = 2, cex = 0.8)
box()

####################################
##Figure 6: Beta_diversity  
####################################

#constrained ordinatiton: consatrained for Description

Fungal_endosphere_CAP <- ordinate(Fungal_data_phyloseq, "CAP", "bray", ~ genotype)
plot_ordination(Fungal_data_phyloseq, Fungal_endosphere_CAP, color = "genotype")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(Fungal_data_phyloseq, Fungal_endosphere_CAP , color = "genotype")
p = p + geom_point(size = 6, alpha = 0.80, stroke =1)
p = p + scale_colour_manual(values = c("#009E73", "#CC79A7"))
p + ggtitle("CAP Fungal endhophytic data, Bray distance samples")
p + theme(panel.border = element_rect(color = "black",fill = NA,  size= 1))
#ANOVA on the axis
anova(Fungal_endosphere_CAP, permutations=5000)

#BC distance
BC <- phyloseq::distance(Fungal_data_phyloseq, "bray")
ad<-adonis2(BC ~ genotype, data= design , permutations = 5000)
ad
#################################
##Figure 7:Alpha diversity 
##################################

FE_alpha <-  estimate_richness(Fungal_data_phyloseq, measures = c("Observed", "Shannon", "Chao1")) 
FE_alpha 

FE_otu_table<-otu_table(Fungal_data_phyloseq)
#Data frame Genotype_Description 
colData = design[colnames(FE_otu_table), ]
rownames(colData)
colData


#Description 
design_genotype  <- as.data.frame(colData[, 1]) 
rownames(design_genotype) <- rownames(colData) 
colnames(design_genotype) <- c("genotype") 
design_genotype  

###########################
####Figure 7:Alpha diversity CHAO1 
##########################
#Chao1 ASVs 
FE_alpha_Chao1 <- as.data.frame(FE_alpha[ ,2]) 
rownames(FE_alpha_Chao1) <- rownames(FE_alpha) 
colnames(FE_alpha_Chao1) <- c("Chao1") 


#Combine the dataset sample description and Chao1 OTUs 
FE_alpha_Chao1_TD <- cbind(design_genotype , FE_alpha_Chao1) 
FE_alpha_Chao1_TD <- as.data.frame(FE_alpha_Chao1_TD) 
FE_alpha_Chao1_TD$genotype


#Order the levels according to a defined order 
FE_alpha_Chao1_TD$genotype <- ordered(FE_alpha_Chao1_TD$genotype, levels=c("Golden_Promise","RGT_planet"))  

#Plotting
with(FE_alpha_Chao1_TD, boxplot(Chao1~genotype, xlab = "Genotypes", ylab = "Number of Taxa", col=c("#009E73", "#CC79A7")))
with(FE_alpha_Chao1_TD, stripchart(Chao1~genotype, xlab = "Genotypes", ylab = "Number of Taxa", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))

#ANOVA 
Chao1_OTUs_stats <- aov(Chao1 ~ genotype, data = FE_alpha_Chao1_TD) 
summary(Chao1_OTUs_stats) 


#ttest
Chao1_OTUs_ttest <- t.test(Chao1 ~ genotype, data = FE_alpha_Chao1_TD) 
summary(Chao1_OTUs_ttest) 

Chao1_OTUs_ttest

##############################################################################
#Supplementary Fig 1: Top ten genus 
#############################################################################
# ## agglomerate at the Family taxonomic rank, this Is a new phyloseq object
FE_phylum <- (x1 <- tax_glom(Fungal_data_phyloseq, taxrank="Phylum") )
# ## How many taxa before/after agglomeration?
ntaxa(FE_phylum); ntaxa(x1)

#agglomerate at family level
FE_family <-  tax_glom(Fungal_data_phyloseq, taxrank="Family")
FE_family

#agglomerate at family level
FE_Genus <-  tax_glom(Fungal_data_phyloseq, taxrank="Genus")
FE_Genus

tax_table(FE_family)[1:11, ]

#transform in to relative abundance
ps_phylum_0 <- transform_sample_counts(FE_Genus, function(x) x / sum(x))
#abundance of all samples plot
plot_bar(ps_phylum_0, fill="Genus")
#merge samples by treatment
ps_phylum_1 <- merge_samples(ps_phylum_0, "genotype")
#transform to relative abudance
ps_phylum_2 <- transform_sample_counts(ps_phylum_1, function(x) x / sum(x))
#plot ra of all phyla
plot_bar(ps_phylum_2, fill="Genus")

sample_data(ps_phylum_2)
tax_table(ps_phylum_2)

#plotting based on sample type
df_phylum <- psmelt(ps_phylum_2)
#write.table(df_phylum, file="df_phylum.txt ", sep="\t")


top_phylum <- df_phylum %>%
  group_by("genotype", Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phylum

#Plotting based on order Supplementary Figure 2
top_phylum <- top_phylum$Phylum[1:5]
df_phylum_0 <- df_phylum %>%
  mutate(Phylum = fct_other(Phylum, top_phylum))
plot_top__phylum <-ggplot(df_phylum_0, aes(Sample, Abundance, fill = fct_reorder(Order, Abundance))) + geom_col() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
plot_top__phylum

#now to add colour blind friendly pallet

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot_top__phylum_colour <- plot_top__phylum + scale_fill_manual(values=cbPalette)

plot_top__phylum_colour

#############################
#### Supplementary Figure 2: Alpha diversity SHANNON 
#############################
#Shannon ASvs 
FE_alpha_Shannon <- as.data.frame(FE_alpha[ ,4]) 
rownames(FE_alpha_Shannon) <- rownames(FE_alpha) 
colnames(FE_alpha_Shannon) <- c("Shannon") 


#Combine the dataset sample description and Shannon OTUs 
FE_alpha_Shannon_TD <- cbind(design_genotype , FE_alpha_Shannon) 
FE_alpha_Shannon_TD <- as.data.frame(FE_alpha_Shannon_TD) 
FE_alpha_Shannon_TD$genotype

#Order the levels according to a defined order 
FE_alpha_Shannon_TD$genotype <- ordered(FE_alpha_Shannon_TD$genotype, levels=c("Golden_Promise","RGT_planet"))  

#Plotting
with(FE_alpha_Shannon_TD, boxplot(Shannon ~ genotype, xlab = "Genotypes", ylab = "Number of Taxa",col=c("#009E73", "#CC79A7")))
with(FE_alpha_Shannon_TD, stripchart(Shannon~genotype, xlab = "Genotypes", ylab = "Number of Taxa", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Shannon_OTUs_stats <- aov(Shannon ~ genotype, data = FE_alpha_Shannon_TD) 
summary(Shannon_OTUs_stats) 

#ttest
Shannon_OTUs_ttest <- t.test(Shannon ~ genotype, data = FE_alpha_Shannon_TD) 
summary(Shannon_OTUs_ttest) 

Shannon_OTUs_ttest
##End
