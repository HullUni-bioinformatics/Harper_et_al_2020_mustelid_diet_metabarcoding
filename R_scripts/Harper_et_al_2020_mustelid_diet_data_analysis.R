#' ---
#' Title: "Using DNA metabarcoding to investigate diet and niche partitioning in the native European otter (Lutra lutra) and invasive American mink (Neovison vison)"
#' Author: "Lynsey Rebecca Harper"
#' Date: "20th November 2020"
#' ---
#' 
#' 
#' European otter and American mink faeces collected at three sites across
#' northern and eastern England, were analysed using DNA metabarcoding 
#' with nested primers that amplify the 12S ribosomal RNA region across 
#' all vertebrates. This project was a proof-of-concept to examine the 
#' diet of both species using molecular scatology. DNA metabarcoding
#' provided inventories of all vertebrate species eaten by otter and mink. 
#' We will assess whether diets overlap or whether niche partitioning 
#' has occurred.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. 
#' 

## Clear memory
rm(list=ls())

## set working directory to the location of the script
# install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()

## On first run, install developmental version of ggmap and legendMap
#devtools::install_github("dkahle/ggmap", force=TRUE)
#devtools::install_github("3wen/legendMap")

## Load required packages
p <- c("plyr","tidyverse","ggpubr","scales","ggmap","legendMap","gtools",
       "reshape2","bipartite","car","multcomp","FSA","lawstat","vegan",
       "betapart","iNEXT")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                          repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies=TRUE)
lapply(p, require, character.only = TRUE)

## Load custom functions
f <- c("CheckResidsFunction.R")
lapply(f, source)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()



#' ---
#' 
#' ## 1) Maps
#' 
#' Plot sampling sites for otter and mink dietary metabarcoding.
#' 

## Import sample metadata
metadata <- read.csv("../Data/metadata.csv", header=TRUE)

## Change labels for sampling locations across River Hull catchment to
## River Hull
metadata$Site <- gsub("River Hull - North catchment", "River Hull", metadata$Site)
metadata$Site <- gsub("River Hull - Tophill Low", "River Hull", metadata$Site)
metadata$Site <- as.factor(metadata$Site)

## Check coordinates plot ok
test <- ggplot(metadata, aes(x=Longitude, y=Latitude)) + geom_point()
test

## Provide API key for Google Maps
register_google(key = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
                account_type = "premium", 
                day_limit = 100000)


######
# UK #
######

## Now we want apply coordinates to a map of the UK. First get the map by 
## downloading the Google terrain map for the UK. The code below will 
## fetch the map from Google, centred on a specific latitude and longitude.

## Find mean latitude and longitude to centre UK map
mean(metadata$Latitude)
mean(metadata$Longitude)

## Get UK map from Google
UK <- get_googlemap(center=c(lon=-0.2927766, lat=53.74583), 
                    zoom=5,
                    color="color",
                    style="feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
ggmap(UK)

## Now plot samples on the map
m1 <- ggmap(UK, extent="panel", crop=T) + 
        geom_point(data=metadata, 
                   aes(x=Longitude, y=Latitude, shape=Read_ID), 
                   size=2) + 
        scale_shape_manual(name="Predator",
                           values=c(15,4,16,8,17,18)) + 
        scale_bar(lon=-13, lat=46, 
                  distance_lon=250, distance_lat=50, 
                  distance_legend=100, dist_unit="km",
                  arrow_length=250, arrow_distance=150, 
                  arrow_north_size=5) + 
        labs(title="A  UK", x="Longitude", y="Latitude") + 
        theme_bw() + 
        theme(plot.title = element_text(face="bold"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              legend.position="right")
m1

## Split dataframe by sampling sites
Malham.coords <- droplevels(subset(metadata, Site == "Malham Tarn"))
Hull.coords <- droplevels(subset(metadata, Site == "River Hull"))
Glaven.coords <- droplevels(subset(metadata, Site == "River Glaven"))


###############
# Malham Tarn #
###############

## Find latitudes and longitudes to centre map of Malham Tarn
mean(Malham.coords$Latitude)
mean(Malham.coords$Longitude)

## Get map of Malham Tarn
Malham <- get_googlemap(center=c(lon=-2.157743, lat=54.09705), 
                        zoom=14,
                        color="color",
                        style="feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
ggmap(Malham)

## Now plot sampling sites on the map
m2 <- ggmap(Malham, extent="panel", crop=T) + 
        geom_point(data=Malham.coords, 
                   aes(x=Longitude, y=Latitude), 
                   size=2, pch=17) + 
        scale_bar(lon=-2.183, lat=54.083,
                  distance_lon=0.5, distance_lat=0.08, 
                  distance_legend=0.15, dist_unit="km",
                  arrow_length=0.5, arrow_distance=0.25, 
                  arrow_north_size=5) + 
        labs(title="B  Malham Tarn", x="Longitude", y="Latitude") + 
        theme_bw() + 
        theme(plot.title = element_text(face="bold"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              legend.position="none")
m2


##############
# River Hull #
##############

## Find latitudes and longitudes to centre map of the River Hull
mean(Hull.coords$Latitude)
mean(Hull.coords$Longitude)

## Get map of the River Hull
Hull <- get_googlemap(center=c(lon=-0.3703279, lat=53.92386),
                      zoom=11,
                      color="color",
                      style="feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
ggmap(Hull)

## Now plot sampling sites on the map
m3 <- ggmap(Hull, extent="panel", crop=T) + 
        geom_point(data=Hull.coords, 
                   aes(x=Longitude, y=Latitude, shape=Read_ID),
                   size=2) + 
        scale_shape_manual(name="Predator",
                           values=c(15,4,16,8,17,18)) + 
        scale_bar(lon=-0.57, lat=53.81, 
                  distance_lon=5, distance_lat=0.6, 
                  distance_legend=1.2, dist_unit="km",
                  arrow_length=3, arrow_distance=2, 
                  arrow_north_size=5) + 
        labs(title="C  River Hull", x="Longitude", y="Latitude") + 
        theme_bw() + 
        theme(plot.title = element_text(face="bold"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              legend.position="none")
m3


################
# River Glaven #
################

## Find latitudes and longitudes to centre map of the River Glaven
mean(Glaven.coords$Latitude)
mean(Glaven.coords$Longitude)

## Get map of the River Glaven
Glaven <- get_googlemap(center=c(lon=1.092444, lat=52.9052), 
                        zoom=12,
                        color="color",
                        style="feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")
ggmap(Glaven)

## Now plot sampling sites on the map
m4 <- ggmap(Glaven, extent="panel", crop=T) + 
        geom_point(data=Glaven.coords, 
                   aes(x=Longitude, y=Latitude, shape=Read_ID), 
                   size=2) + 
        scale_shape_manual(name="Predator",
                           values=c(15,4,16,8,17,18)) + 
        scale_bar(lon=0.995, lat=52.847, 
                  distance_lon=2, distance_lat=0.3, 
                  distance_legend=0.6, dist_unit="km",
                  arrow_length=1.5, arrow_distance=1, 
                  arrow_north_size=5) + 
        labs(title="D  River Glaven", x="Longitude", y="Latitude") + 
        theme_bw() + 
        theme(plot.title = element_text(face="bold"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              legend.position="none")
m4


############
# ALL MAPS #
############

g1 <- ggarrange(m1, m2, m3, m4, nrow=2, ncol=2, 
                common.legend=TRUE, legend="bottom")

#ggsave(filename="Figures/FigS1.png", 
#       plot = g1, width = 10, height = 10, dpi = 300, units = "in")



###################
# REFINE METADATA #
###################

## Refine sample metadata for downstream analyses
## Remove coordinates and predator ID upon field collection as these 
## are no longer required
metadata <- metadata[,-c(4:5,9)]

## Rename Read_ID column
colnames(metadata)[7] <- "Predator"

## Remove samples that contained an inadequate number of predator reads
## for predator assignment, or multiple predators with <90% of predator
## reads (see "Appendix2_predator_assignment.xlsx")
metadata <- metadata[which(!grepl("Inadequate|Multiple", metadata$Predator)),]

## Drop all unused factor levels
metadata <- droplevels(metadata)



#' ---
#' 
#' ## 2) Metabarcoding dataset
#' 
#' The raw output from metaBEAT is not particularly suited to applying
#' statistical analysis in R. Files must be manipulated for downstream 
#' analyses.
#' 

## Import metabarcoding data
ass.raw.2018 <- read.csv("../Data/Library1_ReferenceDatabaseBLAST.csv", header=TRUE)
ass.raw.2019 <- read.csv("../Data/Library2_ReferenceDatabaseBLAST.csv", header=TRUE)
unass.raw.2018 <- read.csv("../Data/Library1_unassignedBLAST.csv", header=TRUE)
unass.raw.2019 <- read.csv("../Data/Library2_unassignedBLAST.csv", header=TRUE)

## Make list of dataframes
raw.data <- list(ass.raw.2018, ass.raw.2019, 
                 unass.raw.2018, unass.raw.2019)
names(raw.data) <- c("ass.raw.2018.formatted", "ass.raw.2019.formatted", 
                     "unass.raw.2018.formatted", "unass.raw.2019.formatted")

## Format dataframes for downstream analyses
raw.data <- lapply(raw.data, function (x) {
        names(x) <- lapply(x[1,], as.character)
        x <- x[-1,]
        rownames(x) <- NULL
        colnames(x) <- gsub("-nc.blast.blast", "", colnames(x))
        colnames(x) <- gsub("-nc.blast", "", colnames(x))
        colnames(x)[1] <- "Assignment"
        x$Assignment <- gsub("s__|g__|f__|o__|c__|p__|k__|sk__|__", "", 
                             x$Assignment)
        return(x)
})

## Unlist dataframes
list2env(raw.data, .GlobalEnv )

## Create new dataframe for calculating the proportional read counts 
## for each taxon downstream.
## Duplicate raw read data
raw1 <- ass.raw.2018.formatted
raw2 <- ass.raw.2019.formatted
raw3 <- unass.raw.2018.formatted
raw4 <- unass.raw.2019.formatted

## Remove unassigned from reference database BLAST results as unassigned
## reads were extracted for BLAST against entire NCBI nucleotide database.
## Also remove last column containing taxonomy.
raw1 <- raw1[-100,-127]
raw2 <- raw2[-98,-154]
raw3 <- raw3[,-127]
raw4 <- raw4[,-154]

## Bind data frames
raw.merged <- smartbind(raw1, raw2, raw3, raw4)

## Replace NAs with 0
raw.merged[is.na(raw.merged)] <- 0

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
raw.merged[,2:278] <- lapply(raw.merged[,2:278], function(x) as.numeric(as.character(x)))
raw.merged <- ddply(raw.merged, .(Assignment), numcolwise(sum))

## Make Assignment column row names
rownames(raw.merged) <- raw.merged$Assignment
raw.merged <- raw.merged[,-1]

## Calculate total number of reads in samples/controls
raw.merged <- rbind(raw.merged, colSums(raw.merged))

## Make new dataframe containing sample ID and total number of reads 
## per sample
raw.total <- raw.merged[174,]
raw.total$Assignment <- "Total"
raw.total <- raw.total[,c(278,1:277)]
rownames(raw.total) <- NULL

## Remove any taxonomic assignments that aren't vertebrate from each dataframe
## using the taxonomy column created during processing with metaBEAT
## Make list of dataframes
formatted.data <- list(ass.raw.2018.formatted, ass.raw.2019.formatted,
                       unass.raw.2018.formatted, unass.raw.2019.formatted)
names(formatted.data) <- c("ass.2018", "ass.2019", "unass.2018", "unass.2019")

## Format dataframes for downstream analyses
formatted.data <- lapply(formatted.data, function (x) {
        x <- x[which(grepl("Chordata", x$taxomomy)),]
        x <- x[,which(!grepl("taxomomy", colnames(x)))]
        return(x)
})

## Unlist dataframes
list2env(formatted.data, .GlobalEnv )

## Bind data frames
merged.df <- smartbind(ass.2018, ass.2019, unass.2018, unass.2019)

## Replace NA values with 0
merged.df[is.na(merged.df)] <- 0

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
merged.df[,2:278] <- lapply(merged.df[,2:278], function(x) as.numeric(as.character(x)))
merged.df <- ddply(merged.df, .(Assignment), numcolwise(sum))

## Remove samples that were analysed with a different extraction kit and for
## a different project
merged.df <- merged.df[,-c(88,90,92,94,96,98,100,102)]
merged.df <- merged.df[,which(!grepl("LIB06|LIB07|toad", 
                                     colnames(merged.df)))]

## Remove empty taxonomic assignments that may be specific to samples 
## extracted with different kit and from different project
merged.df <- merged.df[rowSums(merged.df[, -1]) > 0,]

## Reset row names for further indexing
rownames(merged.df) <- NULL

## Export as .csv file
#write.csv(merged.df, "../Data/Otter_metabarcoding_merged_20200616.csv", 
#          row.names=FALSE)



#' --- 
#' 
#' ## 3) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed. Use NBN atlas to check species
#'    occurrence records and ensure they match with sampling locations. 
#'    This is also a good source for checking current taxonomy.
#' 2. Any Genus or Family assignments containing only one species in 
#'    the UK must be changed to that species.
#' 3. Likewise, for any species assignment which is the only species in 
#'    the UK and also has genus/family reads, read counts from all
#'    assignments must be merged.
#'  
#' A record of these changes will be kept as these changes are made.
#'      

## Inspect new dataframe
summary(merged.df)
head(merged.df)
names(merged.df)
str(merged.df)

#'
#' Now, remove spurious assignments, i.e. vertebrates not found in the study 
#' area, and bacteria/invertebrates that should not have been amplified by the 
#' vertebrate primers. Spurious assignments are as follows:
#' 
#' - *Clupea harengus*, row 29
#' - *Pan*, row 80
#' - *Umbra pygmaea*, row 123
#' - *Zenaida macroura*, row 127
#' 

## Remove spurious assignments
true.assign <- merged.df[-c(29,80,123,127),]

## Reset row names of data frame for further indexing
rownames(true.assign) <- NULL

## Remove underscore from species names
true.assign$Assignment <- gsub("_", " ", true.assign$Assignment)

#'
#' Now, change genus/family level assignments where only one UK species 
#' is contained within that genus/family to species level, or unreliable 
#' species assignments to genus/family level. Also, reassign *Canis lupus* 
#' and *Sus scrofa* to domestic subspecies. Assignments to be changed are 
#' as follows:
#' 
#' - *Anas carolinensis* = *Anas*, row 6
#' - *Canis lupus* = *Canis lupus familiaris*, row 19
#' - Cichlidae = *Maylandia zebra*, row 28
#' - Cottidae = *Cottus gobio*, row 34
#' - *Cottus* = *Cottus gobio*, row 35
#' - Felidae = *Felis catus*, row 42
#' - Hominidae = *Homo sapiens*, row 51
#' - *Melagris* = *Meleagris gallopavo*, row 64
#' - Primates = *Homo sapiens*, row 94
#' - *Rhamphochromis esox* = *Maylandia zebra*, row 100
#' - *Sus scrofa* = *Sus scrofa domesticus*, row 112
#' - *Vulpes* = *Vulpes vulpes*, row 122
#' 
#' These species level assignments already exist in the metabarcoding
#' dataset, therefore we can simply merge the genus/family read counts
#' with the species read counts.
#' 

## Rename the genus/family assignments 
true.assign$Assignment <- as.character(true.assign$Assignment)
true.assign[6, "Assignment"] <- "Anas"
true.assign[19, "Assignment"] <- "Canis lupus familiaris"
true.assign[c(28,100), "Assignment"] <- "Maylandia zebra"
true.assign[34:35, "Assignment"] <- "Cottus gobio"
true.assign[42, "Assignment"] <- "Felis catus"
true.assign[c(51,94), "Assignment"] <- "Homo sapiens"
true.assign[64, "Assignment"] <- "Meleagris gallopavo"
true.assign[112, "Assignment"] <- "Sus scrofa domesticus"
true.assign[122, "Assignment"] <- "Vulpes vulpes"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
true.assign <- ddply(true.assign, .(Assignment), numcolwise(sum))

## Create separate dataframes for fecal samples and controls
samples <- true.assign[,which(!grepl("Negative|Positive|pos|neg",
                                     colnames(true.assign)))]

controls <- true.assign[,which(grepl("Assignment|Negative|Positive|pos|neg", 
                                     colnames(true.assign)))]



#' --- 
#' 
#' ## 4) Contamination
#' 
#' Examine how much contamination occured in PCR positive and negative 
#' controls.
#' 

## Duplicate control data
contamination <- controls

## Make Assignment column row names
rownames(contamination) <- contamination$Assignment
contamination <- contamination[,-1]

## Subset raw totals for PCR positive and negative controls
to.keep <- as.character(colnames(contamination))
control.total <- raw.total[(colnames(raw.total) %in% to.keep)]

## Add total read counts for each control to contamination dataframe
contamination <- rbind(contamination, control.total)

## Calculate proportional read counts for each taxon detected in each
## control
contamination <- contamination/c(contamination[114,])

## Replace NA values with 0
contamination[is.na(contamination)] <- 0

## Remove row containing total proportional read counts from dataframe
contamination <- contamination[-114,]

## Make row names first column in dataframe
contamination <- tibble:::rownames_to_column(contamination, "Assignment")

## Remove empty taxonomic assignments
contamination <- contamination[rowSums(contamination[, -1]) > 0,]

## Reset row names for further indexing
rownames(contamination) <- NULL

## Transpose dataframe
contaminants <- setNames(data.frame(t(contamination[,-1])), contamination[,1])

## Make row names first column in dataframe
contaminants <- tibble:::rownames_to_column(contaminants, "ID")

## Create column specifying type of negative control
contaminants$Type <- ifelse(grepl("Negative|neg", contaminants$ID), "Negative",
                            ifelse(grepl("Positive|pos", contaminants$ID), "Positive",
                                   "Other"))

## Move to start of dataframe
contaminants <- contaminants[,c(1,32,2:31)]

## Melt dataframe for plotting
contaminants <- melt(contaminants, c("ID","Type"))

## Rename columns
colnames(contaminants)[3:4] <- c("Assignment", "PRC")

## Only keep assignments that have reads
contaminants <- filter(contaminants, PRC > 0)

## Plot contamination found in PCR controls
hm0 <- ggplot(contaminants, aes(x=ID, 
                                y=fct_rev(as_factor(Assignment)), 
                                fill=PRC)) + 
        geom_tile(colour="black") + 
        scale_fill_gradient(name="Proportional\nread counts", 
                            limits=c(0,1),
                            breaks=c(0,0.25,0.50,0.75,1),
                            low="white", high="red",
                            guide=guide_colourbar(frame.colour="black",
                                                  ticks.colour="black")) + 
        labs(x="PCR controls", y="Taxonomic assignment") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"), 
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour="black"),
              text = element_text(size=20),
              legend.key.size = unit(1, 'lines')) + 
        facet_grid(. ~ Type, scales="free",space="free")
hm0

#ggsave(filename="Figures/FigS2.png", 
#       plot = hm0, width = 13, height = 15, dpi = 300, units = "in")



#' ---
#' 
#' ## 5) Remove potential false positives
#' 
#' Now the data is in a form that can be manipulated easily, it must be 
#' filtered to remove contaminants and potential false positives. We will
#' calculate the maximum frequency of cichlid DNA contamination across 
#' all faecal samples and apply this value as our false positive sequence
#' threshold. PCR negative controls are ineffective for calculating
#' sequence thresholds as negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially.
#' 

## Subset raw totals for samples
to.keep <- as.character(colnames(samples))
sample.total <- raw.total[(colnames(raw.total) %in% to.keep)]

## Bind total read counts to refined dataset
cichlid <- rbind(samples, sample.total)

## Make Assignment column row names
cichlid <- tibble:::column_to_rownames(cichlid, "Assignment")

## Create new dataframe containing the proportional read counts for each
## taxon detected in each sample
cichlid.freq  <- cichlid/c(cichlid[114,])
cichlid.freq[is.na(cichlid.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
cichlid.freq$Threshold <- apply(cichlid.freq, 1, max)

## Check Threshold column has been created properly
head(cichlid.freq[,215:217])

## Combine read frequencies with taxonomic assignment
cichlid.freq <- tibble:::rownames_to_column(cichlid.freq, "Assignment")

## Print contamination threshold based on max. level of cichlid
## DNA contamination
max(cichlid.freq[57,218])   # 0.01122693

## In this scenario, any assignments <1.12% total reads in faecal 
## samples would be considered a false positive. To examine the effect of 
## this threshold on the data, replace all assignments that are less than 
## or equal to this value for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
cichlid.test <- cichlid.freq
cichlid.test[cichlid.test <= 0.01122693] <- 0

## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples
cichlid.test <- cichlid.test[-114,-218]
rownames(cichlid.test) <- NULL
cichlid.conversion <- rbind(cichlid.test, sample.total)

## Now convert frequencies back to read counts
cichlid.conversion <- tibble:::column_to_rownames(cichlid.conversion, "Assignment")
cichlid.FP <- cichlid.conversion*c(cichlid.conversion[114,])

## Remove total row, reset row names, recreate Assignment column
cichlid.FP <- cichlid.FP[-114,]
cichlid.FP <- tibble:::rownames_to_column(cichlid.FP, "Assignment")

## Check whether cichlid DNA has been removed from eDNA samples
subset(cichlid.FP, cichlid.FP$Assignment=="Maylandia zebra")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- cichlid[-114,]
temp$Assignment <- cichlid.FP$Assignment
temp <- temp[,c(217,1:216)]
temp[,2:217][temp[,2:217] > 0] <- 1
cichlid.FP[,2:217][cichlid.FP[,2:217] > 0] <- 1
cichlid1 <- data.frame(colSums(temp[,2:217]))
cichlid2 <- data.frame(colSums(cichlid.FP[,2:217]))
cichlid.compare <- cbind(cichlid1, cichlid2)

## Calculate proportion of information retained
cichlid.compare$proportion <- cichlid.compare[,2]/cichlid.compare[,1]*100
cichlid.compare[is.na(cichlid.compare)] <- 0

## Examine mean and range of proportion of taxa retained
range(cichlid.compare$proportion)
mean(cichlid.compare$proportion)
## This would result in up to 94% taxa being removed. On average, 39.66% 
## species detections are retained.

## Recreate dataframe with the cichlid false positive threshold applied 
## Convert frequencies back to read counts, remove total row, reset row 
## names, recreate Assignment column
cichlid.FP <- cichlid.conversion*c(cichlid.conversion[114,])
cichlid.FP <- cichlid.FP[-114,]
cichlid.FP <- tibble:::rownames_to_column(cichlid.FP, "Assignment")

## Remove empty taxonomic assignments
cichlid.FP <- cichlid.FP[rowSums(cichlid.FP[, -1]) > 0,]

## Export as .csv file
#write.csv(cichlid.FP, "../Data/Otter_metabarcoding_FP-threshold-applied_20200616.csv",
#          row.names=FALSE)



#' ---
#'
#' 6) Compare data with and without false positive threshold
#' 

## Duplicate dataframes with and without false positive threshold applied
samples.NT <- cichlid.freq[-114,-218]
samples.TA <- cichlid.test

## Add new column stating whether threshold applied or not
samples.NT$Threshold <- "Pre-threshold"
samples.TA$Threshold <- "Post-threshold"

## Melt dataframes for plotting
samples.NT <- melt(samples.NT, c("Assignment","Threshold"))
samples.TA <- melt(samples.TA, c("Assignment","Threshold"))

## Rename columns
colnames(samples.NT)[3:4] <- c("ID", "PRC")
colnames(samples.TA)[3:4] <- c("ID", "PRC")

## Keep only assignments with more than 0 read counts
samples.NT <- filter(samples.NT, PRC > 0)
samples.TA <- filter(samples.TA, PRC > 0)

## Get list of taxa removed by false positive threshold
missing.taxa <- unique(samples.NT$Assignment)[!unique(samples.NT$Assignment) %in% unique(samples.TA$Assignment)]

## Combine dataframes
compare.threshold <- rbind(samples.NT, samples.TA)

## Create new factor specifying order that facets are to be plotted
compare.threshold$fThreshold <- factor(compare.threshold$Threshold,
                                       levels=c("Pre-threshold", 
                                                "Post-threshold"))

## Order dataframe by taxonomic assignment
compare.threshold <- compare.threshold[order(compare.threshold$Assignment),]

## Create custom y axis labels (add asterisk to taxa removed by threshold)
labels <- unique(compare.threshold$Assignment)
labels[labels %in% missing.taxa] <- paste("*", labels[labels %in% missing.taxa])
labels

## Plot taxonomic assignments before threshold application
hm1 <- ggplot(compare.threshold, aes(x=ID, 
                                     y=fct_rev(as_factor(Assignment)), 
                                     fill=PRC)) + 
        geom_tile(colour="black") + 
        scale_fill_gradient(name="Proportional\nread counts    ", 
                            limits=c(0,1), 
                            breaks=c(0,0.25,0.5,0.75,1),
                            low="white", high="blue",
                            guide=guide_colourbar(frame.colour="black",
                                                  ticks.colour="black")) + 
        scale_y_discrete(labels=rev(labels)) + 
        labs(x="Samples", y="Taxonomic assignment") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"), 
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              axis.title.x = element_text(size = 30, margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(size = 30, margin = unit(c(0, 5, 0, 0), "mm")),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size=16, colour="black"),
              strip.text.x = element_text(size=20),
              legend.key.size = unit(3, "lines"),
              legend.key.width = unit(5, "lines"),
              legend.title = element_text(size=20),
              legend.text = element_text(size=20),
              legend.position = "bottom",
              legend.direction = "horizontal") + 
        facet_grid(. ~ fThreshold, scales="free", space="free")
hm1

#ggsave(filename="Figures/FigS3.png", 
#       plot = hm1, width = 25, height = 30, dpi = 300, units = "in")



#' ---
#'
#' ## 7) Spurious assigments
#'
#' Remove vertebrate assignments above species level, excluding the 
#' genera and families which have no species-level assignments or known 
#' issues in the current UK vertebrate reference database. Also remove 
#' the PCR positive control and domestic species.
#' 

## Remove any assignments above species level with exceptions
spp.df <- cichlid.FP[which(grepl(" |Anas|Aythya|Laridae|Microtus|Percidae|Tringa", 
                                 cichlid.FP$Assignment)),]

## Reset row names of data frame for further indexing
rownames(spp.df) <- NULL

## Correct names for genera and families
spp.df$Assignment <- as.character(spp.df$Assignment)
spp.df[3, "Assignment"] <- "Anas spp."
spp.df[6, "Assignment"] <- "Aythya spp."
spp.df[22, "Assignment"] <- "Laridae spp."
spp.df[25, "Assignment"] <- "Microtus spp."
spp.df[32, "Assignment"] <- "Pelophylax ridibundus"
spp.df[34, "Assignment"] <- "Percidae spp."
spp.df[50, "Assignment"] <- "Tringa spp."

## Merge read counts by taxonomic assignment
spp.df <- ddply(spp.df, .(Assignment), numcolwise(sum))

## Remove domestic species from samples
spp.df <- spp.df[-c(8,10,20,46),]

## Reset row names for further indexing
rownames(spp.df) <- NULL



#' ---
#'
#' ## 8) Tidy data for downstream analyses
#'
#' We need to remove samples that have no predator reads (as predator 
#' cannot be confirmed), and samples where there are reads for multiple 
#' predators but no single predator is associated with >90% of the total 
#' predator reads. These samples were determined in an excel spreadsheet 
#' ("Appendix2_predator_assignment.xlsx):
#' 
#' - LIB01-PS-Box01-T05
#' - LIB01-PS-Box01-T11
#' - LIB01-PS-Box01-T12
#' - LIB02-TL02
#' - LIB02-TL04
#' - LIB02-TL08
#' - LIB02-TL11
#' - LIB02-TL15
#' - LIB02-TL17
#' - LIB04-TL42
#' - O3
#' - TL104
#' - TL113
#' - TL115
#' - TL117
#' - TL120
#' - TL122
#' - TL126
#' - TL136
#' - TL93 
#'

## Remove samples with problematic predator ID
filtered.df <- spp.df[,which(!grepl("LIB01-PS-Box01-T05|LIB01-PS-Box01-T11|LIB01-PS-Box01-T12|LIB02-TL02|LIB02-TL04|LIB02-TL08|LIB02-TL11|LIB02-TL15|LIB02-TL17|LIB04-TL42|O3|TL104|TL113|TL115|TL117|TL120|TL122|TL126|TL136|TL93",
                                   colnames(spp.df)))]

## Create new dataframe containing only sample ID and associated predator
predators <- data.frame(metadata$Sequencing_ID, metadata$Predator)
colnames(predators) <- c("Sample","Predator")

## Create summary dataframe
summary.stats <- setNames(data.frame(t(filtered.df[,-1])), filtered.df[,1])
summary.stats <- tibble::rownames_to_column(summary.stats, "Sample")
summary.stats <- merge(predators, summary.stats, by="Sample")
summary.stats <- summary.stats[,-1]
summary.stats <- ddply(summary.stats, .(Predator), numcolwise(sum))
summary.stats[,2:47][summary.stats[,2:47] > 0] <- 1
summary.stats <- summary.stats[,which(!grepl("Lutra|Mustela|Neovison|Vulpes",
                                             colnames(summary.stats)))]

## Calculate total number of prey taxa for each predator
summary.stats$Total <- rowSums(summary.stats[,2:43])

## Now remove samples associated with fox and polecat as sample sizes 
## are too small to be included for downstream analyses
## Create vector of sample IDs
to.remove <- droplevels(subset(metadata, select="Sequencing_ID", 
                               Predator %in% c("Fox","Polecat")))
final.df <- filtered.df[!(colnames(filtered.df) %in% to.remove$Sequencing_ID)]
metadata <- droplevels(subset(metadata, !Predator %in% c("Fox","Polecat")))

## Remove fox and polecat from final dataframe as these would now be 
## considered contamination of otter and mink samples
final.df <- final.df[which(!grepl("Mustela|Vulpes", final.df$Assignment)),]

## Export as .csv file
#write.csv(final.df, "../Data/Final_dataset_for_analyses_20200616.csv", 
#          row.names=FALSE)



#' ---
#'
#' ## 9) Basic summaries
#'
#' Summarise otter and mink in terms of proportional read counts and 
#' detection rate across samples. Plot as pie chart, bar plot, and heat 
#' map. For pie charts, we will produce one chart including predator 
#' reads and a chart with prey reads only.
#' 

###############################################################
# Pie chart: read counts by vertebrate group (inc. predators) #
###############################################################

## Make row names of final dataframe the same as Assignment column
rownames(final.df) <- final.df$Assignment

## Create vectors containing taxa belonging to different vertebrate groups
amphibian <- c("Bufo bufo", "Pelophylax ridibundus", "Rana temporaria")
fish <- c("Abramis brama","Anguilla anguilla","Barbatula barbatula",
          "Carassius carassius","Cottus gobio","Esox lucius",
          "Gasterosteus aculeatus","Gobio gobio","Lampetra fluviatilis",
          "Oncorhynchus mykiss","Percidae spp.","Phoxinus phoxinus",
          "Platichthys flesus","Pungitius pungitius","Rutilus rutilus",
          "Salmo trutta","Scardinius erythrophthalmus","Thymallus thymallus",
          "Tinca tinca")
bird <- c("Alectoris rufa","Anas spp.","Aythya spp.","Columba oenas",
          "Fulica atra","Gallinula chloropus","Garrulus glandarius",
          "Laridae spp.","Phalacrocorax carbo","Phasianus colchicus",
          "Scolopax rusticola","Sturnus vulgaris","Tringa spp.")
mammal <- c("Arvicola amphibius","Lepus europaeus","Microtus spp.",
            "Myodes glareolus","Neomys fodiens","Oryctolagus cuniculus",
            "Rattus norvegicus")

## Create vectors of sample IDs associated with each predator
otter <- droplevels(subset(metadata, select="Sequencing_ID", 
                           Predator == "Otter"))
mink <- droplevels(subset(metadata, select="Sequencing_ID", 
                          Predator == "Mink"))

## Subset final dataframe for samples belonging to otter and mink
group.otter <- final.df[colnames(final.df) %in% otter$Sequencing_ID]
group.mink <- final.df[colnames(final.df) %in% mink$Sequencing_ID]

## Make row names first column in dataframe
group.otter <- tibble::rownames_to_column(group.otter, "Assignment")
group.mink <- tibble::rownames_to_column(group.mink, "Assignment")

## Remove empty taxonomic assignments
group.otter <- group.otter[rowSums(group.otter[, -1]) > 0,]
group.mink <- group.mink[rowSums(group.mink[, -1]) > 0,]

## Remove mink reads from otter data and vice versa
group.otter <- group.otter[which(!grepl("Neovison", group.otter$Assignment)),]
rownames(group.otter) <- NULL

group.mink <- group.mink[which(!grepl("Lutra", group.mink$Assignment)),]
rownames(group.mink) <- NULL

## Create new column containing the group each taxon belongs to
group.otter$Group <- factor(ifelse(group.otter$Assignment %in% amphibian, "Amphibian",
                                      ifelse(group.otter$Assignment %in% fish, "Fish",
                                             ifelse(group.otter$Assignment %in% bird, "Bird",
                                                    ifelse(group.otter$Assignment %in% mammal, "Mammal",
                                                           "Predator")))))
group.otter <- group.otter[,c(1,173,2:172)]

group.mink$Group <- factor(ifelse(group.mink$Assignment %in% amphibian, "Amphibian",
                                   ifelse(group.mink$Assignment %in% fish, "Fish",
                                          ifelse(group.mink$Assignment %in% bird, "Bird",
                                                 ifelse(group.mink$Assignment %in% mammal, "Mammal",
                                                        "Predator")))))
group.mink <- group.mink[,c(1,21,2:20)]

## Remove Assignment column
group.otter <- group.otter[,-1]
group.mink <- group.mink[,-1]

## Merge read counts by vertebrate group
group.otter <- ddply(group.otter, .(Group), numcolwise(sum))
group.mink <- ddply(group.mink, .(Group), numcolwise(sum))

## Calculate total read counts for each group
group.otter$otter_total <- rowSums(group.otter[,2:172])
group.mink$mink_total <- rowSums(group.mink[,2:20])

## Remove samples
group.otter <- group.otter[,-c(2:172)]
group.mink <- group.mink[,-c(2:20)]

## Combine totals for each predator
group.temp <- smartbind(group.otter, group.mink)
group.temp[is.na(group.temp)] <- 0
group.temp <- ddply(group.temp, .(Group), numcolwise(sum))

## Rename columns
colnames(group.temp) <- c("Group", "Otter", "Mink")

## Remove Group column from dataframe
group.prop <- group.temp[,-1]

## Calculate total number of reads for each predator
group.prop <- rbind(group.prop, colSums(group.prop))

## Create new dataframe containing the proportion of reads in each group
group.prop <- (group.prop/c(group.prop[6,]))*100

## Add Group column back to dataframe
group.prop <- cbind(group.temp$Group, group.prop[1:5,])

## Rename columns
colnames(group.prop)[1] <- "Group"

## Otter pie chart
pie1 <- ggplot(group.prop, aes(x="", y=Otter, fill=Group))
pie1 <- pie1 + labs(title="Otter (n = 171)")
pie1 <- pie1 + geom_bar(width = 1, stat = "identity", colour="black")
pie1 <- pie1 + coord_polar("y", start=0)
pie1 <- pie1 + geom_text(aes(label=percent(Otter/100)),
                         position = position_stack(vjust = 0.5),
                         cex=5)
pie1 <- pie1 + scale_fill_manual(values=c("limegreen","goldenrod1",
                                          "dodgerblue","hotpink",
                                          "grey60"))
pie1 <- pie1 + theme_bw()
pie1 <- pie1 + theme(panel.border = element_blank(),
                     panel.grid=element_blank(),
                     plot.title = element_text(face="bold", hjust=0.5, colour="black"),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20),
                     legend.key.size = unit(1, 'lines'))
pie1

## Mink pie chart
pie2 <- ggplot(group.prop, aes(x="", y=Mink, fill=Group))
pie2 <- pie2 + labs(title="Mink (n = 19)")
pie2 <- pie2 + geom_bar(width = 1, stat = "identity", colour="black")
pie2 <- pie2 + coord_polar("y", start=0)
pie2 <- pie2 + geom_text(aes(label=percent(Mink/100)), 
                         position = position_stack(vjust = 0.5),
                         cex=5)
pie2 <- pie2 + scale_fill_manual(values=c("limegreen","goldenrod1",
                                          "dodgerblue","hotpink",
                                          "grey60"))
pie2 <- pie2 + theme_bw()
pie2 <- pie2 + theme(panel.border = element_blank(),
                     panel.grid=element_blank(),
                     plot.title = element_text(face="bold", hjust=0.5, colour="black"),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20),
                     legend.key.size = unit(1, 'lines'))
pie2


###############################################################
# Pie chart: read counts by vertebrate group (exc. predators) #
###############################################################

## Create vectors of sample IDs associated with each predator
otter <- droplevels(subset(metadata, select="Sequencing_ID", 
                           Predator == "Otter"))
mink <- droplevels(subset(metadata, select="Sequencing_ID", 
                          Predator == "Mink"))

## Create copies of original dataframes for otter and mink
group.otter <- final.df[colnames(final.df) %in% otter$Sequencing_ID]
group.mink <- final.df[colnames(final.df) %in% mink$Sequencing_ID]

## Make row names first column in dataframe
group.otter <- tibble::rownames_to_column(group.otter, "Assignment")
group.mink <- tibble::rownames_to_column(group.mink, "Assignment")

## Remove empty taxonomic assignments
group.otter <- group.otter[rowSums(group.otter[, -1]) > 0,]
group.mink <- group.mink[rowSums(group.mink[, -1]) > 0,]

## Remove predators from dataframes
group.otter <- group.otter[which(!grepl("Lutra|Neovison", group.otter$Assignment)),]
rownames(group.otter) <- NULL

group.mink <- group.mink[which(!grepl("Lutra|Neovison", group.mink$Assignment)),]
rownames(group.mink) <- NULL

## Create new column containing the group each species belongs to
group.otter$Group <- factor(ifelse(group.otter$Assignment %in% amphibian, "Amphibian",
                                   ifelse(group.otter$Assignment %in% fish, "Fish",
                                          ifelse(group.otter$Assignment %in% bird, "Bird",
                                                 ifelse(group.otter$Assignment %in% mammal, "Mammal",
                                                        "NA")))))
group.otter <- group.otter[,c(1,173,2:172)]

group.mink$Group <- factor(ifelse(group.mink$Assignment %in% amphibian, "Amphibian",
                                  ifelse(group.mink$Assignment %in% fish, "Fish",
                                         ifelse(group.mink$Assignment %in% bird, "Bird",
                                                ifelse(group.mink$Assignment %in% mammal, "Mammal",
                                                       "NA")))))
group.mink <- group.mink[,c(1,21,2:20)]

## Remove species-level assignment
group.otter <- group.otter[,-1]
group.mink <- group.mink[,-1]

## Merge read counts by group
group.otter <- ddply(group.otter, .(Group), numcolwise(sum))
group.mink <- ddply(group.mink, .(Group), numcolwise(sum))

## Calculate total read counts for each group
group.otter$otter_total <- rowSums(group.otter[,2:172])
group.mink$mink_total <- rowSums(group.mink[,2:20])

## Remove samples
group.otter <- group.otter[,-c(2:172)]
group.mink <- group.mink[,-c(2:20)]

## Combine totals for each predator
group.temp <- smartbind(group.otter, group.mink)
group.temp[is.na(group.temp)] <- 0
group.temp <- ddply(group.temp, .(Group), numcolwise(sum))

## Rename columns
colnames(group.temp) <- c("Group", "Otter", "Mink")

## Remove Group column from dataframe
group.prop <- group.temp[,-1]

## Calculate total number of reads for each predator
group.prop <- rbind(group.prop, colSums(group.prop))

## Create new dataframe containing the proportion of reads in each group
group.prop <- (group.prop/c(group.prop[5,]))*100

## Add Group column back to dataframe
group.prop <- cbind(group.temp$Group, group.prop[1:4,])

## Rename columns
colnames(group.prop)[1] <- "Group"

## Otter pie chart
pie3 <- ggplot(group.prop, aes(x="", y=Otter, fill=Group))
pie3 <- pie3 + labs(title="")
pie3 <- pie3 + geom_bar(width = 1, stat = "identity", colour="black")
pie3 <- pie3 + coord_polar("y", start=0)
pie3 <- pie3 + geom_text(aes(label=percent(Otter/100)),
                         position = position_stack(vjust = 0.5),
                         cex=5)
pie3 <- pie3 + scale_fill_manual(values=c("limegreen","goldenrod1",
                                          "dodgerblue","hotpink"))
pie3 <- pie3 + theme_bw()
pie3 <- pie3 + theme(panel.border = element_blank(),
                     panel.grid=element_blank(),
                     plot.title = element_text(face="bold", hjust=0.5, colour="black"),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20),
                     legend.key.size = unit(1, 'lines'))
pie3

## Mink pie chart
pie4 <- ggplot(group.prop, aes(x="", y=Mink, fill=Group))
pie4 <- pie4 + labs(title="")
pie4 <- pie4 + geom_bar(width = 1, stat = "identity", colour="black")
pie4 <- pie4 + coord_polar("y", start=0)
pie4 <- pie4 + geom_text(aes(label=percent(Mink/100)), 
                         position = position_stack(vjust = 0.5),
                         cex=5)
pie4 <- pie4 + scale_fill_manual(values=c("limegreen","goldenrod1",
                                          "dodgerblue","hotpink"))
pie4 <- pie4 + theme_bw()
pie4 <- pie4 + theme(panel.border = element_blank(),
                     panel.grid=element_blank(),
                     plot.title = element_text(face="bold", hjust=0.5, colour="black"),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     text = element_text(size=20),
                     legend.key.size = unit(1, 'lines'))
pie4

## Plot pie charts alongside each other
g2 <- ggarrange(pie1, pie2, pie3, pie4, ncol=2, nrow=2,
                common.legend=TRUE, legend="bottom")

#ggsave(filename="Figures/Fig1.svg", 
#       plot = g2, width = 8, height = 8, dpi = 300, units = "in")


####################################################
# Heatmap: proportional read counts for each taxon #
####################################################

## Remove predators from dataframe
spp.temp <- final.df[which(!grepl("Lutra|Neovison", final.df$Assignment)),]

## Remove Assignment column
spp.temp <- spp.temp[,-1]

## Calculate total prey read counts for each sample
spp.temp <- rbind(spp.temp, colSums(spp.temp))

## Calculate proportional read counts
spp.prop <- spp.temp/c(spp.temp[43,])

## Remove row containing total proportional read counts
spp.prop <- spp.prop[-43,]

## Make row names first column in dataframe
spp.prop <- tibble:::rownames_to_column(spp.prop, "Assignment")

## Create column specifying which vertebrate group taxa belong to
spp.prop$Group <- ifelse(spp.prop$Assignment %in% amphibian, "Amphibian",
                         ifelse(spp.prop$Assignment %in% fish, "Fish",
                                ifelse(spp.prop$Assignment %in% bird, "Bird", "Mammal")))

## Move new column to start of dataframe
spp.prop <- spp.prop[,c(1,192,2:191)]

## Melt dataframe
prc.plot <- melt(spp.prop, id=c("Assignment","Group"))

## Rename columns
colnames(prc.plot)[3:4] <- c("Sample","PRC")

## Create dataframe containing sample ID, predator, and site information
id.pred.site <- data.frame(metadata$Sequencing_ID, 
                           metadata$Site,
                           metadata$Predator)
colnames(id.pred.site) <- c("Sample","Site","Predator")

## Merge predator information with proportional read count data
prc.plot <- merge(prc.plot, id.pred.site, by="Sample")

## Reorder dataframe
prc.plot <- prc.plot[,c(1,5:6,3,2,4)]

## Only keep taxa that have reads
prc.plot <- filter(prc.plot, PRC > 0)

## Sort dataframe by Group then Assignment
prc.plot <- with(prc.plot, prc.plot[order(Group, Assignment),])

## Order dataframe by Site then Sample
prc.plot <- prc.plot[order(prc.plot$Site, prc.plot$Sample),]
prc.plot$Sample <- factor(prc.plot$Sample, levels=unique(prc.plot$Sample))

## Create factor to order prey taxa by vertebrate group
prc.plot$fAssignment <- factor(prc.plot$Assignment, 
                               levels= c("Bufo bufo","Pelophylax ridibundus",
                                         "Rana temporaria","Abramis brama",
                                         "Anguilla anguilla","Barbatula barbatula",
                                         "Carassius carassius","Cottus gobio", 
                                         "Esox lucius","Gasterosteus aculeatus",
                                         "Gobio gobio","Lampetra fluviatilis",
                                         "Oncorhynchus mykiss","Percidae spp.",
                                         "Phoxinus phoxinus","Platichthys flesus",
                                         "Pungitius pungitius","Rutilus rutilus",
                                         "Salmo trutta","Scardinius erythrophthalmus",
                                         "Thymallus thymallus","Tinca tinca",
                                         "Alectoris rufa","Anas spp.",
                                         "Aythya spp.","Columba oenas",
                                         "Fulica atra","Gallinula chloropus",
                                         "Garrulus glandarius","Laridae spp.",
                                         "Phalacrocorax carbo","Phasianus colchicus",
                                         "Scolopax rusticola","Sturnus vulgaris",
                                         "Tringa spp.","Arvicola amphibius",
                                         "Lepus europaeus","Microtus spp.",
                                         "Myodes glareolus","Neomys fodiens",
                                         "Oryctolagus cuniculus","Rattus norvegicus"))

## Plot:
hm3 <- ggplot(prc.plot, aes(x=Sample, 
                            y=fct_rev(as_factor(fAssignment)), 
                            fill=PRC)) + 
        geom_tile(aes(colour=Group), size=0.5) + 
        scale_colour_manual(values=c("limegreen","goldenrod1",
                                     "dodgerblue","hotpink")) + 
        scale_fill_gradient(name="Proportional\nread counts", 
                            limits=c(0,1), 
                            breaks=c(0,0.25,0.5,0.75,1),
                            low="white", high="grey30",
                            guide=guide_colourbar(frame.colour="black",
                                                  ticks.colour="black")) + 
        labs(x="Sample", y="Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour="black", face = "italic"),
              text = element_text(size=20),
              legend.key.size = unit(1, 'lines')) + 
        facet_grid(Site ~ Predator, scales="free", space="free")
hm3

#ggsave(filename="Figures/FigS6.png", 
#       plot = hm3, width = 15, height = 15, dpi = 300, units = "in")


############################################
# Barplot: sample occupancy for each taxon #
############################################

## Duplicate original dataframe
spp.SO <- final.df

## Remove predators from data
spp.SO <- spp.SO[which(!grepl("Lutra|Neovison", spp.SO$Assignment)),]

## Reset row names
rownames(spp.SO) <- NULL

## Transpose dataframe
spp.SO <- setNames(data.frame(t(spp.SO[,-1])), spp.SO[,1])

## Make rownames first column in dataframe
spp.SO <- tibble::rownames_to_column(spp.SO, "Sample")

## Merge predator with read count data
spp.SO <- merge(spp.SO, id.pred.site, by="Sample")
spp.SO <- spp.SO[,c(1,44:45,2:43)]

## Count the number of samples belonging to each predator
sum(spp.SO$Predator == "Otter")
sum(spp.SO$Predator == "Mink")

## Create dataframe for each predator containing prey taxa and the total 
## number of samples they were found in
otter.SO <- subset(spp.SO, Predator=="Otter")
rownames(otter.SO) <- otter.SO$Sample
otter.SO <- otter.SO[,-c(1:3)]
otter.SO[otter.SO > 0] <- 1
otter.SO <- data.frame(colnames(otter.SO), colSums(otter.SO))
colnames(otter.SO) <- c("Assignment", "Otter_samples")
otter.SO$Otter_total <- 171

mink.SO <- subset(spp.SO, Predator=="Mink")
rownames(mink.SO) <- mink.SO$Sample
mink.SO <- mink.SO[,-c(1:3)]
mink.SO[mink.SO > 0] <- 1
mink.SO <- data.frame(colnames(mink.SO), colSums(mink.SO))
colnames(mink.SO) <- c("Assignment", "Mink_samples")
mink.SO$Mink_total <- 19

## Create columns containing proportion of predator samples that an 
## assignment was detected
otter.SO$Otter <- (otter.SO$Otter_samples/otter.SO$Otter_total)*100
mink.SO$Mink <- (mink.SO$Mink_samples/mink.SO$Mink_total)*100

## Manipulate dataframes for plotting
otter.plot <- otter.SO[,-3]
mink.plot <- mink.SO[,-3]

## Create columns specifying which vertebrate group prey taxa belong to
otter.plot$Group <- ifelse(otter.plot$Assignment %in% amphibian, "Amphibian",
                           ifelse(otter.plot$Assignment %in% fish, "Fish",
                                  ifelse(otter.plot$Assignment %in% bird, "Bird", "Mammal")))

mink.plot$Group <- ifelse(mink.plot$Assignment %in% amphibian, "Amphibian",
                           ifelse(mink.plot$Assignment %in% fish, "Fish",
                                  ifelse(mink.plot$Assignment %in% bird, "Bird", "Mammal")))

## Melt dataframes for plotting
otter.plot <- melt(otter.plot, id=c("Assignment","Group","Otter_samples"))
mink.plot <- melt(mink.plot, id=c("Assignment","Group","Mink_samples"))

## Rename columns
colnames(otter.plot)[3:5] <- c("Count","Predator","Proportion")
colnames(mink.plot)[3:5] <- c("Count","Predator","Proportion")

## Combine otter and mink data
SO.plot <- rbind(otter.plot, mink.plot)

## Only keep prey taxa found in diet of a predator
SO.plot <- filter(SO.plot, Proportion > 0)

## Create factor to order prey taxa by vertebrate group
SO.plot$fAssignment <- factor(SO.plot$Assignment, 
                              levels= c("Bufo bufo","Pelophylax ridibundus",
                                        "Rana temporaria","Abramis brama",
                                        "Anguilla anguilla","Barbatula barbatula",
                                        "Carassius carassius","Cottus gobio", 
                                        "Esox lucius","Gasterosteus aculeatus",
                                        "Gobio gobio","Lampetra fluviatilis",
                                        "Oncorhynchus mykiss","Percidae spp.",
                                        "Phoxinus phoxinus","Platichthys flesus",
                                        "Pungitius pungitius","Rutilus rutilus",
                                        "Salmo trutta","Scardinius erythrophthalmus",
                                        "Thymallus thymallus","Tinca tinca",
                                        "Alectoris rufa","Anas spp.",
                                        "Aythya spp.","Columba oenas",
                                        "Fulica atra","Gallinula chloropus",
                                        "Garrulus glandarius","Laridae spp.",
                                        "Phalacrocorax carbo","Phasianus colchicus",
                                        "Scolopax rusticola","Sturnus vulgaris",
                                        "Tringa spp.","Arvicola amphibius",
                                        "Lepus europaeus","Microtus spp.",
                                        "Myodes glareolus","Neomys fodiens",
                                        "Oryctolagus cuniculus","Rattus norvegicus"))

## Create custom facet labels for plotting
predator.names <- c("Otter" = "Otter (n = 171)",
                    "Mink" = "Mink (n = 19)")

## Plot sample occupancy
p1 <- ggplot(SO.plot, aes(x=fAssignment, y=Proportion, fill=Group)) + 
        geom_bar(stat="identity", colour="black") + 
        scale_y_continuous(limits=c(0,100),
                           breaks=seq(0,100,10)) + 
        scale_fill_manual(name="Group",
                          values=c("limegreen","goldenrod1",
                                   "dodgerblue","hotpink")) + 
        labs(x="Taxonomic assignment", y="Samples (%)") + 
        geom_text(aes(label=Count, y=Proportion), 
                  stat="identity", vjust = -.5) + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", angle = 60, vjust=1, hjust=1),
              axis.text.y = element_text(colour="black"),
              text = element_text(size=16),
              legend.position="bottom") + 
        facet_grid(Predator ~ ., labeller=labeller(Predator=predator.names))
p1



####################################################
# Barplot: sample occupancy for each taxon by site #
####################################################

## Duplicate final dataframe
site.SO <- final.df

## Remove predators from dataframe
site.SO <- site.SO[which(!grepl("Lutra|Neovison", site.SO$Assignment)),]

## Reset row names
rownames(site.SO) <- NULL

## Transpose dataframe
site.SO <- setNames(data.frame(t(site.SO[,-1])), site.SO[,1])

## Make rownames first column in dataframe
site.SO <- tibble::rownames_to_column(site.SO, "Sample")

## Merge predator information with read count data
site.SO <- merge(site.SO, id.pred.site, by="Sample")
site.SO <- site.SO[,c(1,44:45,2:43)]

## Count the number of samples belonging to each predator
sum(site.SO$Predator == "Otter")
sum(site.SO$Predator == "Mink")

## Create dataframe for each predator at each site containing prey taxa 
## and total number of samples they were found in
otter.RG <- subset(site.SO, Predator=="Otter" & Site=="River Glaven")
otter.MT <- subset(site.SO, Predator=="Otter" & Site=="Malham Tarn")
otter.RH <- subset(site.SO, Predator=="Otter" & Site=="River Hull")

mink.RG <- subset(site.SO, Predator=="Mink" & Site=="River Glaven")
mink.RH <- subset(site.SO, Predator=="Mink" & Site=="River Hull")

## Calculate how many samples each prey taxa was found in
rownames(otter.RG) <- otter.RG$Sample
otter.RG <- otter.RG[,-c(1:3)]
otter.RG[otter.RG > 0] <- 1
otter.RG <- data.frame(colnames(otter.RG), colSums(otter.RG))
colnames(otter.RG) <- c("Assignment", "Count")
otter.RG$Otter_total <- 36

rownames(otter.MT) <- otter.MT$Sample
otter.MT <- otter.MT[,-c(1:3)]
otter.MT[otter.MT > 0] <- 1
otter.MT <- data.frame(colnames(otter.MT), colSums(otter.MT))
colnames(otter.MT) <- c("Assignment", "Count")
otter.MT$Otter_total <- 25

rownames(otter.RH) <- otter.RH$Sample
otter.RH <- otter.RH[,-c(1:3)]
otter.RH[otter.RH > 0] <- 1
otter.RH <- data.frame(colnames(otter.RH), colSums(otter.RH))
colnames(otter.RH) <- c("Assignment", "Count")
otter.RH$Otter_total <- 110

rownames(mink.RG) <- mink.RG$Sample
mink.RG <- mink.RG[,-c(1:3)]
mink.RG[mink.RG > 0] <- 1
mink.RG <- data.frame(colnames(mink.RG), colSums(mink.RG))
colnames(mink.RG) <- c("Assignment", "Count")
mink.RG$Mink_total <- 3

rownames(mink.RH) <- mink.RH$Sample
mink.RH <- mink.RH[,-c(1:3)]
mink.RH[mink.RH > 0] <- 1
mink.RH <- data.frame(colnames(mink.RH), colSums(mink.RH))
colnames(mink.RH) <- c("Assignment", "Count")
mink.RH$Mink_total <- 16

## Create columns containing proportion of predator samples that prey 
## taxa were detected
otter.RG$RG <- (otter.RG$Count/otter.RG$Otter_total)*100
otter.MT$MT <- (otter.MT$Count/otter.MT$Otter_total)*100
otter.RH$RH <- (otter.RH$Count/otter.RH$Otter_total)*100

mink.RG$RG <- (mink.RG$Count/mink.RG$Mink_total)*100
mink.RH$RH <- (mink.RH$Count/mink.RH$Mink_total)*100

## Manipulate dataframes for plotting
otter.RG <- otter.RG[,-3]
otter.MT <- otter.MT[,-3]
otter.RH <- otter.RH[,-3]

mink.RG <- mink.RG[,-3]
mink.RH <- mink.RH[,-3]

## Create columns specifying which group prey taxa belong to
otter.RG$Group <- ifelse(otter.RG$Assignment %in% amphibian, "Amphibian",
                           ifelse(otter.RG$Assignment %in% fish, "Fish",
                                  ifelse(otter.RG$Assignment %in% bird, "Bird", "Mammal")))

otter.MT$Group <- ifelse(otter.MT$Assignment %in% amphibian, "Amphibian",
                         ifelse(otter.MT$Assignment %in% fish, "Fish",
                                ifelse(otter.MT$Assignment %in% bird, "Bird", "Mammal")))

otter.RH$Group <- ifelse(otter.RH$Assignment %in% amphibian, "Amphibian",
                         ifelse(otter.RH$Assignment %in% fish, "Fish",
                                ifelse(otter.RH$Assignment %in% bird, "Bird", "Mammal")))

mink.RG$Group <- ifelse(mink.RG$Assignment %in% amphibian, "Amphibian",
                         ifelse(mink.RG$Assignment %in% fish, "Fish",
                                ifelse(mink.RG$Assignment %in% bird, "Bird", "Mammal")))

mink.RH$Group <- ifelse(mink.RH$Assignment %in% amphibian, "Amphibian",
                          ifelse(mink.RH$Assignment %in% fish, "Fish",
                                 ifelse(mink.RH$Assignment %in% bird, "Bird", "Mammal")))

## Add column specifying predator
otter.RG$Predator <- "Otter"
otter.MT$Predator <- "Otter"
otter.RH$Predator <- "Otter"

mink.RG$Predator <- "Mink"
mink.RH$Predator <- "Mink"

## Melt dataframes for plotting
otter.RG <- melt(otter.RG, id=c("Assignment","Group","Predator","Count"))
otter.MT <- melt(otter.MT, id=c("Assignment","Group","Predator","Count"))
otter.RH <- melt(otter.RH, id=c("Assignment","Group","Predator","Count"))

mink.RG <- melt(mink.RG, id=c("Assignment","Group","Predator","Count"))
mink.RH <- melt(mink.RH, id=c("Assignment","Group","Predator","Count"))

## Combine data
site.plot <- rbind(otter.RG, otter.MT, otter.RH, mink.RG, mink.RH)

## Rename columns
colnames(site.plot)[5:6] <- c("Site","Proportion")

## Only keep prey taxa found in diet of a predator
site.plot <- filter(site.plot, Proportion > 0)

## Create factor to order prey taxa by vertebrate group
site.plot$fAssignment <- factor(site.plot$Assignment, 
                                levels= c("Bufo bufo","Pelophylax ridibundus",
                                          "Rana temporaria","Barbatula barbatula",
                                          "Cottus gobio","Gasterosteus aculeatus",
                                          "Gobio gobio","Phoxinus phoxinus",
                                          "Pungitius pungitius","Carassius carassius", 
                                          "Lampetra fluviatilis","Percidae spp.",
                                          "Platichthys flesus","Rutilus rutilus",
                                          "Scardinius erythrophthalmus","Abramis brama",
                                          "Anguilla anguilla","Esox lucius",
                                          "Oncorhynchus mykiss","Salmo trutta",
                                          "Thymallus thymallus","Tinca tinca",
                                          "Alectoris rufa","Anas spp.",
                                          "Aythya spp.","Columba oenas",
                                          "Fulica atra","Gallinula chloropus",
                                          "Garrulus glandarius","Laridae spp.",
                                          "Phalacrocorax carbo","Phasianus colchicus",
                                          "Scolopax rusticola","Sturnus vulgaris",
                                          "Tringa spp.","Arvicola amphibius",
                                          "Lepus europaeus","Microtus spp.",
                                          "Myodes glareolus","Neomys fodiens",
                                          "Oryctolagus cuniculus","Rattus norvegicus"))

## Unabbreviate site names
site.plot$Site <- gsub("RG", "River Glaven", site.plot$Site)
site.plot$Site <- gsub("MT", "Malham Tarn", site.plot$Site)
site.plot$Site <- gsub("RH", "River Hull", site.plot$Site)

## Plot sample occupancy by site
p2 <- ggplot(site.plot, aes(x=fAssignment, y=Proportion, fill=fAssignment)) + 
        geom_bar(stat="identity", colour="black") + 
        scale_y_continuous(limits=c(0,100),
                           breaks=seq(0,100,10)) + 
        scale_fill_manual(breaks=c("Bufo bufo","Pelophylax ridibundus",
                                   "Rana temporaria","Barbatula barbatula",
                                   "Cottus gobio","Gasterosteus aculeatus",
                                   "Gobio gobio","Phoxinus phoxinus",
                                   "Pungitius pungitius","Carassius carassius", 
                                   "Lampetra fluviatilis","Percidae spp.",
                                   "Platichthys flesus","Rutilus rutilus",
                                   "Scardinius erythrophthalmus","Abramis brama",
                                   "Anguilla anguilla","Esox lucius",
                                   "Oncorhynchus mykiss","Salmo trutta",
                                   "Thymallus thymallus","Tinca tinca",
                                   "Alectoris rufa","Anas spp.",
                                   "Aythya spp.","Columba oenas",
                                   "Fulica atra","Gallinula chloropus",
                                   "Garrulus glandarius","Laridae spp.",
                                   "Phalacrocorax carbo","Phasianus colchicus",
                                   "Scolopax rusticola","Sturnus vulgaris",
                                   "Tringa spp.","Arvicola amphibius",
                                   "Lepus europaeus","Microtus spp.",
                                   "Myodes glareolus","Neomys fodiens",
                                   "Oryctolagus cuniculus","Rattus norvegicus"),
                          values=c("limegreen","limegreen",
                                   "limegreen","lightskyblue",
                                   "lightskyblue","lightskyblue",
                                   "lightskyblue","lightskyblue",
                                   "lightskyblue","deepskyblue3",
                                   "deepskyblue3","deepskyblue3",
                                   "deepskyblue3","deepskyblue3",
                                   "deepskyblue3","deepskyblue4",
                                   "deepskyblue4","deepskyblue4",
                                   "deepskyblue4","deepskyblue4",
                                   "deepskyblue4","deepskyblue4",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","goldenrod1",
                                   "goldenrod1","hotpink",
                                   "hotpink","hotpink",
                                   "hotpink","hotpink",
                                   "hotpink","hotpink")) + 
        labs(x="Taxon\n\n", y="Frequency of occurrence (%)") + 
        geom_text(aes(label=Count, y=Proportion), stat="identity", vjust = -.5) + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", face="italic", angle = 60, vjust=1, hjust=1),
              axis.text.y = element_text(colour="black"),
              text = element_text(size=16),
              legend.position="none") + 
        facet_grid(Predator ~ Site, scales="free", space="free")
p2

#ggsave(filename="Figures/Fig3.svg", 
#       plot = p2, width = 16, height = 10, dpi = 300, units = "in")



#' ---
#'
#' ## 10) Trophic network
#'
#' Use 'bipartite' to construct a trophic network and obtain network
#' metrics, specifically those relating to the predators otter and 
#' mink.
#' 

## Transpose data into site x taxon matrix
diet.df <- setNames(data.frame(t(final.df[,-1])), final.df[,1])

## Duplicate data and convert to presence/absence
diet.pa <- diet.df
diet.pa[diet.pa > 0] <- 1
diet.pa <- tibble::rownames_to_column(diet.pa, "Sequencing_ID")

## Remove predators from dataframe
diet.pa <- diet.pa[,which(!grepl("Lutra|Neovison", colnames(diet.pa)))]

## Duplicate metadata
diet.metadata <- metadata
rownames(diet.metadata) <- NULL

## Combine metadata with presence-absence data
diet.pa <- merge(diet.metadata, diet.pa, by="Sequencing_ID")

## Remove redundant columns
diet.pa <- diet.pa[,-c(1:6)]

## Collapse data by species
interactions <- ddply(diet.pa, .(Predator), numcolwise(sum))

## Transpose dataframe
## Make first column row names
interactions <- tibble:::column_to_rownames(interactions, "Predator")
network <- data.frame(t(interactions))

## Create vectors that contain order species are to be plotted in
seq.high <- c("Otter","Mink")
seq.low <- c("Bufo bufo","Pelophylax ridibundus",
             "Rana temporaria","Barbatula barbatula",
             "Cottus gobio","Gasterosteus aculeatus",
             "Gobio gobio","Phoxinus phoxinus",
             "Pungitius pungitius","Carassius carassius", 
             "Lampetra fluviatilis","Percidae spp.",
             "Platichthys flesus","Rutilus rutilus",
             "Scardinius erythrophthalmus","Abramis brama",
             "Anguilla anguilla","Esox lucius",
             "Oncorhynchus mykiss","Salmo trutta",
             "Thymallus thymallus","Tinca tinca",
             "Alectoris rufa","Anas spp.",
             "Aythya spp.","Columba oenas",
             "Fulica atra","Gallinula chloropus",
             "Garrulus glandarius","Laridae spp.",
             "Phalacrocorax carbo","Phasianus colchicus",
             "Scolopax rusticola","Sturnus vulgaris",
             "Tringa spp.","Arvicola amphibius",
             "Lepus europaeus","Microtus spp.",
             "Myodes glareolus","Neomys fodiens",
             "Oryctolagus cuniculus","Rattus norvegicus")

## Create vector of colours for lower trophic levels
low.trophic <- c("limegreen","limegreen", "limegreen",
                 "lightskyblue","lightskyblue","lightskyblue",
                 "lightskyblue","lightskyblue","lightskyblue",
                 "deepskyblue3","deepskyblue3","deepskyblue3",
                 "deepskyblue3","deepskyblue3","deepskyblue3",
                 "deepskyblue4","deepskyblue4","deepskyblue4",
                 "deepskyblue4","deepskyblue4","deepskyblue4",
                 "deepskyblue4","goldenrod1","goldenrod1",
                 "goldenrod1","goldenrod1","goldenrod1",
                 "goldenrod1","goldenrod1","goldenrod1",
                 "goldenrod1","goldenrod1","goldenrod1",
                 "goldenrod1","goldenrod1","hotpink",
                 "hotpink","hotpink","hotpink",
                 "hotpink","hotpink","hotpink")

## Create vector for legend keys
name.legend1=c("Amphibian","Fish","Bird","Mammal")
col.legend1=c("limegreen","lightskyblue","goldenrod1","hotpink")

name.legend2=c("Small","Medium","Large")
col.legend2=c("lightskyblue","deepskyblue3","deepskyblue4")

## Plot network
plotweb(sortweb(network, sort.order="seq",
                sequence=list(seq.higher=seq.high, seq.lower=seq.low)),
        method="normal", labsize=1, high.spacing=0.2, low.spacing=0.01, 
        arrow="up.center", col.interaction=c("grey80"), 
        col.high="black", col.low=low.trophic,
        text.rot=90, x.lim=c(0,1.34), y.lim=c(-1,3.5))
legend(x=0.5, y=3.5, xpd=T, 
       legend=name.legend1, fill=col.legend1, bty="n", 
       border="black", x.intersp=0.5, text.width=0.07)
legend(x=0.8, y=3.5, xpd=T, 
       legend=name.legend2, fill=col.legend2, bty="n", 
       border="black", x.intersp=0.5, text.width=0.07)

## Plot number of prey taxa consumed by each predator
visweb(network, type="diagonal", square="interaction",
       text="interaction", frame=TRUE, 
       labsize=1, plotsize=30)

## Network-level analysis
networklevel(network)

#'
#' Report: 
#' 
#' mean.number.of.shared.partners.HL = number of shared prey species (8)
#' 
#' Linkage density = the weighted diversity of interactions per species 
#' (7.70508904)
#' 
#' Weighted connectance = the weighted realised proportion of possible 
#' links, calculated as quantitative linkage density divided by the 
#' number of species in the network (0.184)
#' 
#' Weighted generality (generality.HL) = the mean effective number of 
#' prey species per predator weighted by their marginal totals (14.333)
#' 
#' Weighted network specialisation index H2’ (H2) = the degree of 
#' specialisation among prey and predators across an entire network, 
#' ranges between 0 (no specialisation) and 1 (complete specialisation) 
#' (0.628)
#' 
#' niche.overlap.HL = mean similarity in interaction pattern between 
#' species of the same level. Values near 0 indicate no common use of 
#' niches, 1 indicates perfect niche overlap (0.267)
#' 
#' C.score.HL = mean (normalised) number of checkerboard combinations 
#' across all species. Values close to 1 indicate that there is evidence 
#' for disaggregation, e.g. through competition. Values close to 0 
#' indicate aggregation of species, i.e. no repelling forces between 
#' species (1.000)
#' 


## Species-level analysis
specieslevel(network, level="higher")

#'
#' Report: 
#' 
#' PDI = Paired Differences Index, ranges between 0 (generalist)
#' and 1 (specialist)
#' 
#' resource.range = value of 0 when all resources are used, but a value 
#' of 1 when only one resource is used. It is, in fact, closer to an 
#' “unused resource range”
#' 
#' species.specificity.index = coefficient of variation of interactions, 
#' normalised to values between 0 and 1. Values of 0 indicate low, those 
#' of 1 a high variability (and hence suggesting low and high specificity).
#' Essentially, lower values indicate neither predator using one resource.
#' 
#' Fisher.alpha = Fisher’s alpha diversity for each species 
#' 
#' partner.diversity = Shannon diversity or per-species generality/
#' vulnerability of the interactions of each species
#' 
#' d = Specialisation of each species based on its discrimination from 
#' random selection of partners
#' 



#' ---
#'
#' ## 11) Alpha and beta diversity of prey communities
#'
#' Statistically compare alpha and beta diversity between predators, and
#' between otter samples from different sites. Data exploration (see below)
#' indicated there were not enough data for other comparisons. We will 
#' use 'vegan' and 'betapart' to examine taxon richness and similarity of 
#' prey communities.
#' 

## Remove predators from dataframe
diet.df <- diet.df[,which(!grepl("Lutra|Neovison", colnames(diet.df)))]

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
diet.df <- diet.df[!sapply(diet.df, function(x) all(x == 0))]
diet.df <- diet.df[!apply(diet.df == 0, 1, all),]

## NB: Two otter samples (MT19-1, I4) and two mink samples (O1, RH03) 
## did not contain any prey taxa

## Check structure of metadata
str(diet.metadata)

## Make all columns factors
diet.metadata[] <- lapply(diet.metadata, factor) 

## Basic summaries:
## Total number of sequences per sample
sum.of.rows <- apply(diet.df, 1, sum)
sort(sum.of.rows, decreasing = TRUE)
sum(sum.of.rows)

## Total number of sequences per taxon across all samples
sum.of.columns <- apply(diet.df, 2, sum)
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Number of samples where each taxon occurs
spec.pres <- apply(diet.df > 0, 2, sum) 
sort(spec.pres, decreasing = TRUE)

## Transform read counts to proportional read counts
diet.total <- decostand(diet.df, method="total", MARGIN=1) # by rows (samples)
head(diet.total[,1:3], n = 3)

## Calculate taxon richness
sample.richness <- specnumber(diet.total)
sample.richness

## Create data frame
alpha <- data.frame(sample.richness)
alpha <- tibble:::rownames_to_column(alpha, "Sequencing_ID")

## Add metadata from external file
alpha <- merge(alpha, diet.metadata, by="Sequencing_ID")

## Reset row names of data frame for further indexing
rownames(alpha) <- NULL

## Check structure of dataframe
str(alpha)


###########################################
# PREDATOR, SITE AND WATERBODY COMPARISON #
###########################################

#=================#
# ALPHA DIVERSITY #
#=================#

## Plot prey richness in samples according to site and waterbody.
p3 <- ggplot(alpha, aes(x=Predator, y=sample.richness)) + 
        geom_jitter(aes(colour=Predator), 
                    cex=2, width=0.2, show.legend=FALSE) + 
        geom_boxplot(alpha=0.7, outlier.shape=NA) + 
        scale_y_continuous(limits=c(0,10),
                              breaks=seq(0,10,2)) + 
        labs(x="Predator", y="Taxon richness") + 
        scale_colour_manual(values=c("limegreen","purple")) + 
        theme_bw() + 
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "none",
              text = element_text(size=20)) + 
        facet_grid(Waterbody ~ Site, scales="free", space="free")
p3

#ggsave(filename="Figures/FigS4.png", 
#       plot = p3, width = 13, height = 13, dpi = 300, units = "in")


## There is not enough data for mink or each type of waterbody to examine 
## differences in otter and mink diet with regard to habitat. Furthermore, 
## as the majority of samples for otter and mink came from rivers, the 
## differences observed will be reflected in the richness comparison between
## predators. 


#######################
# PREDATOR COMPARISON #
#######################

#=================#
# ALPHA DIVERSITY #
#=================#

## Check if data are normally distributed
## Histogram:
hist(alpha$sample.richness)  # roughly symmetrical

## Quantile-quantile plot (Q-Q plot):
qqnorm(alpha$sample.richness)
qqline(alpha$sample.richness, lty=2)  # observations don't tail away much

## Quantitative tests for normal distribution
## Shapiro-Wilk test for normality:
shapiro.test(alpha$sample.richness)  # P = 1.884e-08

## Kolmogorov-Smirnov test:
ks.test(alpha$sample.richness, pnorm)  # P < 2.2e-16

## Check whether any common data transformations improve data distribution
## Natural log:
hist(log(alpha$sample.richness))
shapiro.test(log(alpha$sample.richness))  # P = 1.303e-09

## log10:
hist(log10(alpha$sample.richness))
shapiro.test(log10(alpha$sample.richness))  # P = 1.303e-09

## Exponential:
hist(exp(alpha$sample.richness))
shapiro.test(exp(alpha$sample.richness))  # P < 2.2e-16

## logit:
hist(log(alpha$sample.richness/(alpha$sample.richness)))
shapiro.test(log(alpha$sample.richness/(alpha$sample.richness)))  # NA

## Square root:
hist(sqrt(alpha$sample.richness))
shapiro.test(sqrt(alpha$sample.richness))  # P = 1.391e-07

## Reciprocal:
hist(1/(alpha$sample.richness))
shapiro.test(1/(alpha$sample.richness))  # P = 8.935e-16

## Taxon richness data is not normally distributed and transformations 
## do not resolve this issue
## Compare variance in taxon richness between different predators
levene.test(alpha$sample.richness,
            alpha$Predator,
            location="mean")  # P = 0.01313

## It is unlikely that the data will conform to all the assumptions of a
## one-way ANOVA, but we will run the model and assess the residuals
anova <- lm(sample.richness ~ Predator, data=alpha)
summary(anova)  # summaries for each factor level
summary.aov(anova)  # summary for factor as a whole
model.tables(aov(anova), "means", se=TRUE)  # mean values for each factor level
TukeyHSD(aov(anova))  # Tukey post-hoc test
plot(TukeyHSD(aov(anova)))  # visualise pairwise comparisons

# Get standardised residuals for model validation
sresid <- resid(anova, type="pearson")

## Check assumption of normal distribution
## Most models are robust to slight deviations from normality in the 
## residuals
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 3.402e-07
## Not normally distributed

## Check assumption of homogeneity of variances
plot(sresid ~ anova$fitted.values)
plot(sresid ~ alpha$Predator)
## No heteroscedasticity present

## Check assumption of no collinearity
## Only one variable being modelled so no collinearity present.

## Check dataset does not contain serial auto-correlation
## This can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(anova)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation

## Check model not biased by unduly influential observations
influence <- influence.measures(anova)
summary(influence)
CD <- cooks.distance(anova)
plot(CD ~ sresid)
## Cook's distance <1 for so observations not exerting strong influence 
## on model parameters

## Several model assumptions are violated and data cannot be transformed to
## conform to a normal distribution, thus a one-way ANOVA is not applicable. 
## In these scenarios, switching to a Generalised Linear Model (GLM) and 
## changing the combination of error family and link-function terms to
## achieve normally distributed residuals is recommended. 

regression <- glm(sample.richness ~ Predator, 
                  family=poisson(link="log"), 
                  data=alpha)
summary(regression)
anova(regression, test = "Chi")
drop1(regression, test = "Chi")
TukeyHSD(aov(regression))

## Check model meets GLM assumptions
## Test for overdispersion
116.32/184
1-pchisq(116.32, df=184)  # not overdispersed

## Plot the fitted data against the observed data
plot(alpha$sample.richness ~ fitted(regression))

## Perform model validation checks to ensure model is good fit to data 
## and making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(regression)
shapiro.test(sresid) # P = 3.026e-07
## Some deviation from normality as residuals are not normally distributed
## therefore model is unreliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ alpha$Predator, pch=20, cex=2, cex.lab=1.5) 
## No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(regression)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(regression)
summary(influence)
CD <- cooks.distance(regression)
plot(CD ~ sresid)
## Cook's distance < 1 so observations not exerting strong influence on 
## model parameters

## Poisson distribution with log, identity, and square root link-functions 
## did not improve distribution of the residuals (code not repeated - 
## link function should be altered to see different results)
## We will use the non-parametric Kruskal-Wallis test instead 
## (NB: this tests for differences between the median values of different 
## groups, not the mean values)
kruskal.test(sample.richness ~ Predator, data=alpha)

## Dunn's test of multiple comparisons
dunnTest(sample.richness ~ Predator,
         data=alpha,
         method="none")

## KW test indicates significant difference in species richness between
## predators. Plot species richness in each sample according to predator:
p4 <- ggplot(alpha, aes(x=Predator, y=sample.richness, colour=Predator)) + 
        geom_jitter(cex=2, 
                    width=0.15,
                    height=0.2) + 
        geom_boxplot(colour="black",
                     alpha=0.7,
                     outlier.shape=NA) + 
        scale_y_continuous(limits=c(0,10),
                           breaks=seq(0,10,1)) + 
        annotate("text", x = c("Mink","Otter"), y = 10, 
                 label = c("a","b"), cex=10) + 
        labs(title=expression(bold("A "~alpha~"Diversity")),
             x="Predator", y="Taxon richness") + 
        scale_colour_manual(values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(size=25, margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(size=25, margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", size=25),
              axis.text.y = element_text(colour="black", size=25),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "none",
              legend.key = element_blank(),
              text = element_text(size=20))
p4


#==================================#
# RAREFACTION/EXTRAPOLATION CURVES #
#==================================#

## Look at whether differences in richness are influenced by sampling
## effort by performing rarefaction and extrapolation curves using the
## iNEXT package.

## Remove samples that contained no prey taxa from metadata: 
## I4, MT19-1 (otter), O1, RH03 (mink)
pred.metadata <- diet.metadata[diet.metadata$Sequencing_ID %in% rownames(diet.df),]
rownames(pred.metadata) <- NULL

## Create vectors of sample IDs associated with each predator.
otter <- droplevels(subset(pred.metadata, select="Sequencing_ID", Predator == "Otter"))
mink <- droplevels(subset(pred.metadata, select="Sequencing_ID", Predator == "Mink"))

## Transform read counts to presence/absence
pred.pa <- decostand(diet.df, method="pa", MARGIN=1) # by rows (samples)

## Create copies of diet.pa for otter and mink
otter.pa <- pred.pa[rownames(pred.pa) %in% otter$Sequencing_ID,]
mink.pa <- pred.pa[rownames(pred.pa) %in% mink$Sequencing_ID,]

## Calculate incidence frequency for each species detected in otter and
## mink samples
otter.richness <- colSums(otter.pa)
mink.richness <- colSums(mink.pa)

## Add sample size to beginning of each richness vector. This is because 
## the first entry of each list for iNEXT must be the total number of 
## sampling units, followed by the species incidence frequencies.
otter.richness <- append(otter.richness, 169, after=0)
mink.richness <- append(mink.richness, 17, after=0)

## Make list of otter and mink samples
richness <- list(otter.richness, mink.richness)
names(richness) <- c("Otter", "Mink")

## Run iNEXT function
## With 300 samples for each predator:
pred.re.300 <- iNEXT(richness, q=0, datatype="incidence_freq", 
                     endpoint=300, knots=60, se=TRUE, conf=0.95,
                     nboot=1000)
pred.re.300

## With 3500 samples for each predator:
pred.re.3500 <- iNEXT(richness, q=0, datatype="incidence_freq",
                      endpoint=3500, knots=700, se=TRUE, conf=0.95,
                      nboot=1000)
pred.re.3500

## Sample-size-based R/E curves
p5a <- ggiNEXT(pred.re.300, type=1) + 
        scale_colour_manual(values=c("limegreen","purple")) +
        scale_fill_manual(values=c("limegreen","purple")) +
        scale_shape_manual(values = c(19,19)) +
        scale_x_continuous(limits=c(0,300), breaks=seq(0,300,50)) +
        scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
        labs(title="B  Rarefaction/Extrapolation (R/E) curves",
             subtitle="i  Sample size-based R/E curve",
             y="Taxon diversity") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="bottom",
              legend.box="vertical",
              legend.key=element_blank(),
              text=element_text(size=20))
p5a

## Sample completeness curves
p5b <- ggiNEXT(pred.re.300, type=2) +
        scale_colour_manual(values=c("limegreen","purple")) +
        scale_fill_manual(values=c("limegreen","purple")) +
        scale_shape_manual(values = c(19,19)) +
        scale_x_continuous(limits=c(0,300), breaks=seq(0,300,50)) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
        labs(title="",
             subtitle="ii  Sample completeness curve") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="none",
              text=element_text(size=20))
p5b

## Coverage-based R/E curves
p5c <- ggiNEXT(pred.re.300, type=3) + 
        scale_colour_manual(values=c("limegreen","purple")) +
        scale_fill_manual(values=c("limegreen","purple")) +
        scale_shape_manual(values = c(19,19)) +
        scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
        scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
        labs(title="",
             subtitle="iii  Coverage-based R/E curve",
             y="Taxon diversity") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="none",
              text=element_text(size=20))
p5c

## Apply the estimateD() function to obtain diversity estimates of order 
## q = 0, 1, 2 for any particular level of sample size (base="size") or
## any specified level of sample coverage (base="coverage")
estimateD(richness, datatype="incidence_freq", base="size",
          level=NULL, conf=0.95)
estimateD(richness, datatype="incidence_freq", base="coverage", 
          level=0.95, conf=0.95)


#==================================#
# BETA DIVERSITY: DATA EXPLORATION #
#==================================#

## Different digestion rates and bias inherent to PCR amplification with 
## universal primers (i.e. variation in number of template-primer mismatches) 
## mean that some prey species will be preferentially amplified and reads may 
## not accurately reflect relative biomass and composition of prey in faecal 
## sample. DNA metabarcoding assessments of diet typically employ occurrence 
## data as these are less affected by taxon recovery bias than relative 
## read abundance data, but more weight is placed on minor prey items as well 
## as potential secondary predation and contaminants.
##
## We will compare the patterns produced by occurrence (i.e. Jaccard 
## dissimilarity) and relative read abundance data (i.e. Bray-Curtis
## dissimilarity) to see whether one index gives clearer insights to
## variables influencing community composition.
##
## Preliminary analysis of all data points indicated that two samples
## are extreme outliers and potentially skew visualisation of community
## dissimilarity. We will present the data with and without these outliers
## for comparison.

## Create copy of diet.df
pred.beta.jac <- diet.df

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
pred.beta.jac <- pred.beta.jac[!sapply(pred.beta.jac, function(x) all(x == 0))]
pred.beta.jac <- pred.beta.jac[!apply(pred.beta.jac == 0, 1, all),]

## Remove samples that are outliers: LIB02-TL01, LIB02-TL07
pred.beta.jac.no <- pred.beta.jac[which(!grepl("LIB02-TL01|LIB02-TL07",
                                               rownames(pred.beta.jac))),]

## Also remove these outlier samples from metadata
pred.metadata.no <- pred.metadata[which(!grepl("LIB02-TL01|LIB02-TL07",
                                               pred.metadata$Sequencing_ID)),]
rownames(pred.metadata.no) <- NULL

## Create occurrence and relative read abundance (i.e. proportional read 
## counts) dataframes
## Transform read counts to proportional read counts
pred.beta.bc <- decostand(pred.beta.jac, method="total", MARGIN=1)
pred.beta.bc.no <- decostand(pred.beta.jac.no, method="total", MARGIN=1)

## Convert read count data to presence-absence
pred.beta.jac <- decostand(pred.beta.jac, method="pa")
pred.beta.jac.no <- decostand(pred.beta.jac.no, method="pa")

## Tranform data into distance matrices
pred.beta.jac <- vegdist(pred.beta.jac, method="jaccard")
pred.beta.jac.no <- vegdist(pred.beta.jac.no, method="jaccard")
pred.beta.bc <- vegdist(pred.beta.bc, method="bray")
pred.beta.bc.no <- vegdist(pred.beta.bc.no, method="bray")

## Compute homogeneity of group dispersions (variances)
pred.bd.jac <- betadisper(pred.beta.jac, pred.metadata$Predator)
pred.bd.jac.no <- betadisper(pred.beta.jac.no, pred.metadata.no$Predator)
pred.bd.bc <- betadisper(pred.beta.bc, pred.metadata$Predator)
pred.bd.bc.no <- betadisper(pred.beta.bc.no, pred.metadata.no$Predator)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
mod.jac <- with(pred.metadata, pred.bd.jac)
mod.jac.no <- with(pred.metadata.no, pred.bd.jac.no)
mod.bc <- with(pred.metadata, pred.bd.bc)
mod.bc.no <- with(pred.metadata.no, pred.bd.bc.no)

mod.jac
mod.jac.no
mod.bc
mod.bc.no

## Compute mean distance to centroid and variance per group
tapply(pred.bd.jac$distances, pred.metadata$Predator, mean)
tapply(pred.bd.jac$distances, pred.metadata$Predator, var)

tapply(pred.bd.jac.no$distances, pred.metadata.no$Predator, mean)
tapply(pred.bd.jac.no$distances, pred.metadata.no$Predator, var)

tapply(pred.bd.bc$distances, pred.metadata$Predator, mean)
tapply(pred.bd.bc$distances, pred.metadata$Predator, var)

tapply(pred.bd.bc.no$distances, pred.metadata.no$Predator, mean)
tapply(pred.bd.bc.no$distances, pred.metadata.no$Predator, var)

## Ordination plot of distance to centroids
plot(pred.bd.jac)
plot(pred.bd.jac.no)
plot(pred.bd.bc)
plot(pred.bd.bc.no)

## Boxplot of distance to centroids
boxplot(pred.bd.jac, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

boxplot(pred.bd.jac.no, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

boxplot(pred.bd.bc, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

boxplot(pred.bd.bc.no, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

## Plots indicate that there is some difference in multivariate spread 
## between predators. Statistically check whether variance is different 
## between predators using standard parametric anova or permutation tests.
anova(pred.bd.jac)
anova(pred.bd.jac.no)
anova(pred.bd.bc)
anova(pred.bd.bc.no)

permutest(pred.bd.jac)
permutest(pred.bd.jac.no)
permutest(pred.bd.bc)
permutest(pred.bd.bc.no)

## Analyse pairwise differences between groups (predators) using 
## parametric Tukey's HSD test.
TukeyHSD(pred.bd.jac)  
TukeyHSD(pred.bd.jac.no) 
TukeyHSD(pred.bd.bc) 
TukeyHSD(pred.bd.bc.no) 

## Ordination of beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
pred.comm.jac <- metaMDS(pred.beta.jac, 
                         dist="jaccard", 
                         k=2,
                         maxit=999,
                         trymax=1000,
                         noshare=TRUE,
                         wascores=TRUE)

pred.comm.jac.no <- metaMDS(pred.beta.jac.no, 
                            dist="jaccard", 
                            k=2,
                            maxit=999,
                            trymax=1000,
                            noshare=TRUE,
                            wascores=TRUE)

pred.comm.bc <- metaMDS(pred.beta.bc, 
                        dist="bray", 
                        k=2,
                        maxit=999,
                        trymax=1000,
                        noshare=TRUE,
                        wascores=TRUE)

pred.comm.bc.no <- metaMDS(pred.beta.bc.no, 
                           dist="bray", 
                           k=2,
                           maxit=999,
                           trymax=1000,
                           noshare=TRUE,
                           wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
pred.comm.jac$stress
stressplot(pred.comm.jac)

pred.comm.jac.no$stress
stressplot(pred.comm.jac.no)

pred.comm.bc$stress
stressplot(pred.comm.bc)

pred.comm.bc.no$stress
stressplot(pred.comm.bc.no)

## Plot site scores as text
ordiplot(pred.comm.jac, display = "sites", type = "text", cex=0.5)
ordiplot(pred.comm.jac.no, display = "sites", type = "text", cex=0.5)
ordiplot(pred.comm.bc, display = "sites", type = "text", cex=0.5)
ordiplot(pred.comm.bc.no, display = "sites", type = "text", cex=0.5)

## Build dataframes with NMDS coordinates and metadata
pred.NMDS1 <- pred.comm.jac$points[,1]
pred.NMDS2 <- pred.comm.jac$points[,2]
pred.jac.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                            NMDS2=pred.NMDS2,
                            Predator = pred.metadata$Predator)

pred.NMDS1 <- pred.comm.jac.no$points[,1]
pred.NMDS2 <- pred.comm.jac.no$points[,2]
pred.jac.no.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                               NMDS2=pred.NMDS2,
                               Predator = pred.metadata.no$Predator)

pred.NMDS1 <- pred.comm.bc$points[,1]
pred.NMDS2 <- pred.comm.bc$points[,2]
pred.bc.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                           NMDS2=pred.NMDS2,
                           Predator = pred.metadata$Predator)

pred.NMDS1 <- pred.comm.bc.no$points[,1]
pred.NMDS2 <- pred.comm.bc.no$points[,2]
pred.bc.no.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                              NMDS2=pred.NMDS2,
                              Predator = pred.metadata.no$Predator)

## Check data
head(pred.jac.NMDS)
head(pred.jac.no.NMDS)
head(pred.bc.NMDS)
head(pred.bc.no.NMDS)

## Plot NMDS
p6a <- ggplot(pred.jac.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + stat_ellipse() + 
        scale_x_continuous(limits=c(-0.5,0.8), breaks=seq(-0.5,0.8,0.1),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.1,0.4), breaks=seq(-0.1,0.4,0.1),
                           labels=scales::number_format(accuracy = 0.1)) + 
        labs(title="A Jaccard dissimilarity", 
             subtitle="(i) All data", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.65, y=0.4, label="stress = 0.0002", cex=5) + 
        scale_colour_manual(name="Predator",
                                 values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=14),
              legend.position = "bottom",
              legend.direction="horizontal",
              legend.box="vertical",
              legend.key=element_blank())
p6a

p6b <- ggplot(pred.bc.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,0.8), breaks=seq(-0.6,0.8,0.1),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.3,0.3), breaks=seq(-0.3,0.3,0.15),
                                labels=scales::number_format(accuracy = 0.01)) + 
        labs(title="B Bray-Curtis dissimilarity", 
             subtitle="(i) All data", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.6, y=0.3, label="stress = 0.0003", cex=5) + 
        scale_colour_manual(name="Predator",
                            values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=14),
              legend.position = "bottom",
              legend.direction="horizontal",
              legend.box="vertical",
              legend.key=element_blank())
p6b

p6c <- ggplot(pred.jac.no.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,1), breaks=seq(-0.6,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.03,0.02), breaks=seq(-0.03,0.02,0.01),
                           labels=scales::number_format(accuracy = 0.01)) + 
        labs(title="", 
             subtitle="(ii) Outliers removed",
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.75, y=0.02, label="stress = 0.043", cex=5) + 
        scale_colour_manual(name="Predator",
                            values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=14),
              legend.position = "none")
p6c

p6d <- ggplot(pred.bc.no.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,1), breaks=seq(-0.6,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.02,0.04), breaks=seq(-0.02,0.04,0.02),
                           labels=scales::number_format(accuracy = 0.01)) + 
        labs(title="", 
             subtitle="(ii) Outliers removed", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.75, y=0.04, label="stress = 0.043", cex=5) + 
        scale_colour_manual(name="Predator",
                            values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=14),
              legend.position = "none")
p6d

## Show all plots together
g4 <- ggarrange(p6a,p6b,p6c,p6d,
                ncol=2, nrow=2,
                common.legend=TRUE, 
                legend="bottom")

#ggsave(filename="Figures/FigS5.png", 
#       plot = g4, width = 13, height = 10, dpi = 300, units = "in")

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
pred.jac.adonis <- adonis(pred.beta.jac ~ Predator, pred.metadata)
pred.jac.no.adonis <- adonis(pred.beta.jac.no ~ Predator, pred.metadata.no)
pred.bc.adonis <- adonis(pred.beta.bc ~ Predator, pred.metadata)
pred.bc.no.adonis <- adonis(pred.beta.bc.no ~ Predator, pred.metadata.no)

## Inspect results:
pred.jac.adonis
pred.jac.no.adonis 
pred.bc.adonis 
pred.bc.no.adonis 


#=======================================#
# BETA DIVERSITY: JACCARD DISSIMILARITY #
#=======================================#

## Create copy of diet.df
pred.beta.jac <- diet.df

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
pred.beta.jac <- pred.beta.jac[!sapply(pred.beta.jac, function(x) all(x == 0))]
pred.beta.jac <- pred.beta.jac[!apply(pred.beta.jac == 0, 1, all),]

## Remove samples that are extreme outliers: LIB02-TL01, LIB02-TL07
pred.beta.jac <- pred.beta.jac[which(!grepl("LIB02-TL01|LIB02-TL07", 
                                            rownames(pred.beta.jac))),]

## Also remove these outlier samples from metadata
pred.metadata <- pred.metadata[which(!grepl("LIB02-TL01|LIB02-TL07",
                                            pred.metadata$Sequencing_ID)),]
rownames(pred.metadata.no) <- NULL

## Convert read count data to presence-absence
pred.beta.jac[pred.beta.jac > 0] <- 1

## Create vectors of sample IDs associated with each predator
otter <- droplevels(subset(pred.metadata, select="Sequencing_ID", Predator == "Otter"))
mink <- droplevels(subset(pred.metadata, select="Sequencing_ID", Predator == "Mink"))

## Create copies of original dataframes for otter and mink
otter.beta <- pred.beta.jac[rownames(pred.beta.jac) %in% otter$Sequencing_ID,]
mink.beta <- pred.beta.jac[rownames(pred.beta.jac) %in% mink$Sequencing_ID,]

## Beta diversity across otter samples
otter.multi <- beta.multi(otter.beta, index.family="jaccard")
print(otter.multi)

## Beta diversity across mink samples
mink.multi <- beta.multi(mink.beta, index.family="jaccard")
print(mink.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for both otter and mink faecal samples.

## Pairwise between-site values of each component of beta diversity
pred.dist <- beta.pair(pred.beta.jac, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
pred.bd.turn <- betadisper(pred.dist$beta.jtu, pred.metadata$Predator)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
mod.turn <- with(pred.metadata, pred.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(pred.bd.turn$distances, pred.metadata$Predator, mean)

## Compute variance per group
tapply(pred.bd.turn$distances, pred.metadata$Predator, var)

## Ordination plot of distances to centroid
plot(pred.bd.turn)

## Boxplot of distances to centroid
boxplot(pred.bd.turn, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

## Plots indicate that there is some difference in multivariate 
## dispersion of turnover partition between predators. Statistically 
## check whether variance is different between predators using standard 
## parametric anova or permutation tests.
anova(pred.bd.turn)     # Significant difference between predators
permutest(pred.bd.turn) # Significant difference between predators

## Analyse pairwise differences between groups (predators) using 
## parametric Tukey's HSD test.
TukeyHSD(pred.bd.turn)  # Significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
pred.comm.turn <- metaMDS(pred.dist$beta.jtu, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          noshare=TRUE,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
pred.comm.turn$stress
stressplot(pred.comm.turn)

## Plot site scores as text
ordiplot(pred.comm.turn, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
pred.NMDS1 <- pred.comm.turn$points[,1]
pred.NMDS2 <- pred.comm.turn$points[,2]
pred.turn.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                             NMDS2=pred.NMDS2,
                             Predator = pred.metadata$Predator)

## Check data
head(pred.turn.NMDS)

## Plot data frame
p7a <- ggplot(pred.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,1), breaks=seq(-0.6,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.03,0.036), breaks=seq(-0.03,0.03,0.015),
                                labels=scales::number_format(accuracy = 0.001)) + 
        labs(title=expression(bold("C "~beta~"Diversity")), 
             subtitle="i  Turnover", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.85, y=0.03, label="stress = 0.043", cex=5) + 
        scale_colour_manual(name="Predator",
                            values=c("limegreen","purple")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position = "bottom",
              legend.key=element_blank(),
              legend.title=element_blank())
p7a

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
pred.turn.adonis <- adonis(pred.dist$beta.jtu ~ Predator, pred.metadata)

## Inspect results:
pred.turn.adonis

## Result is significant. There is a substantial difference in taxon 
## replacement (i.e. turnover) between predator diet. Therefore, taxa
## consumed by one predator are substituted by species for another 
## predator.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
pred.bd.nest <- betadisper(pred.dist$beta.jne, pred.metadata$Predator)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
pred.nest <- with(pred.metadata, pred.bd.nest)
pred.nest

## Compute mean distance to centroid per group
tapply(pred.bd.nest$distances, pred.metadata$Predator, mean)

## Compute variance per group
tapply(pred.bd.nest$distances, pred.metadata$Predator, var)

## Ordination plot of distances to centroid
plot(pred.bd.nest)

## Boxplot of distances to centroid
boxplot(pred.bd.nest, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

## Plots indicates that there is no difference in multivariate dispersion
## between predators. Statistically check whether variance is different 
## between predators using standard parametric anova or permutation tests.
anova(pred.bd.nest)     # No significant difference between predators
permutest(pred.bd.nest) # No significant difference between predators

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(pred.bd.nest)  # No significant difference between predators

## Ordination of beta diversity partitioned by nestedness:
pred.comm.nest <- metaMDS(pred.dist$beta.jne, 
                          dist="jaccard",
                          k=2,
                          maxit=999,
                          trymax=1000,
                          noshare=TRUE,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
pred.comm.nest$stress
stressplot(pred.comm.nest)

## plot site scores as text
ordiplot(pred.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
pred.NMDS1 <- pred.comm.nest$points[,1]
pred.NMDS2 <- pred.comm.nest$points[,2]
pred.nest.NMDS <- data.frame(NMDS1=pred.NMDS1, 
                             NMDS2=pred.NMDS2,
                             Predator = pred.metadata$Predator)

## Check data
head(pred.nest.NMDS)

## Plot data frame
p7b <- ggplot(pred.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.3),
                           labels=scales::number_format(accuracy = 0.1)) + 
        labs(title="", 
             subtitle="ii  Nestedness-resultant",
             x="NMDS1", y="NMDS2") + 
        annotate("text", x=0.5, y=0.6, label="stress = 0.176", cex=5) + 
        scale_colour_manual(name="diet set",
                            values=c("limegreen","purple")) + 
        guides(colour=guide_legend(nrow=2, byrow=TRUE)) + 
        theme(panel.background = element_rect(fill = 'white'), 
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position="none")
p7b

## Statistically check difference in nestedness of communities
## Look at PermANOVA using adonis(), considered more robust than anosim()
pred.nest.adonis <- adonis(pred.dist$beta.jne ~ Predator, pred.metadata)

## Inspect results
## no summary() or plot() diets included
pred.nest.adonis

## Result is not significant. There is no substantial difference in taxon
## loss or gain (i.e. nestedness) between predators across sites.


## 3. TOTAL BETA DIVERSITY
pred.bd.total <- betadisper(pred.dist$beta.jac, pred.metadata$Predator)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(pred.metadata, pred.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(pred.bd.total$distances, pred.metadata$Predator, mean)

## Compute variance per group
tapply(pred.bd.total$distances, pred.metadata$Predator, var)

## Ordination plot of distance to centroids
plot(pred.bd.total)

## Boxplot of distance to centroids
boxplot(pred.bd.total, xlab="Predator", xaxt="n", bty="n")
axis(side=1, at=c(1:2), labels=c("Mink","Otter"))

## Plots indicates that there is some difference in multivariate dispersion 
## between predators. Statistically check whether variance is different 
## between predators using standard parametric anova or permutation tests.
anova(pred.bd.total)     # Significant difference between predators
permutest(pred.bd.total) # Significant difference between predators

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(pred.bd.total)  # Significant difference between predators

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
pred.comm.total <- metaMDS(pred.dist$beta.jac, 
                           dist="jaccard",
                           k=2,
                           maxit=999,
                           trymax=1000,
                           noshare=TRUE,
                           wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
pred.comm.total$stress
stressplot(pred.comm.total)

## plot site scores as text
ordiplot(pred.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
pred.NMDS1 <- pred.comm.total$points[,1]
pred.NMDS2 <- pred.comm.total$points[,2]
pred.total.NMDS <- data.frame(NMDS1=pred.NMDS1,
                              NMDS2=pred.NMDS2,
                              Predator=pred.metadata$Predator)

## Check data
head(pred.total.NMDS)

## Plot data frame
p7c <- ggplot(pred.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Predator)) + 
        geom_point(cex=2, alpha=0.3) + stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,1), breaks=seq(-0.6,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.02,0.03), breaks=seq(-0.02,0.03,0.01),
                           labels=scales::number_format(accuracy = 0.001)) +
        labs(title="",
             subtitle=expression(bold("iii  Total"~beta~"Diversity")),
             x="NMDS1", y="NMDS2") + 
        annotate("text", x=0.85, y=0.03, label="stress = 0.043", cex=5) + 
        scale_colour_manual(name="Predator",
                            values=c("limegreen","purple")) + 
        guides(colour=guide_legend(nrow=2, byrow=TRUE)) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position="none")
p7c

## Statistically check difference in total beta diversity of ponds
## Look at PermANOVA using adonis(), considered more robust than anosim()
pred.total.adonis <- adonis(pred.dist$beta.jac ~ Predator, pred.metadata)

## Inspect results
## no summary() or plot() diets included
pred.total.adonis

## Again result is significant. There is substantial variation in overall 
## prey composition of samples belong to different predators.

## Plot all diversity results
g5 <- ggarrange(p4,
                ggarrange(p5a,p5b,p5c, nrow=3, ncol=1, align="hv",
                          common.legend=TRUE, legend="bottom"), 
                ggarrange(p7a,p7b,p7c, nrow=3, ncol=1, align="hv",
                          common.legend=TRUE, legend="bottom"), 
                nrow=1, ncol=3)

#ggsave(filename="Figures/Fig5.png", 
#       plot = g5, width = 23, height = 15, dpi = 300, units = "in")



#########################
# OTTER SITE COMPARISON #
#########################

#=================#
# ALPHA DIVERSITY #
#=================#

## Subset alpha dataframe for otter samples
otter.alpha <- droplevels(subset(alpha, Predator=="Otter"))

## Kruskal-Wallis test:
kruskal.test(sample.richness ~ Site, data=otter.alpha)

## Dunn's test of multiple comparisons
dunnTest(sample.richness ~ Site,
         data=otter.alpha,
         method="bh")

## Plot prey richness in otter samples according to site and waterbody
p8 <- ggplot(otter.alpha, aes(x=Site, y=sample.richness, colour=Site)) + 
        geom_jitter(cex=2, 
                    width=0.15,
                    height=0.2) + 
        geom_boxplot(colour="black",
                     alpha=0.7,
                     outlier.shape=NA) + 
        annotate("text", x = c("Malham Tarn",
                               "River Glaven",
                               "River Hull"), y = 10, 
                 label = c("a","a","b"), cex=10) + 
        scale_colour_manual(values=c("grey30","goldenrod2","dodgerblue3")) + 
        scale_y_continuous(limits=c(0,10),
                           breaks=seq(0,10,1)) + 
        labs(title=expression(bold("A"~alpha~"Diversity")), 
             x="Site", y="Taxon richness") + 
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(size=25, margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(size=25, margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", size=20),
              axis.text.y = element_text(colour="black", size=24),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "none",
              text = element_text(size=20))
p8


#==================================#
# RAREFACTION/EXTRAPOLATION CURVES #
#==================================#

## Look at whether differences in richness are influenced by sampling
## effort using iNEXT package. Create vectors of sample IDs associated 
## with each site.
MT <- droplevels(subset(otter.alpha, select="Sequencing_ID", Site == "Malham Tarn"))
RG <- droplevels(subset(otter.alpha, select="Sequencing_ID", Site == "River Glaven"))
RH <- droplevels(subset(otter.alpha, select="Sequencing_ID", Site == "River Hull"))

## Create copies of pred.pa for each site
MT.pa <- pred.pa[rownames(pred.pa) %in% MT$Sequencing_ID,]
RG.pa <- pred.pa[rownames(pred.pa) %in% RG$Sequencing_ID,]
RH.pa <- pred.pa[rownames(pred.pa) %in% RH$Sequencing_ID,]

## Calculate incidence frequency for each species detected in otter 
## samples from each site
MT.richness <- colSums(MT.pa)
RG.richness <- colSums(RG.pa)
RH.richness <- colSums(RH.pa)

## Add sample size to beginning of each richness vector. This is because 
## the first entry of each list for iNEXT must be the total number of 
## sampling units, followed by the species incidence frequencies.
MT.richness <- append(MT.richness, 24, after=0)
RG.richness <- append(RG.richness, 35, after=0)
RH.richness <- append(RH.richness, 110, after=0)

## Make list of otter and mink samples
site.richness <- list(MT.richness, RG.richness, RH.richness)
names(site.richness) <- c("Malham Tarn", "River Glaven", "River Hull")

## Run iNEXT function:
## With 300 samples:
site.re.300 <- iNEXT(site.richness, q=0, datatype="incidence_freq", 
                     endpoint=300, knots=60, se=TRUE, conf=0.95,
                     nboot=1000)
site.re.300

## With 3500 samples:
site.re.3500 <- iNEXT(site.richness, q=0, datatype="incidence_freq", 
                     endpoint=3500, knots=700, se=TRUE, conf=0.95,
                     nboot=1000)
site.re.3500

## Sample-size-based R/E curves
p9a <- ggiNEXT(site.re.300, type=1) + 
        scale_colour_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits=c(0,300), breaks=seq(0,300,50)) +
        scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
        labs(title="B  Rarefaction/Extrapolation (R/E) curves", 
             subtitle="i  Sample size-based R/E curve",
             y="Taxon diversity") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="bottom",
              legend.box="vertical",
              legend.key=element_blank(),
              text=element_text(size=20))
p9a

## Sample completeness curves
p9b <- ggiNEXT(site.re.300, type=2) +
        scale_colour_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits=c(0,300), breaks=seq(0,300,50)) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
        labs(title="",
             subtitle="ii  Sample completeness curve") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="none",
              text=element_text(size=20))
p9b

## Coverage-based R/E curves
p9c <- ggiNEXT(site.re.300, type=3) + 
        scale_colour_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
        scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10)) +
        labs(title="",
             subtitle="iii  Coverage-based R/E curve",
             y="Taxon diversity") +
        theme(panel.background=element_rect(fill="white"),
              panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(),
              axis.line.x=element_line(colour="black", size=0.5, linetype="solid"),
              axis.line.y=element_line(colour="black", size=0.5, linetype="solid"),
              axis.title.x=element_text(margin=unit(c(8, 0, 0, 0), "mm")),
              axis.title.y=element_text(margin=unit(c(0, 5, 0, 0), "mm")),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black"),
              plot.title=element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle=element_text(face="bold", hjust=0, colour="black", margin=unit(c(2, 0, 0, 0), "mm")),
              legend.position="none",
              text=element_text(size=20))
p9c

## Apply the estimateD() function to obtain diversity estimates of order 
## q = 0, 1, 2 for any particular level of sample size (base="size") or
## any specified level of sample coverage (base="coverage")
estimateD(site.richness, datatype="incidence_freq", base="size",
          level=NULL, conf=0.95)
estimateD(site.richness, datatype="incidence_freq", base="coverage", 
          level=0.95, conf=0.95)



#######################################
# BETA DIVERSITY BY SITE (OTTER ONLY) #
#######################################

## Create copy of diet.df
otter.beta.jac <- diet.df

## Subset data for otter samples
otter.beta.jac <- otter.beta.jac[rownames(otter.beta.jac) %in% otter$Sequencing_ID,]

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
otter.beta.jac <- otter.beta.jac[!sapply(otter.beta.jac, function(x) all(x == 0))]
otter.beta.jac <- otter.beta.jac[!apply(otter.beta.jac == 0, 1, all),]

## Remove samples that are extreme outliers: LIB02-TL07, LIB04-TL57
otter.beta.jac <- otter.beta.jac[which(!grepl("LIB02-TL07|LIB04-TL57",
                                              rownames(otter.beta.jac))),]

## Remove extreme outliers from metadata
otter.metadata <- droplevels(subset(pred.metadata, Predator == "Otter"))
otter.metadata <- otter.metadata[which(!grepl("LIB02-TL07|LIB04-TL57",
                                              otter.metadata$Sequencing_ID)),]
rownames(otter.metadata) <- NULL

## Convert read count data to presence-absence
otter.beta.jac[otter.beta.jac > 0] <- 1

## Subset otter dataframe for samples from each site
MT.beta <- otter.beta.jac[rownames(otter.beta.jac) %in% MT$Sequencing_ID,]
RG.beta <- otter.beta.jac[rownames(otter.beta.jac) %in% RG$Sequencing_ID,]
RH.beta <- otter.beta.jac[rownames(otter.beta.jac) %in% RH$Sequencing_ID,]

## Beta diversity across Malham Tarn samples
MT.multi <- beta.multi(MT.beta, index.family="jaccard")
print(MT.multi)

## Beta diversity across River Glaven samples
RG.multi <- beta.multi(RG.beta, index.family="jaccard")
print(RG.multi)

## Beta diversity across River Hull samples
RH.multi <- beta.multi(RH.beta, index.family="jaccard")
print(RH.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for otter samples from all sites.

## Pairwise between-site values of each component of beta diversity
site.dist <- beta.pair(otter.beta.jac, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
site.bd.turn <- betadisper(site.dist$beta.jtu, otter.metadata$Site)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
site.turn <- with(otter.metadata, site.bd.turn)
site.turn

## Compute mean distance to centroid per group
tapply(site.bd.turn$distances, otter.metadata$Site, mean)

## Compute variance per group
tapply(site.bd.turn$distances, otter.metadata$Site, var)

## Ordination plot of distances to centroid
plot(site.bd.turn)

## Boxplot of distances to centroid
boxplot(site.bd.turn, xlab="Site", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("Malham Tarn","River Glaven","River Hull"))

## Plots indicate that there is some difference in multivariate dispersions
## between sites. Statistically check whether turnover is different 
## between sites using standard parametric anova or permutation tests.
anova(site.bd.turn)     # Significant difference between sites
permutest(site.bd.turn) # Significant difference between sites

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(site.bd.turn)  # Significant difference between sites

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
site.comm.turn <- metaMDS(site.dist$beta.jtu, 
                          distance="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          noshare=TRUE,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
site.comm.turn$stress
stressplot(site.comm.turn)

## Plot site scores as text
ordiplot(site.comm.turn, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
site.NMDS1 <- site.comm.turn$points[,1]
site.NMDS2 <- site.comm.turn$points[,2]
site.turn.NMDS <- data.frame(NMDS1=site.NMDS1, 
                             NMDS2=site.NMDS2,
                             Site = otter.metadata$Site)

## Check data
head(site.turn.NMDS)

## Plot data frame
p10a <- ggplot(site.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Site)) + 
        geom_point(cex=2, alpha=0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.2,1), breaks=seq(-0.2,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.02,0.03), breaks=seq(-0.02,0.03,0.01),
                           labels=scales::number_format(accuracy = 0.01)) + 
        labs(title=expression(bold("C "~beta~"Diversity")),
             subtitle="i  Turnover", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.85, y=0.03, label="stress = 0.077", cex=5) + 
        scale_colour_manual(name="Site",
                            values=c("grey30","goldenrod2","dodgerblue3")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position = "bottom",
              legend.key=element_blank(),
              legend.title=element_blank())
p10a

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
site.turn.adonis <- adonis(site.dist$beta.jtu ~ Site, otter.metadata)

## Inspect results:
site.turn.adonis

## Result is significant. There is a substantial difference in species 
## replacement (i.e. turnover) between sites. Therefore, species 
## consumed at one site are substituted by species at a different site.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
site.bd.nest <- betadisper(site.dist$beta.jne, otter.metadata$Site)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
site.nest <- with(otter.metadata, site.bd.nest)
site.nest

## Compute mean distance to centroid per group
tapply(site.bd.nest$distances, otter.metadata$Site, mean)

## Compute variance per group
tapply(site.bd.nest$distances, otter.metadata$Site, var)

## Ordination plot of distances to centroid
plot(site.bd.nest)

## Boxplot of distances to centroid
boxplot(site.bd.nest, xlab="Site", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("Malham Tarn","River Glaven","River Hull"))

## Plots indicate that there is some difference in multivariate dispersions
## between sites. Statistically check whether nestedness is different 
## between sites using standard parametric anova or permutation tests.
anova(site.bd.nest)     # Significant difference between sites
permutest(site.bd.nest) # Significant difference between sites

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(site.bd.nest)  # Significant difference between sites

## Ordination of beta diversity partitioned by nestedness:
## The metaMDS function automatically transforms data and checks solution
## robustness
site.comm.nest <- metaMDS(site.dist$beta.jne, 
                          distance="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          noshare=TRUE,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
site.comm.nest$stress
stressplot(site.comm.nest)

## Plot site scores as text
ordiplot(site.comm.nest, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
site.NMDS1 <- site.comm.nest$points[,1]
site.NMDS2 <- site.comm.nest$points[,2]
site.nest.NMDS <- data.frame(NMDS1=site.NMDS1, 
                             NMDS2=site.NMDS2,
                             Site = otter.metadata$Site)

## Check data
head(site.nest.NMDS)

## Plot data frame
p10b <- ggplot(site.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Site)) + 
        geom_point(cex=2, alpha=0.3) + stat_ellipse() + 
        scale_x_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.3),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.3),
                           labels=scales::number_format(accuracy = 0.1)) + 
        labs(title="", 
             subtitle="ii  Nestedness-resultant", 
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.45, y=0.6, label="stress = 0.183", cex=5) + 
        scale_colour_manual(name="Site",
                            values=c("grey30","goldenrod2","dodgerblue3")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position = "none")
p10b

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
site.nest.adonis <- adonis(site.dist$beta.jne ~ Site, otter.metadata)

## Inspect results:
site.nest.adonis

## Result is not significant. There is no difference in species loss or 
## gain (i.e. nestedness) between sites.


## 3. TOTAL BETA DIVERSITY
site.bd.total <- betadisper(site.dist$beta.jac, otter.metadata$Site)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
site.total <- with(otter.metadata, site.bd.total)
site.total

## Compute mean distance to centroid per group
tapply(site.bd.total$distances, otter.metadata$Site, mean)

## Compute variance per group
tapply(site.bd.total$distances, otter.metadata$Site, var)

## Ordination plot of distances to centroid
plot(site.bd.total)

## Boxplot of distances to centroid
boxplot(site.bd.total, xlab="Site", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("Malham Tarn","River Glaven","River Hull"))

## Plots indicate that there is some difference in multivariate dispersions 
## between sites. Statistically check whether variance is different 
## between sites using standard parametric anova or permutation tests.
anova(site.bd.total)     # Significant difference between sites
permutest(site.bd.total) # Significant difference between sites

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(site.bd.total)  # Significant difference between sites

## Ordination of beta diversity partitioned by total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
site.comm.total <- metaMDS(site.dist$beta.jac, 
                           distance="jaccard", 
                           k=2,
                           maxit=999,
                           trymax=1000,
                           noshare=TRUE,
                           wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
site.comm.total$stress
stressplot(site.comm.total)

## Plot site scores as text
ordiplot(site.comm.total, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
site.NMDS1 <- site.comm.total$points[,1]
site.NMDS2 <- site.comm.total$points[,2]
site.total.NMDS <- data.frame(NMDS1=site.NMDS1, 
                             NMDS2=site.NMDS2,
                             Site = otter.metadata$Site)

## Check data
head(site.total.NMDS)

## Plot data frame
p10c <- ggplot(site.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Site)) + 
        geom_point(cex=2, alpha=0.3) + stat_ellipse() + 
        scale_x_continuous(limits=c(-0.2,1), breaks=seq(-0.2,1,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.01,0.02), breaks=seq(-0.01,0.02,0.01), 
                           labels=scales::number_format(accuracy = 0.01)) + 
        labs(title="", 
             subtitle=expression(bold("iii  Total"~beta~"Diversity")),
             x="NMDS1",y="NMDS2") + 
        annotate("text", x=0.85, y=0.02, label="stress = 0.077", cex=5) + 
        scale_colour_manual(name="Site",
                            values=c("grey30","goldenrod2","dodgerblue3")) + 
        theme(panel.background = element_rect(fill = 'white'),
              axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_text(colour="black"),
              plot.title = element_text(face="bold", hjust=0, colour="black"),
              plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size=20),
              legend.position = "none")
p10c

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
site.total.adonis <- adonis(site.dist$beta.jac ~ Site, otter.metadata)

## Inspect results:
site.total.adonis

## Again result is significant. There is substantial variation in overall 
## prey composition of otter samples from different sites.

## Plot all diversity results
g6 <- ggarrange(p8,
                ggarrange(p9a,p9b,p9c, nrow=3, ncol=1, align="hv",
                          common.legend=TRUE, legend="bottom"), 
                ggarrange(p10a,p10b,p10c, nrow=3, ncol=1, align="hv",
                          common.legend=TRUE, legend="bottom"), 
                nrow=1, ncol=3)

#ggsave(filename="Figures/Fig6.png", 
#       plot = g6, width = 23, height = 15, dpi = 300, units = "in")


#################
# END OF SCRIPT #
#################

