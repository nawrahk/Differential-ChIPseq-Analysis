---
#title: "Differential ChIP-seq/ATAC-seq Analysis"
#author: "Nawrah Khader"
---
  
###Installation of required packages

#install.packages('tibble')
#install.packages('tidyverse')
#install.packages('tidyr')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install('DiffBind')
#BiocManager::install('DESeq2')
#BiocManager::install('RColorBrewer')
#BiocManager::install('ComplexHeatmap')
#BiocManager::install('circlize')
#BiocManager::install('digest')
#BiocManager::install('cluster')




###Load packages 
library(tibble)
library(tidyverse)
library(DiffBind)
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)




args <- commandArgs(TRUE)

###These correspond to the bash command above:

#argument 1 refers to the sampledata sheet that is required to import your samples and the information pertaining      to them
samplefile <- args[1]

#argument 2 is the name for your output
name <- args [2]




#Import the files associated with your samples (bam and bed [or narrowpeak])
samples <- read.csv(file = samplefile)
chip <- dba(sampleSheet = samples)

#Confirms that the files have been correctly imported and prints plot as a pdf
pdf(file = paste('chip_samples_plot', name='.pdf', sep=''))
plot(chip)
dev.off()




#Overlapping venn diagram analysis for consensus peaks, and the consenesus argument will vary based on how you set up your sampledata.csv. If your samples contain different treatments, it will be DBA_TREATMENT and if it has different conditions then it's  DBA_CONDITION
chip_consensus <- dba.peakset(chip, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=0.66)
chip_consensus <- dba(chip_consensus, mask=chip_consensus$masks$Consensus, minOverlap=1)

#Print the consensus peaks to a tab delimited file
consensus_peaks <- dba.peakset(chip_consensus, bRetrieve=TRUE, DataType = DBA_DATA_FRAME, writeFile = "chip_consensuspeaks.tsv")

#Construct and print the Venn diagram plot as a pdf. Selecting the mask as consensus ensures that you select the peaks in the DBA object that are consensus peaks in your samples.
pdf(file = paste('chip_consensus_vennplot', name='.pdf', sep=''))
chip_venn <- dba.plotVenn (chip_consensus, mask=chip_consensus$masks$Consensus)
dev.off()




#Count features from the dataset, additional parameters added to reduce memory usage
chipcounts <- dba.count(chip, bUseSummarizeOverlaps=TRUE, bParallel = TRUE)

#Normalize counts based on sequencing depth
chipcounts_norm <- dba.normalize(chipcounts)

#Visualize and print the correlation heatmap and principal component analysis (PCA) plot as pdfs to see how your samples and their respective replicates cluster.
pdf(file=paste('normalized_chip_counts_plots', name='.pdf', sep=''))
dba.plotHeatmap(chipcounts_norm)

#For the PCA plot, as before, the attributes will vary based on how you set up your sampledata.csv. 
dba.plotPCA(chipcounts_norm, attributes = DBA_TISSUE, label = DBA_ID)

dev.off()


#Before running the main differential analysis, you need to specify the groups (DBA_CONDTION or DBA_TREATMENT) you are comparing. This again relies on how your sampledata.csv is set up.
chipcontrast <- dba.contrast(chipcounts_norm, categories = c(DBA_TISSUE, DBA_CONDITION), minMembers = 2)


#This is the main function of DiffBind, which enables the identification of statistical significant differential binding sites between  sample groups. By default, it uses DESeq2, but you can specify a different method by writing "method=DBA_ALL_METHODS" or "method=DBA_EDGER"
chipdiff <- dba.analyze(chipcontrast)


#Generate and export report containing results from the differential analysis.
chipdiff_report <- dba.report(chipdiff)

#Arrange in order of decreasing fold change
chipdiff_report <- as.data.frame(chipdiff_report) %>%
  arrange(desc(Fold))

#Export
chipdiff_report_df <- dba.peakset(chipdiff_report, DataType=DBA_DATA_FRAME, writeFile = paste('DiffBindReport_', name, '.csv', sep=""))


###Customized volcano plot using Enhanced Volcano based on the differential analysis and specified cut-offs. Peaks meeting specified cut-offs are colored.

#specify fold change cut-off
FC <- 1.5
#specify adjusted p value cut-off
p <- 0.01

#Peaks that don't meet fold change or adjusted pvalue (FDR) cut-offs
keyvals <- rep('grey75', nrow(chipdiff_report))
names(keyvals) <- rep('NS', nrow(chipdiff_report))

#Peaks that don't meet adjusted pvalue (FDR) cut-off but meet fold change cut-off
keyvals[which(abs(chipdiff_report$Fold) > FC & chipdiff_report$FDR > p)] <- 'grey50'
names(keyvals)[which(abs(chipdiff_report$Fold) > FC & chipdiff_report$FDR > p)] <- 'log2FoldChange'

#Peaks that don't meet fold change cut-off but meet adjusted pvalue (FDR) cut-off 
keyvals[which(abs(chipdiff_report$Fold) < FC & chipdiff_report$FDR < p)] <- 'grey25'
names(keyvals)[which(abs(chipdiff_report$Fold)  < FC & chipdiff_report$FDR < p)] <- 'adjusted p-value'

#Peaks that meet negative fold change cut-off and meet adjusted pvalue (FDR) cut-off 
keyvals[which(chipdiff_report$Fold < -FC & chipdiff_report$FDR < p)] <- 'darkorange'
names(keyvals)[which(chipdiff_report$Fold  < -FC & chipdiff_report$FDR < p)] <- 'Signif. depleted'

#Peaks that meet positive fold change cut-off and meet adjusted pvalue (FDR) cut-off 
keyvals[which(chipdiff_report$Fold > FC & chipdiff_report$FDR < p)] <- 'red'
names(keyvals)[which(chipdiff_report$Fold > FC & chipdiff_report$FDR < p)] <- 'Signif. enriched'

#Ensures unique values
unique(keyvals)
unique(names(keyvals))

#Plotting the volcano plot with specific parameters. This can be adjusted as per parameters on the EnhancedVolcano manual.
library(EnhancedVolcano)
p <- EnhancedVolcano(chipdiff_report,
                     lab = NA,
                     x = "Fold",
                     y = "FDR",
                     title = 'ChIPseq Diffbind',
                     pCutoff = 0.01,
                     FCcutoff = 1.0,
                     xlim = c(-10,10),
                     pointSize = 1,
                     colCustom=keyvals)

#Export volcano plot
pdf(file = paste('DiffBind_ChIPseq_EnhancedVolcano', name='.pdf', sep=''))
p
dev.off()




###The following allows you to extract the normalized read counts (RPKM) from the differential analysis (DBA) object we created above. 

#Creating the counts dataframe and selecting the first three columns gives you the chr, start, end (peaks).

peak_df <- chipdiff$peaks[[1]] %>% select(1:3)

for (x in 1:5) {
  current_peak_set <- chipdiff$peaks[[x]] %>% select(5)
  #naming the counts columns based on the names of your samples
  timepoint_name <- colnames(chipdiff[["class"]])[[x]]
  #select the counts (RPKM) in the dba object
  peak_df[timepoint_name] <- current_peak_set$RPKM
}

#Columns 4-9 in this dataframe correspond to the counts columns per replicate. For this example, there were two conditions, with three replicates each which is why it extends to 9 [3 + (2 conditions x 3 reps = 6) = 9]. This number will change for more/less replicates and conditions.
##The following command selects the counts columns from the above dataframe, scales the counts, and writes them into the original dataframe. 
peak_df <- peak_df %>% 
  mutate_at(4:9, funs(c(scale(.))))


#Write the dataframe consisting of normalized read counts for each sample (and their replicates) at all peaks
write_csv(peak_df, "chipseq_peak_RPKM.csv")


#Select the significant peaks that meet the fold change and adjusted pvalue cut-off
chip_report_filtered <- chipdiff_report %>% 
  dplyr::filter (FDR < 0.01 & (Fold >= 1.5 | Fold <= -1.5))

#Select only the peak coordinates and combine the chr, start, end columns into one column
chip_report_filtered_peaks <- chip_report_filtered %>%
  dplyr::select(1:3) %>%
  unite(col = 'peak', c('Chr', 'Start', 'End'), sep='.')


#Tidy the dataframe consisting of counts so the chr, start, end columns are combined into one column
peak_df_new <- peak_df %>%
  unite(col = 'peak', c('seqnames', 'start', 'end'), sep='.')

###Creating the heatmap using the counts from the Diffbind analysis and subsetting by the significantly differential peaks obtained. This is largely based off of the "Simple Tutorial for a ComplexHeatmap by Kevin Blighe"

#Subset the counts dataframe by significant peaks
chip_matrix <- peak_df_new[sig_peaks, ]

#Filter any rows where all read values are 0
chip_matrix <- dplyr::filter_all(chip_matrix, any_vars(. != 0))

#Scale the count matrix
chip_scaled <- t(scale(t(chip_matrix)))

#Remove any NA values
chip_no_NA <- na.omit(chip_scaled)


#Perform clustering analysis using pamclusters *N.B.only 65536 rows are allowed, if you exceed then, then you have to try a different clustering method.

#Select the number of k clusters
pamClusters <- cluster::pam(chip_no_NA, k = 3)
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)

#Fix order of the clusters to have 1 to 4, top to bottom. The order of the clusters can be changed
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 3', 'Cluster 3'))



#Choose colors for your heatmap 
myCol <- colorRampPalette(c('#5b9ad5', 'white', '#ed792c'))(100)
myBreaks <- seq(-3, 3, length.out = 100)


##Constructing the heatmap
chip_heatmap <- Heatmap(chip_no_NA,
                        
                        #split the counts rows according to the PAM clusters
                        split = pamClusters$clustering,
                        cluster_row_slices = T,
                        
                        name = 'Z-\nscore',
                        
                        col = colorRamp2(myBreaks, myCol),
                        
                        # parameters for the colour-bar that represents gradient of enrichment
                        heatmap_legend_param = list(
                          color_bar = 'continuous',
                          legend_direction = 'vertical',
                          legend_width = unit(8, 'cm'),
                          legend_height = unit(5.0, 'cm'),
                          title_position = 'topcenter',
                          title_gp=gpar(fontsize = 12, fontface = 'bold'),
                          labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                        
                        # row (peak) parameters
                        cluster_rows = TRUE,
                        show_row_dend = TRUE,
                        row_title = 'Statistically significant peaks',
                        row_title_side = 'left',
                        row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                        row_title_rot = 90,
                        show_row_names = F,
                        row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                        row_names_side = 'left',
                        row_dend_width = unit(25,'mm'),
                        
                        # column (sample) parameters
                        cluster_columns = T,
                        column_split = 4,
                        show_column_dend = TRUE,
                        column_title = '',
                        column_title_side = 'bottom',
                        column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                        column_title_rot = 0,
                        show_column_names = TRUE,
                        column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                        column_names_max_height = unit(10, 'cm'),
                        column_dend_height = unit(25,'mm'),
                        
                        # cluster methods for rows and columns
                        clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                        clustering_method_columns = 'ward.D2',
                        clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                        clustering_method_rows = 'ward.D2')

chip_heatmap 


#Export heatmap as a PDF

pdf(file = paste('chip_differentialpeaks_heatmap', name='.pdf', sep=''))
draw(chip_heatmap)
dev.off()


##The following allows you to obtain the peaks (and the counts for all replicates) that are associated with each cluster

#Obtaining the cluster number
chip_heatmap_clus <- as.data.frame(chip_heatmap@matrix_param[["row_split"]]) %>% 
  rownames_to_column(var = "Row")

#Obtaining the peak location and scaled counts for all replicates
chip_heatmap_counts <- as.data.frame(chip_heatmap@matrix) %>% 
  rownames_to_column(var = "Row")

#Merging the dataframe that has cluster numbers with the dataframe containing peak location and scaled counts for all replicates 
chip_heatmap_countswclus <- merge(hmap_counts, hmap_clus, by = "Row") %>% 
  dplyr::select (-Row)

write_csv(chip_heatmap_countswclus, "chipseq_peak_scaledcounts_clus.csv")

