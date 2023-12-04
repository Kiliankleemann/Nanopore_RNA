BiocManager::install('rtracklayer')
######--------LOADING LIBRARIES--------########
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(limma)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(KEGGREST)
library(RColorBrewer)
library(ggrepel)
library(fgsea)
library(clusterProfiler)
library(pheatmap)
library(BiocGenerics)
library(tximport)
library(S4Vectors)
library(DESeq2)
library(openxlsx)
library(biomaRt)
library(tidyverse)
library(ggplot2)
library(edgeR)
library(clusterProfiler)
library(ggpubr)
require(reshape)
library(splitstackshape)





##---------------------------------FILE PREPARATION TRANSCRIPTS----------------------------------##
# ## useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
# ensembl_hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# ## getBM is a function used to retrieve information from the BioMart database
# ensembl_df <- getBM(attributes = c("refseq_mrna", "external_gene_name"), mart = ensembl_hs_mart) %>%
#   as_tibble()
# 
# # Build tx2gene table from ensembl_df
# tx2gene <- ensembl_df %>%
#   dplyr::select(refseq_mrna, external_gene_name) %>%
#   rename(c('transcript_id' = "refseq_mrna" , "Gene" = "external_gene_name" )) %>%
#   data.frame()

##### ---------------------------------IMPORT DATA TRANSCRIPTS---------------------------------- #####
cDNA_transcript_count <- as.data.frame(read.csv(paste0("/media/kilian/OS/Melani_THP1_August_2023/cDNA_ANALYSIS/stringtie_output_cDNA_specific_L1/transcript_count_matrix.csv"), sep = ","))
# cDNA_transcript_count_edit <- inner_join(tx2gene,cDNA_transcript_count,'transcript_id')
# cDNA_transcript_count_edit$Gene_transcript <- paste0(cDNA_transcript_count_edit$Gene,"_",RNA_transcript_count_edit$transcript_id)

##### ---------------------------------NORMALIZATION OF COUNTS---------------------------------- #####
cDNA_transcript_count$S1_SUM <- sum(cDNA_transcript_count$S1_cDNA)
cDNA_transcript_count$WT_SUM <- sum(cDNA_transcript_count$WT_cDNA)

cDNA_transcript_count$S1_norm_counts <- (cDNA_transcript_count$S1_cDNA)/(cDNA_transcript_count$S1_SUM)*1000000
cDNA_transcript_count$WT_norm_counts <- (cDNA_transcript_count$WT_cDNA)/(cDNA_transcript_count$WT_SUM)*1000000



##### ---------------------------------CALCULATE LOG2 FC & FILTER ---------------------------------- #####
cDNA_transcript_count$SUM <- cDNA_transcript_count$S1_cDNA + cDNA_transcript_count$WT_cDNA
cDNA_transcript_count_filter <- cDNA_transcript_count %>% filter(SUM > 100)
cDNA_transcript_count_filter$FC <- (cDNA_transcript_count_filter$S1_norm_counts+1) /(cDNA_transcript_count_filter$WT_norm_counts+1)
cDNA_transcript_count_filter$LOG2_FC <- log(cDNA_transcript_count_filter$FC,2)
cDNA_transcript_count_filter_top <- cDNA_transcript_count_filter %>% filter(LOG2_FC > 2 | LOG2_FC < -2)

#Save dataframe
write.xlsx(cDNA_transcript_count_filter_top, paste0('/media/kilian/OS/Melani_THP1_August_2023/cDNA_ANALYSIS/stringtie_output_cDNA_specific_L1/cDNA_transcript_top.xlsx'))



##### ---------------------------------IMPORT DATA TRANSCRIPTS---------------------------------- #####
#cDNA_transcript_count_filter_top$transcript_id <- gsub('\\/','', cDNA_transcript_count_filter_top$transcript_id)
goi <- cDNA_transcript_count_filter_top %>% head(20) %>% dplyr::pull('transcript_id')
goi_data <-  cDNA_transcript_count_filter_top %>% filter(transcript_id %in% goi)

goi_data2 <- dplyr::rename(goi_data,
                           'S1_KO-' = dplyr::contains('S1_norm'),
                           'WT-' = dplyr::contains('WT_norm'))

goi_data2 <- goi_data2 %>% dplyr::select('transcript_id','S1_KO-','WT-') 

goi_data2

lineWidth = 1
pointSize = 14

dir.create(paste('cDNA_ANALYSIS/plots/barplots/L1_telescope'))

for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(transcript_id == i) %>% pivot_longer(
    cols = -(1:1),
    values_to = c("DS_Counts"),
    names_to = c("condition", "replicate"),
    names_sep  = "-")
  
  barplot <- ggplot(goi_data3, aes(x = condition, y = DS_Counts, fill = condition)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    #geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_point(aes(x = condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Normalized Counts ±SE")) +
    ggtitle(paste0(i))+
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(size = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=20, face="bold.italic"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      # legend.key.height = unit(0.1, "cm"),
      # legend.key.width = unit(0.2, "cm"),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    #stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('cDNA_ANALYSIS/plots/barplots/L1_telescope/', i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}



##### ---------------------------------IMPORT GENE LOCI INFO---------------------------------- #####
rmsk_repeats <- as.data.frame(read.csv(paste0("/media/kilian/OS/Melani_THP1_August_2023/reference/L1_transcripts_telescope.gtf"), sep = "\t", header = T))
rmsk_repeats <- rmsk_repeats[,-c(2:3,6:8)]

colnames(rmsk_repeats) <-  c('Chromosome', 'Start', 'End', 'Name')
rmsk_repeats_edit<- cSplit(rmsk_repeats,'Name', ';')
rmsk_repeats_edit <- rmsk_repeats_edit[,-c(4:5,8:9)]


colnames(rmsk_repeats_edit) <-  c('Chromosome', 'Start', 'End', 'Locus','Category','L1base_ID')
rmsk_repeats_edit$Locus <- gsub('locus ','',rmsk_repeats_edit$Locus)
rmsk_repeats_edit$Category <- gsub('category ','',rmsk_repeats_edit$Category)
rmsk_repeats_edit$L1base_ID <- gsub('l1base_id ','',rmsk_repeats_edit$L1base_ID)


#Merging significant TEs with location information
colnames_cDNA_transcript_count_filter_top <- colnames(cDNA_transcript_count_filter_top)[-1]
colnames(cDNA_transcript_count_filter_top) <-  c('Locus', colnames_cDNA_transcript_count_filter_top)

merged_L1_final <- dplyr::left_join(cDNA_transcript_count_filter_top, rmsk_repeats_edit, 'Locus')

write.xlsx(merged_L1_final, paste('cDNA_ANALYSIS/results/L1_filtered_Log2_counts.xlsx'))

##### ---------------------------------VIZUALISATION OF GENE LOCI INFO---------------------------------- #####
#VOLCANO PLOT
##Theme
theme_scatter <- theme(aspect.ratio = 1, 
                       panel.background = element_blank(),
                       panel.border=element_rect(fill=NA),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background=element_blank(),
                       axis.title  = element_text(size = 20, colour = "black"),
                       axis.text.x=element_text(colour="black",size =20),
                       axis.text.y=element_text(colour="black",size =20),
                       axis.ticks=element_line(colour="black"),
                       plot.margin=unit(c(1,1,1,1),"line"),
                       plot.title = element_text(colour="black",size =20))

###SELECT LABELS
select_UP <- merged_L1_final %>% arrange(desc(LOG2_FC)) %>% head(10) %>% dplyr::pull('Locus')
select_DOWN <- merged_L1_final %>% arrange(LOG2_FC) %>% head(10) %>% dplyr::pull('Locus')


#Select data of labels
select_UP_data <- merged_L1_final %>% filter(Locus %in% c(select_UP))
select_DOWN_data <- merged_L1_final %>% filter(Locus %in% c(select_DOWN))

#Select dimensions
col_max <- max(merged_L1_final$LOG2_FC)
Volcano_prefix <- 'L1_elements_full_length'

scatter_volcano <- ggplot(merged_L1_final, aes(x= LOG2_FC, y= SUM, label = Locus, size = SUM)) +
  geom_point(data = merged_L1_final,aes(size = SUM,fill = LOG2_FC),  shape = 21) +
  geom_point(data = select_UP_data , aes(fill = LOG2_FC,size = SUM), shape = 21, color = 'black') +
  geom_point(data = select_DOWN_data , aes(fill = LOG2_FC,size = SUM), shape = 21, color = 'black') +
  scale_fill_gradient2(limits=c(-col_max, col_max), low="navyblue", mid="whitesmoke", high = "firebrick", na.value = 'navy') +
  theme_scatter +
  #xlim(-15,15)+
  #ylim(0,9)+
  geom_vline(xintercept = 0.1,colour="grey", linetype = "longdash") +
  geom_vline(xintercept = -0.1,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = 1.3, colour="grey", linetype = "longdash") +
  labs(title = str_wrap(paste0(Volcano_prefix),60), x = "Log2FC", y = "SUM_count") +
  geom_text_repel(data = select_UP_data,
                  color = 'black', size = 5, fontface = 'italic',
                  min.segment.length = 0,
                  nudge_y = 5,
                  #nudge_x = 14 - select_UP_data$log2FoldChange,
                  #nudge_x = 2,
                  max.overlaps = Inf,
                  box.padding = unit(0.5, 'mm'),
                  point.padding = unit(0.5, 'mm'),
                  segment.linetype = 5) +
  geom_text_repel(data = select_DOWN_data,
                  color = 'black',size = 5, fontface = 'italic',
                  min.segment.length =0,
                  nudge_y = 5,
                  #nudge_x = -14 + abs(select_DOWN_data$log2FoldChange),
                  #nudge_x = -2,
                  max.overlaps = Inf,
                  point.padding = unit(0.5, 'mm'),
                  segment.linetype = 5)
scatter_volcano 

Volcano_prefix
ggsave(paste0('cDNA_ANALYSIS/plots/volcanos/','L1_Stringtie_output', '.pdf'),
       plot = scatter_volcano,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)

#DONUT PLOT 
#Piechart showing families 
merged_L1_final$UP_or_DOWN[merged_L1_final$LOG2_FC>0] <- 'UP'
merged_L1_final$UP_or_DOWN[merged_L1_final$LOG2_FC<0] <- 'DOWN'

donut_prefix <- 'L1_overall_change'
donut_data <- merged_L1_final

hsize <- 4

donut_data <- donut_data %>% 
  group_by(UP_or_DOWN) %>% # Variable to be transformed
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc)) %>% 
  mutate(x = hsize)


# Basic Donut of TE elements
piechart <- ggplot(donut_data, aes(x=hsize, y= perc, fill=UP_or_DOWN, label = labels)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  geom_label(aes(label = labels),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  labs(title = paste(donut_prefix))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  scale_fill_manual(values = c('blue','firebrick'))+
  #scale_fill_brewer(palette = "Paired") +
  xlim(c(0.2, hsize + 0.5))


ggsave(paste0('cDNA_ANALYSIS/plots/donut/',donut_prefix, '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

#Family level 
donut_data_UP <- merged_L1_final %>% filter(UP_or_DOWN == 'UP')
donut_data_DOWN <-  merged_L1_final %>% filter(UP_or_DOWN == 'DOWN')

donut_prefix <- 'donut_data_DOWN'
donut_TEclass <- 'L1_category'
donut_data <- donut_data_DOWN

hsize <- 4

donut_data <- donut_data %>% 
  group_by(Category) %>% # Variable to be transformed
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc)) %>% 
  mutate(x = hsize)


# Basic Donut of TE elements
piechart <- ggplot(donut_data, aes(x=hsize, y= perc, fill=Category, label = labels)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  geom_label(aes(label = labels),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  labs(title = paste(donut_prefix))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  scale_fill_brewer(palette = "Paired") +
  xlim(c(0.2, hsize + 0.5))


ggsave(paste0('cDNA_ANALYSIS/plots/donut/',donut_prefix, '_',donut_TEclass, '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)



