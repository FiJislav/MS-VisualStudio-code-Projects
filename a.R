library(dendextend)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library("grid") 
library("heatmaply") 
library("ggdendro") 
library(gridExtra) 



input_list<-read_xlsx("D:/1_DATA/Rcoding/Heatmap_try/Book3_fin.xlsx", col_names = T) 


df <- data.frame(input_list)

df_fasta_short <- df %>%
  mutate(df, fastaShort = sub(" OS=.*", "", df$Fasta_headers)) #
iBAQ_exp <- format(df_fasta_short$iBAQ, scientific = TRUE, digits = 3)
df_fasta_short <- cbind(iBAQ_exp, df_fasta_short)



df_fasta_short_sorted<- arrange(df_fasta_short,Mol_weight)

#df_fasta_short_sorted_normalized <- df_fasta_short_sorted %>% mutate_at(7:46, funs((.-min(.))/max(.-min(.)))) #to normalize on columns


df_fasta_short_sorted_united <- unite(df_fasta_short_sorted, col='fastaShort_MW_iBAQ_rank', c('fastaShort', 'Mol_weight', 'iBAQ_exp', 'rank'), sep='   ', remove = T)


normrows<-t(apply(df_fasta_short_sorted_united[,7:46], 1, function(x)(x-min(x))/(max(x)-min(x))))

df_fasta_short_sorted_united[, 7:46] <- normrows

df_fasta_short_sorted_united<-df_fasta_short_sorted_united %>% drop_na("FR1","FR2","FR3","FR4","FR5","FR6","FR7","FR8","FR9","FR10","FR11","FR12","FR13","FR14","FR15","FR16","FR17","FR18","FR19","FR20","FR21","FR22","FR23","FR24","FR25","FR26","FR27","FR28","FR29","FR30","FR31","FR32","FR33","FR34","FR35","FR36","FR37","FR38","FR39","FR40")

# Run clustering
dendro_matrix <- as.matrix(df_fasta_short_sorted_united[, -c(1:6)])
rownames(dendro_matrix) <- df_fasta_short_sorted_united$fastaShort_MW_iBAQ_rank
dendro_dendro <- as.dendrogram(hclust(d = dist(x = dendro_matrix)))

# Create dendro
dendro_plot <- ggdendrogram(data = dendro_dendro, rotate = T)+
  theme(axis.text.y = element_text(size = 4.5))+
  theme(axis.text.x = element_blank())
# Preview the plot
print(dendro_plot)


ggsave("dendro_Glyco_eploooremsVisualCode.svg", device="svg", dpi=1000, width=20, height=41, bg='transparent',units = "cm")


df_for_heatmap <- pivot_longer(data = df_fasta_short_sorted_united,
cols = -c(Fasta_headers, fastaShort_MW_iBAQ_rank, Sequence_coverage, Intensity, iBAQ_peptides, iBAQ),
                           names_to = "Fractions",
                           values_to = "iBAQ_values")

df_for_heatmap<-df_for_heatmap %>% drop_na("iBAQ_values")


df_for_heatmap$Fractions<-factor(df_for_heatmap$Fractions, levels=c("FR1","FR2","FR3","FR4","FR5","FR6","FR7","FR8","FR9","FR10","FR11","FR12","FR13","FR14","FR15","FR16","FR17","FR18","FR19","FR20","FR21","FR22","FR23","FR24","FR25","FR26","FR27","FR28","FR29","FR30","FR31","FR32","FR33","FR34","FR35","FR36","FR37","FR38","FR39","FR40"))

heatmap_plot <- ggplot(data = df_for_heatmap, aes(x = Fractions, y = fastaShort_MW_iBAQ_rank))+
  geom_tile(aes(fill = iBAQ_values))+
  scale_fill_gradient2(low = "black",high = "red")+
  theme(axis.text.y = element_text(size = 2.5))+
  theme(axis.text.x = element_text(size = 4))+
  ylab("Identified proteins")+
  xlab("SEC fractions")+
  guides(fill = guide_colourbar(title = "normalized iBAQ"))+
  theme(panel.border = element_rect(fill=NA),axis.text.x = element_text(angle=90))
print(heatmap_plot)

dendro_order <- order.dendrogram(dendro_dendro)

df_for_heatmap$fastaShort_MW_iBAQ_rank <- factor(x = df_for_heatmap$fastaShort_MW_iBAQ_rank,
                               levels = df_fasta_short_sorted_united$fastaShort_MW_iBAQ_rank[dendro_order], 
                               ordered = TRUE)


heatmap_plot <- ggplot(data = df_for_heatmap, aes(x = Fractions, y = fastaShort_MW_iBAQ_rank))+
  geom_tile(aes(fill = iBAQ_values))+
  #scale_fill_viridis(option = "viridis")+   #vitidis_scale
  scale_fill_gradient2(low = "black",high = "red")+
  theme(axis.text.y = element_text(size= 4.5), legend.position = "top")+
  theme(axis.text.x = element_text(size = 6))+
  ylab("Identified proteins")+
  xlab("SEC fractions")+
  guides(fill = guide_colourbar(title = "normalized iBAQ"))+
  theme(panel.border = element_rect(fill=NA),axis.text.x = element_text(angle=90))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
print(heatmap_plot)
ggsave("D:/1_DATA/Rcoding/Heatmap_try/heatmap_Glyco_unknosdwn.svg", device="svg", dpi=1000, width=20, height=41, bg='transparent',units = "cm") 