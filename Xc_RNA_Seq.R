library(edgeR)
library(magrittr)
library(pheatmap)
library(ggplot2)



#Set directory to where output of FeatureCounts is located
setwd("/Users/julius/Desktop/Hind_Cucurbita_RNA-Seq_Cleaned_Data/FeatureCountsv2/")


#---------------------------------------------

###import FeatureCounts results, change headers and row names, and get rid of extra columns
read.counts <- read.table ( "featureCounts_results.txt", header = TRUE )
row.names(read.counts) = read.counts$Geneid

read.counts = read.counts [ , -c (1:6) ]
orig_names = names(read.counts)

names(read.counts) = c("Control 1", "Control 2", "Control 3", "Hrp 1", "Hrp 2", "Hrp 3")
sample_names= c("Control 1", "Control 2", "Control 3", "Hrp 1", "Hrp 2", "Hrp 3")

#---------------------------------------------

#Put FeatureCounts data into an EdgeR data structure

sample_info.edger <- factor (c( rep( "PS" , 3) , rep( "XVM" , 3) ) )
sample_info.edger<- relevel(sample_info.edger , ref = "PS" )


#---------------------------------------------


#Carry out DGEList
y =  DGEList(counts = read.counts , group = sample_info.edger )
y = calcNormFactors(y)

#Get rid of genes in y that have no expression in any of the samples
keep = rowSums(cpm(y$counts) > 0) >= 1
y =y[keep,]


#Lets us know which samples are C or E

group =c(0,0,0,1,1,1)
group_matrix = model.matrix(~group)
y = estimateDisp(y, design = group_matrix)

plotMDS(y, method="bcv", col=group_matrix)
plotBCV(y)

#Run DE Tests

fit = (glmFit (y, group_matrix, dispersion = y$trended.dispersion))
de = glmLRT(fit,coef = 2) 


out = topTags(de, n = nrow(y), adjust.method = "BH",sort.by="PValue", p.value=0.05)

#Get top 30 DEGs
temp = order(-abs(out$table$logFC))[1:30]

deGenes=decideTestsDGE(de, adjust.method="BH", p.value=0.05)



#Write out Tables into csv and tsv files

#write.csv(out, file="Xcuc_toptags_latest.csv")
#write.table(out, file='Xcuc_toptags.tsv', quote=FALSE, sep='\t', col.names = NA)
#write.csv(de$table, file="Xcuc_All_Data.csv")


#---------------------------------------------

### VISUALIZATIONS ###

#---------------------------------------------

theme_update(plot.title = element_text(hjust = 0.5))

#Generate Bar graph

# Get a list of DGE and if they are Up or Down Regulated, from decidetestDGE

Bar_graph_data=summary(deGenes)
Bar_graph_data = as.data.frame.matrix(Bar_graph_data)
Bar_graph_data$type=c("Downregulated","Not Significant","Upregulated")
Bar_graph_data = Bar_graph_data[-c(2),]
Bar_graph_data = Bar_graph_data[c(2,1),]

cols <- c("Downregulated"="#00BFC4","Upregulated"="#F8766D")

ggplot(data=Bar_graph_data, aes(x=type, y=group, fill=type, width = 0.7)) + 
  geom_bar(stat="identity",show.legend = FALSE) + 
  geom_text(size=5, aes(label=group),position = position_stack(vjust = 0.5)) + 
  labs( y="Number of DEGs",x ="") +
  theme(axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=12)) +
        scale_fill_manual(values = cols)


#---------------------------------------------

#Generate Smear Plot


de_w_Avelog= de
ave_Log_CPM_col = aveLogCPM(y)

#Add a column for Ave Log CPM
de_w_Avelog$table$ave_Log_CPM=ave_Log_CPM_col

DEGs_Up <- de_w_Avelog$table$PValue < 0.05 & de_w_Avelog$table$logFC > 1.5
DEGs_volcano_Up = de_w_Avelog$table
DEGs_volcano_Up= DEGs_volcano_Up[DEGs_Up,]


DEGs_Down <- de_w_Avelog$table$PValue < 0.05 & de_w_Avelog$table$logFC < -1.5
DEGs_volcano_Down = de_w_Avelog$table
DEGs_volcano_Down= DEGs_volcano_Down[DEGs_Down,]


smear = ggplot(de_w_Avelog$table,aes(x=ave_Log_CPM,y=logFC)) + 
  geom_point(size=1,shape=1) + 
  geom_point(data=DEGs_volcano_Down,color="#00BFC4",size=1, shape=1) + 
  geom_point(data=DEGs_volcano_Up,color="#F8766D",size=1, shape=1) + 
  labs(title="Smear Plot", y="LogFC", x="Average LogCPM")

smear + geom_hline(yintercept=1.5, linetype="dashed", color = "dark grey") + 
  geom_hline(yintercept=-1.5, linetype="dashed", color = "dark grey")

#---------------------------------------------

#Generate Volcano Plot

volcano = ggplot(de$table,aes(x=logFC,y=-log10(PValue))) + 
  geom_point(size=1,shape=1) + 
  geom_point(data=DEGs_volcano_Down,color="#00BFC4",size=1, shape=1) + 
  geom_point(data=DEGs_volcano_Up,color="#F8766D",size=1, shape=1) + 
  labs( y="-log(P-value)", x="LogFC") 

volcano + geom_vline(xintercept=1.5, linetype="dashed", color = "dark grey") + 
  geom_vline(xintercept=-1.5, linetype="dashed", color = "dark grey") +
  geom_hline(yintercept=1.30102999566, linetype="dashed", color = "dark grey") +
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=12))


#---------------------------------------------

#Generate Heat Maps 


# get all genes that don't have NA
tmpy = y$counts
keep = which(apply(y$counts, 1, sd) > 0)
tmpy = tmpy[keep,]


#Colors for annotation
ann_colors = list(
  Treatment = c(Control="#1B9E77", Hrp="#FFB6C1"))  #GREEN AND PINK

#---------------------------------------------

#This is for individual rep annotation
my_sample_col <- data.frame(sample = rep(c("Control", "Hrp"), c(3,3)))
row.names(my_sample_col) <- colnames(y)
names(my_sample_col) <-c("Treatment")

pheatmap(tmpy, scale= "row", 
         clustering_distance_rows= "correlation",
         show_rownames = F,
         annotation_col = my_sample_col,
         annotation_colors = ann_colors)

#---------------------------------------------

#Select only genes with logFC > 1.5
selY <- y[rownames(out$table)[out$table$FDR<.1 & abs(out$table$logFC)>5],]
pheatmap(selY$counts, scale= "row", 
         main="Heatmap of highly changed genes", 
         clustering_distance_rows= "correlation", 
         annotation_col = my_sample_col,show_rownames = F)

