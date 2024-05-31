setwd("D:/submission/others/energy/NC31/Figure5/NC31")
library(microeco)
#?microtable # 展示所有的功能和细节描述
library(magrittr)
library(ggplot2)
library(randomcoloR)
theme_set(theme_bw())

#import data
otu_table_16S_1<-read.csv("DATA-ORIGIN.csv",row.names = 1)
sample_info_16S_1<-read.csv("SampleID-ORIGIN.csv",row.names = 1)
taxonomy_table_16S_1<-read.csv("tax-ORIGIN.csv",row.names = 1)
dataset <- microtable$new(sample_table = sample_info_16S_1, 
                          otu_table = otu_table_16S_1,
                          tax_table = taxonomy_table_16S_1)

#group
dataset$sample_table$Group <-  factor(dataset$sample_table$Group, levels = c("HF","NY15","HeNa1610","FXJWS23M45","NC31"))

#touch data
dataset$tidy_dataset()
dataset$cal_abund()
class(dataset$taxa_abund)
#创建丰度表格
#dir.create("taxa_abund")
#dataset$save_abund(dirpath = "taxa_abund")


#color
color_compaired <- c("#09f9f5", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", "#e86502", "#931635", "#373bbf", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f","#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#a1ce4c", "#ef3bb6", "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")


#alpha_dicersity
dataset$tidy_dataset()
print(dataset)
dataset$rarefy_samples(sample.size = 10000)
dataset$sample_sums() %>% range
dataset$cal_alphadiv(PD = FALSE)
dir.create("alpha_diversity")
dataset$save_alphadiv(dirpath = "alpha_diversity")


#phylum
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)

#class
t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 15)
t1$plot_box(group = "Group")


#genus
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40,groupmean = "Group")
t1$plot_bar(legend_text_italic = FALSE, color_values = c(color_compaired),xtext_size = 6, others_color = "grey70", use_alluvium = T)



#venn
dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()


