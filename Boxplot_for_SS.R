library(ggplot2)
library(ggpubr)
library(gridExtra) 
library(dplyr)
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("gridExtra")
install.packages("dplyr")

ss_codon <- read.table("/Your/work/path/yeast.extract_SS_relativeASA.txt",
                       header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

colnames(ss_codon) <-c("RelativePosition","AminoAcid","SecondaryStructure","RelativeASA", "ReadsCount")
correlation<-ss_codon[,c(3, 4,5)]

correlationlog<-mutate(correlation,ReadsCount=ReadsCount+0.01)
correlationlog2<- mutate(correlationlog,ReadsCount = log10(ReadsCount) )
colnames(correlationlog2) <-c("Secondary_Structure","RelativeASA", "Log_Normalized_ReadsCount")

correlationlog_ss<-correlationlog2[,c(1, 3)]
colnames(correlationlog_ss) <-c("Secondary_Structure", "Log_Normalized_ReadsCount")
correlationlog_ss$Secondary_Structure<-factor(correlationlog_ss$Secondary_Structure,
                                              levels = c("-","B","E","G","H","I","S","T"),
                                              labels = c("Coil","B-bridge","Beta-sheet","3-helix","Alpha-helix","5-helix","Bend","Turn"))

correlationlog_ss_use= dplyr::filter(correlationlog_ss, (Secondary_Structure == "Coil" | Secondary_Structure == "Beta-sheet" | Secondary_Structure == "Alpha-helix"))


my_comparisons = list(c("Coil", "Beta-sheet"), c("Coil", "Alpha-helix"))
                

p4<-ggplot(correlationlog_ss_use, aes(x = Secondary_Structure, y = Log_Normalized_ReadsCount,col=Secondary_Structure, fill = Secondary_Structure)) +
  geom_boxplot(alpha=0.8)+stat_compare_means(comparisons = my_comparisons)+labs(title = "Rat")+ 
  theme(plot.title=element_text(hjust=0.5))+ylab("Logarithm Scaled Footprint")+ xlab("Protein Secondary Structure")


p5<-p4+theme(axis.text = element_text(size = 15,face="bold"))+theme(axis.title = element_text(size = 10,face="bold"))+theme(legend.text = element_text(size = 15,face="bold"))+ theme(legend.title = element_text(size = 15,face="bold"))+theme_bw()+ theme(panel.grid=element_blank())

p6<- p5+theme(legend.position = "none")

p7 <- p6 + theme(
  panel.border = element_rect(linewidth = 1, color = "black") )

png("/Your/work/path/Plot/yeast.SS.png", width = 1500, height = 1000, res = 200)

print(p7)
dev.off()



# compare the association between secondary structures
coil_reads <-correlationlog_ss [correlationlog_ss$Secondary_Structure=="Coil",]
coil <- coil_reads[,c(2)]

beta_sheet_reads <-correlationlog_ss [correlationlog_ss$Secondary_Structure=="Beta-sheet",]
beta_sheet <-beta_sheet_reads[,c(2)]

alpha_helix_reads <-correlationlog_ss [correlationlog_ss$Secondary_Structure=="Alpha-helix",]
alpha_helix <-alpha_helix_reads[,c(2)]

turn_reads <- alpha_helix_reads <-correlationlog_ss [correlationlog_ss$Secondary_Structure=="Turn",]
turn <-turn_reads[,c(2)]

wilcox.test(coil,alpha_helix,paired = FALSE,alternative = 'two.sided')
wilcox.test(coil,alpha_helix,paired = FALSE,alternative = 'great')
wilcox.test(coil,alpha_helix,paired = FALSE,alternative = 'less')

wilcox.test(coil,beta_sheet,paired = FALSE,alternative = 'two.sided')
wilcox.test(coil,beta_sheet,paired = FALSE,alternative = 'great')
wilcox.test(coil,beta_sheet,paired = FALSE,alternative = 'less')

wilcox.test(coil,turn,paired = FALSE,alternative = 'two.sided')
wilcox.test(coil,turn,paired = FALSE,alternative = 'great')
wilcox.test(coil,turn,paired = FALSE,alternative = 'less')








