
if (T) {
  dir.create("scripts")
  #dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("00_00_origin_datas/GEO",recursive = T)
  dir.create("00_00_origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RisktypeProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_mutibarplot=function(df,xlab='group',leg.title='',cols=pal_d3()(10)[5:6]){
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
##########
genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)


###TCGA#############
#######
tcga.tumor.fpkm=read.delim('00_origin_datas/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt',check.names = F)
tcga.tumor.fpkm[1:5,1:5]
tcga.tumor.fpkm=tcga.tumor.fpkm[which(tcga.tumor.fpkm[, 1]!=''),] 
tcga.tumor.fpkm=tcga.tumor.fpkm[!duplicated(tcga.tumor.fpkm$Hugo_Symbol),]
rownames(tcga.tumor.fpkm)=tcga.tumor.fpkm$Hugo_Symbol
tcga.tumor.fpkm=tcga.tumor.fpkm[,-c(1,2)]
dim(tcga.tumor.fpkm)
tcga.normal.fpkm=read.delim('00_origin_datas/brca_tcga_pan_can_atlas_2018/normals/data_mrna_seq_v2_rsem_normal_samples.txt',check.names = F)
tcga.normal.fpkm[1:5,1:5]
tcga.normal.fpkm=tcga.normal.fpkm[which(tcga.normal.fpkm[, 1]!=''),] 
tcga.normal.fpkm=tcga.normal.fpkm[!duplicated(tcga.normal.fpkm$Hugo_Symbol),]
rownames(tcga.normal.fpkm)=tcga.normal.fpkm$Hugo_Symbol
tcga.normal.fpkm=tcga.normal.fpkm[,-c(1,2)]

com.genes=intersect(rownames(tcga.tumor.fpkm),rownames(tcga.normal.fpkm))
length(com.genes)
#20486

tcga.fpkm=cbind(tcga.tumor.fpkm[com.genes,],tcga.normal.fpkm[com.genes,])
dim(tcga.fpkm)
range(tcga.fpkm)
tcga.fpkm=log2(tcga.fpkm+1)


tcga.type=data.frame(Samples=c(colnames(tcga.tumor.fpkm),colnames(tcga.normal.fpkm)),
                     Type=rep(c('Tumor','Normal'),c(length(colnames(tcga.tumor.fpkm)),length(colnames(tcga.normal.fpkm)))))
rownames(tcga.type)=tcga.type$Samples
table(tcga.type$Type)

#####
tcga.cli=read.delim('00_origin_datas/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018_clinical_data.tsv',check.names = F)
head(tcga.cli)
table(tcga.cli$`Cancer Type Detailed`)
tcga.cli=data.frame(Samples=tcga.cli$`Sample ID`,Age=tcga.cli$`Diagnosis Age`,
                    Stage=tcga.cli$`Neoplasm Disease Stage American Joint Committee on Cancer Code`,
                    T.stage=tcga.cli$`American Joint Committee on Cancer Tumor Stage Code`,
                    N.stage=tcga.cli$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`,
                    M.stage=tcga.cli$`American Joint Committee on Cancer Metastasis Stage Code`,
                    OS=str_split_fixed(tcga.cli$`Overall Survival Status`,':',2)[,1],
                    OS.time=as.numeric(tcga.cli$`Overall Survival (Months)`/12*365),
                    DSS=str_split_fixed(tcga.cli$`Disease-specific Survival status`,':',2)[,1],
                    DSS.time=as.numeric(tcga.cli$`Months of disease-specific survival`/12*365),
                    PFS=str_split_fixed(tcga.cli$`Progression Free Status`,':',2)[,1],
                    PFS.time=as.numeric(tcga.cli$`Progress Free Survival (Months)`/12*365))
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli[tcga.cli$OS.time>30 & tcga.cli$OS.time<3650,]

table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='STAGE X']=NA
tcga.cli$Stage=gsub('STAGE ','',tcga.cli$Stage)
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)

table(tcga.cli$T.stage)
tcga.cli$T.stage[tcga.cli$T.stage=='TX']=NA
tcga.cli$T.stage=gsub('[ABCD]','',tcga.cli$T.stage)

table(tcga.cli$N.stage)
tcga.cli$N.stage[tcga.cli$N.stage=='NX']=NA
tcga.cli$N.stage[tcga.cli$N.stage=='N0 (MOL+)']='N0'
tcga.cli$N.stage[tcga.cli$N.stage=='N1MI']='N1'
tcga.cli$N.stage=gsub('[ABCI()+-]','',tcga.cli$N.stage)
tcga.cli$N.stage=gsub(' ','',tcga.cli$N.stage)

table(tcga.cli$M.stage)
tcga.cli$M.stage[tcga.cli$M.stage=='MX']=NA
tcga.cli$M.stage=gsub('[ABCI()+-]','',tcga.cli$M.stage)
tcga.cli$M.stage=gsub(' ','',tcga.cli$M.stage)

######
tcga.sample=intersect(tcga.cli$Samples,colnames(tcga.fpkm))
length(tcga.sample)
#994
range(tcga.fpkm)
tcga.exp=tcga.fpkm[intersect(rownames(tcga.fpkm),mrna_genecode$SYMBOL),tcga.sample]
dim(tcga.exp)
tcga.cli=tcga.cli[tcga.sample,]
dim(tcga.cli)

# table(tcga.cli$OS.time>30,tcga.cli$OS.time<3650)

#########GSE20685#####
load('00_origin_datas/GEO/GSE20685.RData')
GSE20685.cli=pData(GSE20685)
head(GSE20685.cli)
GSE20685.cli=data.frame(Samples=GSE20685.cli$geo_accession,
                          Age=GSE20685.cli$`age at diagnosis:ch1`,
                          T.stage=GSE20685.cli$`t_stage:ch1`,
                          N.stage=GSE20685.cli$`n_stage:ch1`,
                          M.stage=GSE20685.cli$`m_stage:ch1`,
                          OS.time=GSE20685.cli$`follow_up_duration (years):ch1`,
                          OS=GSE20685.cli$`event_death:ch1`)
rownames(GSE20685.cli)=GSE20685.cli$Samples
head(GSE20685.cli)
GSE20685.cli=crbind2DataFrame(GSE20685.cli)
GSE20685.cli$OS.time=GSE20685.cli$OS.time*365
GSE20685.cli$T.stage=paste0('T',GSE20685.cli$T.stage)
GSE20685.cli$N.stage=paste0('N',GSE20685.cli$N.stage)
GSE20685.cli$M.stage=paste0('M',GSE20685.cli$M.stage)

GSE20685.df=exprs(GSE20685)
GSE20685.df[1:5,1:5]
range(GSE20685.df)
GSE20685.exp=exp_probe2symbol_v2(datExpr = GSE20685.df,GPL = 'GPL570')
dim(GSE20685.exp)


####01.##########
dir.create('01_GBP_gene')
GBP.gene <- read.table("01_GBP_gene/GBP_Gene.txt",sep = "\t")
GBP.gene <- GBP.gene$V1
length(GBP.gene)
#7
GBP.gene <- intersect(GBP.gene,rownames(tcga.exp))
######1.1######################
GBP.gene.df <- data.frame(t(tcga.fpkm[GBP.gene,tcga.type$Samples]),group =tcga.type$Type)
GBP.gene.df <- melt(GBP.gene.df)
dodge_width <-1.5
ggplot(GBP.gene.df,aes(x=variable,y=value,fill=group))+
  geom_violin(position = position_dodge(dodge_width), trim = FALSE,show.legend = T) +
  geom_boxplot(width = 0.3, position = position_dodge(dodge_width), outlier.shape = NA,show.legend = F) +
  scale_fill_manual(values = c("skyblue","#FF8000"))+
  ggpubr::stat_compare_means(aes(group=group), label = 'p.signif', method = 't.test') +
  facet_wrap(~variable, scales = 'free', nrow = 3, ncol = 3)
ggsave("01_GBP_gene/expression.GBP.gene.tumor_normal.pdf",height = 10,width = 11)

#####1.2############
tcga.maf=getTCGAMAFByCode('BRCA')
tcga.type.use=tcga.type
table(tcga.type.use$Type)
colnames(tcga.type.use)[1]='Tumor_Sample_Barcode'
tcga.type.use$Tumor_Sample_Barcode=substr(tcga.type.use$Tumor_Sample_Barcode,1,12)
tcga.type.use.Normal=tcga.type.use[which(tcga.type.use$Type=='Normal'),]
tcga.type.use.Tumor=tcga.type.use[which(tcga.type.use$Type=='Tumor'),]
write.table(tcga.type.use.Normal,file='01_GBP_gene/tcga.type.Normal.txt')
write.table(tcga.type.use.Tumor ,file='01_GBP_gene/tcga.type.Tumor.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.type.use.Normal$Tumor_Sample_Barcode))
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = '01_GBP_gene/tcga.type.Normal.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.type.use.Tumor $Tumor_Sample_Barcode))
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = '01_GBP_gene/tcga.type.Tumor.txt')
tcga.maf2@clinical.data


#type.color=c('rosybrown1','lightgoldenrod1')
pdf('01_GBP_gene/normal_mut.pdf',height = 5,width =6)
oncoplot(maf=tcga.maf1,
         #clinicalFeatures = 'type',
         genes = GBP.gene,
         sortByAnnotation = T,
         #annotationColor = list(type=c(Normal=as.character(type.color[1])))
)
dev.off()
pdf('01_GBP_gene/tumor_mut.pdf',height = 5,width =6)
oncoplot(maf=tcga.maf2,
         #clinicalFeatures = 'type',
         genes = GBP.gene,
         sortByAnnotation = T,
         #annotationColor = list(type=c(Tumor =as.character(type.color[2])))
)
dev.off()
##################1.3############
GBP.score=ssGSEAScore_by_genes(tcga.fpkm,genes = GBP.gene)
GBP.score <- as.data.frame(t(GBP.score))
p1a <- my_violin(GBP.score ,group=tcga.type$Type,ylab='GBP scores')
ggsave('01_GBP_gene/p1a.pdf',height = 6,width = 6)
###################1.4############
#######
GBP.score.cli=data.frame(GBP.score=GBP.score[tcga.cli$Samples,],tcga.cli[,c('OS','OS.time')])
GBP.score.cli$group=ifelse(GBP.score.cli$GBP.score>median(GBP.score.cli$GBP.score),'High','Low')
GBP.score.cli <- as.data.frame(GBP.score.cli)
GBP.score.cli$OS <- as.numeric(GBP.score.cli$OS)
str(GBP.score.cli)
GBP.score.roc=ggplotTimeROC(GBP.score.cli$OS.time,
                            GBP.score.cli$OS,
                            GBP.score.cli$GBP.score,mks = c(1,3,5))
GBP.score.roc
GBP.score.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                       data = GBP.score.cli),
                           data=GBP.score.cli,
                           conf.int = T,pval = T,risk.table = T,
                           fun = "pct",size = 1,surv.median.line = 'hv',
                           title='TCGA-BRCA',legend.title='GBP score',
                           legend.labs = c('High','Low'),
                           linetype = c("solid", "dashed","strata")[1],
                           #palette = risktype.col,
                           ylab='Overall Survival(OS)',
                           legend=c(0.85,0.8),#标签位置
                           ggtheme = theme_bw(base_size = 12))
GBP.score.km.OS=mg_merge_plot(GBP.score.km.OS$plot,GBP.score.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

GBP.score.km.OS
ggsave("01_GBP_gene/GBP.score.km.OS.pdf",GBP.score.km.OS,height = 6,width = 6)

####02WGCNA###########
dir.create("02_WGCNA")

library(WGCNA)
allowWGCNAThreads(nThreads = 36)#
enableWGCNAThreads(nThreads = 36)# 

my_mad <- function(x){mad(x,na.rm = TRUE)} #
wgcna_exp=t(tcga.exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga.exp)
# 18448   500
#
tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]

#tpm_T2=tcga.exp.HBV[which(apply(tcga.exp.HBV,1,sd)>0.5),]
#tpm_T2=(2^tpm_T2-1)
range(tpm_T2)
pdf('02_WGCNA//wgcna.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.2,
                                 minModuleSize=60)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))
#16
write.csv(tpm_T2.module$Modules,file = "02_WGCNA//WGCNA_Modules.csv",col.names = T,row.names = T)
pdf('02_WGCNA/2.pdf',height = 6,width = 10)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
writeMatrix(tpm_T2.module$Modules,outpath = '02_WGCNA/tcga.wgcna.module.genes.txt')

pdf('02_WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02_WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module genes",xlab = "", sub = "")
dev.off()


############### 
tcga_cli_use <-data.frame(tcga.cli[tcga.cli$Samples,c(1,8)],GBP.score=GBP.score[tcga.cli$Samples,])
head(tcga_cli_use)
tcga_cli_use=as.data.frame(tcga_cli_use[,-c(1,2)])
colnames(tcga_cli_use)="GBP.score"
tcga_cli_use.part=tcga_cli_use
str(tcga_cli_use.part)
#tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))


spms=tcga_cli_use.part
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms[]
                  ,use = 'pairwise.complete.obs')

modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02_WGCNA/5.pdf',width = 6,height =12)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))
head(geneTraitSignificance)
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

####lightcyan####
module="lightcyan"
column= match(module, modNames)


moduleGenes <- c(tpm_T2.module$Modules[,'mergedColors']==module)
greenyellow.gene=names(which(moduleGenes))
length(greenyellow.gene)
# 81
writeMatrix(greenyellow.gene,'02_WGCNAlightcyan.gene.txt',header = F)
table(tpm_T2.module$Modules[,'mergedColors'])

pdf('02_WGCNA/6_lightcyan.pdf',height = 6,width = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, ]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for GBP.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = '#7BD3EA',lwd=2)
dev.off()

####magenta####
module="magenta"
column= match(module, modNames)


moduleGenes <- c(tpm_T2.module$Modules[,'mergedColors']==module)
magenta.gene=names(which(moduleGenes))
length(magenta.gene)
#1496
writeMatrix(magenta.gene,'02_WGCNAlightcyan.gene.txt',header = F)
table(tpm_T2.module$Modules[,'mergedColors'])

pdf('02_WGCNA/6_magenta.pdf',height = 6,width = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, ]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for GBP.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = 'magenta',lwd=2)
dev.off()
####
wgcna.gene <- c(greenyellow.gene,magenta.gene)
length(wgcna.gene)
#1577
  ####03_###########
dir.create('03_diff_genes')
tcga.limma=mg_limma_DEG(exp =tcga.fpkm[,tcga.type$Samples],
                        group=tcga.type$Type,ulab='Tumor',dlab = 'Normal')
tcga.limma$Summary

#######
tcga.Normal_Tumor.degs=tcga.limma$DEG
tcga.Normal_Tumor.degs=tcga.Normal_Tumor.degs[abs(tcga.Normal_Tumor.degs$logFC)>log2(2) & tcga.Normal_Tumor.degs$adj.P.Val<0.05,]
dim(tcga.Normal_Tumor.degs)
write.csv(tcga.Normal_Tumor.degs,'03_diff_genes/tcga.Normal_Tumor.degs.csv')

######3.1####
p_cutoff <-0.05 
fc_cutoff <- log2(2)
degs_dat=tcga.limma$DEG
degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                            ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))

#library(ggbreak)
library(ggplot2)
# library(ggprism)

col=c("hotpink","sienna","grey")
ylab='-log10 (adj.PVal)'
xlab='log2 (FoldChange)'
leg.pos='right'
plot2a<- ggplot(degs_dat, aes(x=logFC, y=-log10(adj.P.Val), color=type)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=col) +
  theme_bw() +
  theme(legend.position = leg.pos) +
  ylab(ylab) +
  xlab(xlab) +
  geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty=3, col="black", lwd=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), lty=3, col="black", lwd=0.5) 
#coord_cartesian(ylim=c(0, 25)) # 
#scale_y_break(c(25,100),#
# space = 0.3,#
# scales = 0.5)+#
#theme_prism(palette = "black_and_white",
# base_fontface = "plain", 
# base_family = "serif", 
# base_size = 16,
# base_line_size = 0.8,
# axis_text_angle = 0)
ggsave('03_diff_genes/fig3a.pdf',plot2a,height = 6,width = 6)
######3.2#######
venny.gene=intersect(wgcna.gene,rownames(tcga.Normal_Tumor.degs))
vennydata=list(WGCNA_module=wgcna.gene,TCGA_DEG=rownames(tcga.Normal_Tumor.degs))
library("VennDiagram")
pdf('03_diff_genes/venn.pdf',height = 7,width = 6)
# 
venn.plot <- venn.diagram(
  x = vennydata,
  filename = NULL,
  output = TRUE,
  fill = c("dodgerblue", "#80C4E9"),  # 
  alpha = 0.90,  # 
  cat.cex = 1.2,  # 
  cat.col = "gray10",  # 
  cat.fontface = "bold"  # 
)

# 
grid.draw(venn.plot)
dev.off()
######3.2 ########
length(venny.gene)
venny.gene.enrichment=mg_clusterProfiler(genes = venny.gene)
head(venny.gene.enrichment$GO_BP)
p3=list()
p3[[1]]=barplot(venny.gene.enrichment$GO_BP)
p3[[2]]=barplot(venny.gene.enrichment$KEGG)
p3cd=mg_merge_plot(p3,labels = c('C','D'),ncol = 2,widths = c(1,1))


write.xlsx(
  list(GO_BP=venny.gene.enrichment$GO_BP,
       KEGG=venny.gene.enrichment$KEGG),'03_diff_genes/venny.gene.enrichment.xlsx',
  overwrite = T)

savePDF('03_diff_genes/Fig3CD.pdf',p3cd,height = 8,width = 15)

####04.########
dir.create('04_model')
##########
tcga.cli$OS <- as.numeric(tcga.cli$OS)
sig.gene.cox=cox_batch(dat = tcga.exp[intersect(venny.gene,rownames(tcga.exp)),tcga.cli$Samples],
                       time = tcga.cli$OS.time,event = tcga.cli$OS)
sig.gene.cox



table(sig.gene.cox$p.value<0.05)
# FALSE  TRUE 
#   181   212 
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.05,])
length(pre.genes)#212
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
write.csv(sig.gene.cox,'04_model/sig.cox.csv')
# #####LASSO####
library(glmnet)
set.seed(2024)
fit1=glmnet(as.matrix(tcga_model_data[,-c(1,2)])
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1)

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,-c(1,2)])
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))
#8
pdf('04_model//LASSO.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(sig.coef), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
##########
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"-0.28*PSME2+-0.157*DACT2+-0.066*PIGR+-0.237*STX11"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)

###########
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(lan), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
fig_ggforest <- survminer::ggforest(cox,data=tcga_model_data)
ggsave("04_model/ggforest.pdf",height =4,width =7 )

####
risktype.col=c('#FFB200',"#704214")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
######KM####
tcga.data.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
tcga.cutoff <- as.numeric(summary(tcga.data.point)[1])
tcga.cutoff
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-BRCA',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),#
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS

tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.risktype.cli,cols =risktype.col,xlab = 'sample',
                         a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(tcga.risktype.cli$Riskscore),labs = '')
##
tcga_expr <- my_mutiboxplot( t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),notch = T,
                             group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                             ylab = 'Gene Expression',fill = 'Risktype',angle = 0,hjust = .5,title = 'TCGA')
ggsave('04_model/TCGA_expr.pdf',tcga_expr ,height = 5,width = 7)

###########
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
########
model.gene.df=data.frame(tcga.cli[,c('OS','OS.time')],t(tcga.exp[names(lan),tcga.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='TCGA',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.8,0.8),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=mg_merge_plot(module.gene.km[[i]]$plot,module.gene.km[[i]]$table,nrow=2,heights = c(3,1),align = 'v')
  
}
tcga.module.km=mg_merge_plot(module.gene.km[1:4],ncol=4,nrow=1)
roc_km_plot=mg_merge_plot(tcga.km.OS,tcga.roc,tcga_expr,ncol=3,nrow=1,labels = LETTERS[4:6],widths = c(1,1,1.5))
fig4 <- mg_merge_plot(roc_km_plot,tcga.module.km,ncol=1,nrow=2,labels = c("","G"))
ggsave("04_model/fig4.pdf",fig4,height =11,width =15 )

######GSE20685##########
GSE20685_model_data <- cbind(GSE20685.cli[, c("OS.time", "OS")],
                             t(GSE20685.exp[names(lan), GSE20685.cli$Samples]))
colnames(GSE20685_model_data) <- gsub('-', '_', colnames(GSE20685_model_data))


setdiff(names(lan),colnames(GSE20685_model_data))

risk.GSE20685=as.numeric(lan%*%as.matrix(t(GSE20685_model_data[GSE20685.cli$Samples,names(lan)])))
GSE20685.risktype.cli=data.frame(GSE20685.cli,Riskscore=risk.GSE20685)
GSE20685.point <- surv_cutpoint(GSE20685.risktype.cli, time = "OS.time", event = "OS",
                                variables = 'Riskscore')
GSE20685.point.cutoff <- as.numeric(summary(GSE20685.point)[1])
GSE20685.point.cutoff

GSE20685.risktype.cli$Risktype=ifelse(GSE20685.risktype.cli$Riskscore>GSE20685.point.cutoff,'High','Low')
GSE20685.roc.OS=ggplotTimeROC(GSE20685.risktype.cli$OS.time,
                              GSE20685.risktype.cli$OS,
                              GSE20685.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
GSE20685.roc.OS
GSE20685.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                      data = GSE20685.risktype.cli),
                          data=GSE20685.risktype.cli,
                          conf.int = T,pval = T,risk.table = T, 
                          fun = "pct",size = 1,surv.median.line = 'hv',
                          title='GSE20685',legend.title='Risktype',
                          legend.labs = c('High','Low'),
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,#ylab='Overall Survival(OS)',
                          legend.position='top',
                          ggtheme = theme_bw(base_size = 12))
GSE20685.km.OS=mg_merge_plot(GSE20685.km.OS$plot,GSE20685.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE20685.km.OS



#####
geo_expr <- my_mutiboxplot( t(GSE20685.exp[names(lan),GSE20685.risktype.cli$Samples]),
                group = GSE20685.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                ylab = 'Gene Expression',fill = 'Risktype',angle = 0,hjust = .5,notch = T,title = 'GSE20685')
fig5=mg_merge_plot(GSE20685.km.OS,GSE20685.roc.OS,geo_expr,ncol=3,nrow=1,labels = LETTERS[1:3],widths = c(1,1,1.5))
ggsave('04_model/fig5.pdf',fig5,height = 5,width = 15)



# library(ggbiplot)
# tcga.pca <- prcomp(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]), scale=T)
# pca_fig1 <- ggbiplot(tcga.pca, scale=1, groups = tcga.risktype.cli$Risktype,
#                      ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
#   scale_color_manual(values =risktype.col) + 
#   theme(legend.direction = 'horizontal', legend.position = 'top',text = element_text(family = 'Times')) +
#   xlab('PCA1') + ylab('PCA2')+ggtitle('TCGA')
# pca_fig1

# GSE20685.pca <- prcomp(t(GSE20685.exp[names(lan),GSE20685.risktype.cli$Samples]), scale=T)
# pca_fig2 <- ggbiplot(GSE20685.pca, scale=1, groups = GSE20685.risktype.cli$Risktype,
#                      ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
#   scale_color_manual(values = risktype.col) + 
#   theme(legend.direction = 'horizontal', legend.position = 'top',text = element_text(family = 'Times')) +
#   xlab('PCA1') + ylab('PCA2')+ggtitle('GSE20685')
# pca_fig2
# 
# 
# mg_merge_plot(pca_fig1,pca_fig2,labels = c('I','J'),align = 'h',common.legend = T)
# ggsave('04_model/pca.pdf',height = 5,width = 10)


#5.DAIVD####
############
dir.create("05_david")
p_cutoff <- 0.05
fc_cutoff <-log2(2) 
risk.limma <- mg_limma_DEG(tcga.exp,group = tcga.risktype.cli$Risktype,dlab = "Low",ulab = "High")
risk.limma$Summary
risk.degs=risk.limma$DEG[which(risk.limma$DEG$adj.P.Val<p_cutoff & abs(risk.limma$DEG$logFC)>fc_cutoff),]
risk.degs$ID=rownames(risk.degs)
dim(risk.degs)
head(risk.degs)
write.csv(risk.degs,file = "05_david/risk.degs.csv")



#
degs.david=read.table('05_david/DAVID.txt',sep = '\t',header = T)
head(degs.david)
degs.david=data.frame(Category=degs.david$Category,
                      ID=stringr::str_split_fixed(degs.david$Term,'~',2)[,1],
                      Term=stringr::str_split_fixed(degs.david$Term,'~',2)[,2],
                      Genes=degs.david$Genes,
                      adj_pval=degs.david$FDR)
dim(degs.david)

library(dplyr)

degs_top10 <- degs.david %>%
  group_by(Category) %>%                           # 
  slice_min(order_by = adj_pval, n = 10) %>%       # 
  ungroup()                                        #                                 # 







####circle_dat()
library(GOplot)
circ=circle_dat(degs_top10,risk.degs)
circ=circ[which(circ$adj_pval<0.01),]
head(circ)
table(circ$category)
# GOTERM_BP_DIRECT GOTERM_CC_DIRECT     KEGG_PATHWAY
# 77              180                7
circ_go=circ[which(circ$category=='GOTERM_BP_DIRECT'|circ$category=='GOTERM_CC_DIRECT'|circ$category=='GOTERM_MF_DIRECT'),]
circ_kegg=circ[which(circ$category=='KEGG_PATHWAY'),]
colnames(circ_kegg)
circ_kegg=data.frame(category=circ_kegg$category,
                     ID=stringr::str_split_fixed(circ_kegg$ID,':',2)[,1],
                     term=stringr::str_split_fixed(circ_kegg$ID,':',2)[,2],
                     circ_kegg[,4:8])
circ=rbind(circ_kegg,circ_go)
head(circ)
write.csv(circ,'05_david/DAVID_result.csv',row.names = F)


circ <- crbind2DataFrame(circ )
table(circ$category,circ$term)
fig6a=GOBubble(subset(circ,category=='GOTERM_BP_DIRECT'),colour =ggsci::pal_aaas()(10)[9],labels=1)#labels=c('cell division','G2/M transition of mitotic cell cycle','mitotic cell cycle','mitotic spindle assembly')
#ggsave("03_diff_genes/fig3c.pdf",height = 6,width = 6)
fig6b=GOCircle(subset(circ_go,category=='GOTERM_CC_DIRECT'))
circ_go_MF <- circ_go %>% subset(category=='GOTERM_MF_DIRECT')
fig6c=GOBubble(subset(circ_go_MF,category=='GOTERM_MF_DIRECT'),colour =ggsci::pal_aaas()(10)[7],labels =1)
fig6d=GOBubble(subset(circ,category=='KEGG_PATHWAY'),colour =ggsci::pal_aaas()(10)[5],labels =1)

pdf('05_david/Fig6.pdf',height = 15,width = 20)
mg_merge_plot(fig6a,fig6b,fig6c,fig6d,ncol=2,nrow=2,labels = LETTERS[1:4])
dev.off()
#06.#####################
dir.create('06_risktype.immu')
##ssgsea####
tcga.ssgsea=immu.ssgsea(exp =tcga.exp,isTCGA = T )
saveRDS(tcga.ssgsea,file ='07.risktype.immu/tcga.ssgsea.RDS')
pdf('06_risktype.immu/Fig7a.pdf',height = 5,width = 15)
mg_PlotMutiBoxplot(tcga.ssgsea[tcga.risktype.cli$Samples,1:22]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   ,group_cols = risktype.col
                   , test_method = c('kruskal.test','wilcox.test')[2]
                   , add = 'boxplot'
                   , ylab = 'Fraction')	
dev.off()

##estimate####
tcga_estimate <- immu_estimate(exp = tcga.exp)
head(tcga_estimate)
my_mutiboxplot(dat = tcga_estimate[tcga.risktype.cli$Samples,],group = tcga.risktype.cli$Risktype,group_cols = risktype.col)
#######
tcga_estimate <- tcga_estimate[,1:3]
tcga_estimate <- cbind(tcga_estimate[tcga.risktype.cli$Samples,],
                       tcga.risktype.cli)

colnames(tcga.risktype.cli)
tcga_StromalScore_cor <- cor_point(x = tcga_estimate$StromalScore,
                                   y = tcga_estimate$Riskscore,
                                   xlab = 'StromalScore',
                                   ylab = 'RiskScore',top_col='#D91656',right_col='#80C4E9'
)
tcga_StromalScore_cor

tcga_ImmuneScore_cor <- cor_point(x = tcga_estimate$ImmuneScore,
                                  y = tcga_estimate$Riskscore,
                                  xlab = 'ImmuneScore',
                                  ylab = 'RiskScore',
                                  top_col='#D91656',right_col='#80C4E9')
tcga_ImmuneScore_cor

tcga_ESTIMATEScore_cor <- cor_point(x = tcga_estimate$ESTIMATEScore,
                                    y = tcga_estimate$Riskscore,
                                    xlab = 'ESTIMATEScore',
                                    ylab = 'RiskScore',
                                    top_col='#D91656',right_col='#80C4E9')
tcga_ESTIMATEScore_cor

# ######

p7b<- cowplot::plot_grid(tcga_StromalScore_cor,
                         tcga_ImmuneScore_cor,
                         tcga_ESTIMATEScore_cor,
                         nrow = 3,labels = c("B",'C','D'))

ggsave("06_risktype.immu/p7b.pdf",p7b,height = 10,width = 5)
##29TMEsignature####
tme.type=readxl::read_excel('00_origin_datas/TME.geneSets.classification.PMID34019806.xlsx',sheet = "Sheet2")
tme.type=data.frame(tme.type,check.names = F,stringsAsFactors = F)
tme.type$Group=factor(tme.type$Group,levels=unique(tme.type$Group))
sort(tme.type$`Process/Signature`)
tme.anno=data.frame(group=tme.type$Group)
rownames(tme.anno)=tme.type$`Process/Signature`
####################TME geneset
tme.genesets=readxl::read_excel('00_origin_datas/29signatrue.PMID34019806.xlsx')
tme.genesets=data.frame(tme.genesets)
head(tme.genesets)
names(table(tme.genesets$Gene.signature))
sort(names(table(tme.genesets$Gene.signature)))
tme.genesets.list=split(x=tme.genesets,f=tme.genesets$Gene.signature)
tme.genesets.list=sapply(tme.genesets.list, function(x){subset(x,select='Gene',drop=TRUE)})
save(tme.genesets.list,file = "06_risktype.immu/tme.genesets.list.Rdata")
##TME geneset  ssgsea####
tcga.tme.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = tcga.exp
                                                ,genelist = tme.genesets.list)
save(tcga.tme.ssgsea,file='06_risktype.immu/tcga.tme.ssgsea.RData')
load('06_risktype.immu/tcga.tme.ssgsea.RData')
writeMatrix(tcga.tme.ssgsea,'06_risktype.immu/tcga.tme.ssgsea.txt')
# pdf('06_risktype.immu/Fig7b.pdf',height = 10,width = 6)
# mg_PlotMutiBoxplot(data = t(tcga.tme.ssgsea[rownames(tme.anno),tcga.risktype.cli$Samples])
#                    ,group = tcga.risktype.cli$Risktype
#                    ,test_method = 'wilcox.test'
#                    ,legend.pos = 'top'
#                    ,ylab = 'score',xangle = 0
#                    ,group_cols = risktype.col
#                    ,add = 'boxplot')+coord_flip()
# dev.off()

immu.dat.RS=cbind(tcga.risktype.cli$Riskscore,
                  tcga_model_data[tcga.risktype.cli$Samples,names(lan)],
                  t(tcga.tme.ssgsea[rownames(tme.anno),tcga.risktype.cli$Samples]))
colnames(immu.dat.RS)[1]='Riskcsore'

cor_res <- Hmisc::rcorr(as.matrix(immu.dat.RS),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0
p.dat=cor_res$P
p.dat[1:5,1:5]
p.dat<-ifelse(p.dat<0.0001,'****',ifelse(p.dat<0.001,'***', ifelse(p.dat<0.01,'**',ifelse(p.dat<0.05,'*',''))))

tme.col=pal_jama()(7)[4:7]
names(tme.col)=names(table(tme.anno$group))
library(ComplexHeatmap)
pdf('06_risktype.immu/Fig7c.pdf',height = 7,width = 10,onefile = F)
pheatmap(cor_res$r[-c(1:5),c(names(lan),'Riskcsore')],
         #scale = 'row',
         #border="white", # 
         color = circlize::colorRamp2(c(-1, 0, 1), c('#3B4992FF', 'white', '#EE0000FF')),
         #main="Heatmap", # 
         #legend_labels = c("-1", "0", "1", "cor\n"),
         #annotation_col = cli_anno,
         annotation_row = tme.anno,
         annotation_colors = list(group=tme.col),
         #cutree_cols = 2, #，
         #cutree_rows =2, # 
         #cellwidth = 6,cellheight = 5, # 
         display_numbers = p.dat[-c(1:5),c(c(names(lan)),'Riskcsore')], # 
         cluster_cols = F, # 
         cluster_rows = F,
         show_rownames = T, #
         show_colnames = T,
         #legend_breaks=c(-5,0,5),
         #angle_col = 45,
         gaps_row = c(8,8+7,8+7+12),
         fontsize_row =10, # 
         fontsize_col = 10)
dev.off()

writeMatrix(cor_res$r[-c(1:5),c(names(lan),'Riskcsore')],'06_risktype.immu/RS_cor_TME.txt')
save.image("BRCA_20241128.Rdata")
