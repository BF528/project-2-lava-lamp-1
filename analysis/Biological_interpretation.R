library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
##7.1
#READ FILES 
fpkm<- read.csv('/project/bf528/project_2/data/fpkm_matrix.csv',sep='\t',header = T)
p0_1<-read.csv('/projectnb/bf528/users/lava-lamp/project-2-lava-lamp-1/samples/cuffdiff_out/genes.fpkm_tracking',sep='\t',header=T)%>%
  select(gene_short_name,P0_FPKM)
#Aggregate by gene short name
p0_1_ag<-aggregate(P0_FPKM~gene_short_name, data = p0_1, FUN = mean, na.rm = TRUE)
id_name<-read.csv('/projectnb/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking',sep='\t',header = T)
id_name<-id_name%>% select(tracking_id, gene_short_name)
fpkm<-merge(fpkm,id_name,by='tracking_id')
#Aggregate by gene short name
fpkm_aggregate<-aggregate(cbind(Ad_1_FPKM,Ad_2_FPKM,P0_2_FPKM,P4_1_FPKM,P4_2_FPKM,P7_1_FPKM,P7_2_FPKM)~gene_short_name, data = fpkm, FUN = mean, na.rm = TRUE)
fpkm_all<-merge(fpkm_aggregate,p0_1_ag,by='gene_short_name')
#Take mean between two duplicates
fpkm_all$P0<- rowMeans(fpkm_all%>% select(P0_FPKM,P0_2_FPKM), na.rm=TRUE)
fpkm_all$P4<- rowMeans(fpkm_all%>% select(P4_1_FPKM,P4_2_FPKM), na.rm=TRUE)
fpkm_all$P7<- rowMeans(fpkm_all%>% select(P7_1_FPKM,P7_2_FPKM), na.rm=TRUE)
fpkm_all$Ad<- rowMeans(fpkm_all%>% select(Ad_1_FPKM,Ad_2_FPKM), na.rm=TRUE)
#Create the final fpkm data frame to plot.
fpkm_final<-fpkm_all%>% select(gene_short_name,P0,P4,P7,Ad)
sarc<-list('Pdlim5','Pygm','Myoz2','Des','Csrp3','Tcap','Cryab')
mito<-list('Mpc1','Prdx3','Acat1','Echs1','Slc25a11','Phyh')
cc<-list('Cdc7','E2f8','Cdk7','Cdc26','Cdc6','E2f1','Cdc27','Bora','Cdc45','Rad51','Aurkb','Cdc23')
x<-list('P0','P4','P7','AD')

temp_1 <- filter(fpkm_final, gene_short_name %in% sarc)%>%
  pivot_longer(cols = c(P0,P4,P7,Ad), names_to = "name", values_to = "value")
sarc_plot<-ggplot(temp_1,aes(x =factor(name,level=x), y = value,
                             color=gene_short_name,group=gene_short_name)) + 
  geom_line()+
  geom_point()+
  labs(title = "Sarcomere", x='in vivo Maturation', y="FPKM")
temp_2 <- filter(fpkm_final, gene_short_name %in% mito)%>%
  pivot_longer(cols = c(P0,P4,P7,Ad), names_to = "name", values_to = "value")
mito_plot<-ggplot(temp_2,aes(x =factor(name,level=x), y = value,color=gene_short_name,group=gene_short_name)) + 
  geom_line()+
  geom_point()+
  labs(title = "Mitochondria", x='in vivo Maturation', y="FPKM")
temp_3 <- filter(fpkm_final, gene_short_name %in% cc)%>%
  pivot_longer(cols = c(P0,P4,P7,Ad), names_to = "name", values_to = "value")
cc_plot<-ggplot(temp_3,aes(x =factor(name,level=x), y = value,color=gene_short_name,group=gene_short_name)) + 
  geom_line()+
  geom_point()+
  labs(title = "Cell Cycle", x='in vivo Maturation', y="FPKM")
grid.arrange(grobs = list(sarc_plot,mito_plot, cc_plot))

##7.3
diff_exp<-read.csv('/projectnb/bf528/users/lava-lamp/project-2-lava-lamp-1/samples/cuffdiff_out/gene_exp.diff',sep='\t',header=T)
#GET TOP 1000 SIGINIFICANT RESULTS
diff_exp<- diff_exp[diff_exp$significant == 'yes', ]
diff_exp<- diff_exp[order(diff_exp$q_value)[1:1000],]
#Get significant genes
significant_genes<-diff_exp$gene
fpkm_significant<-filter(fpkm_final,gene_short_name%in% significant_genes)
row.names(fpkm_significant)<-fpkm_significant$gene_short_name
heatmap_colors<-colorRampPalette(c("blue","black", "yellow"), interpolate = "spline")
gene_heatmap <- heatmap(as.matrix(fpkm_significant[, 2:5]),
                      scale = "row", 
                      col = heatmap_colors(10),
                      xlab = "Sample name", 
                      ylab = "Genes")
legend(x="right", legend=c('<10%','10%-20%','20%-30%','30%-40%','40%-50%','50%-60%',
                           '60%-70%','70%-80%','80%-90%','90%-100%'),fill=heatmap_colors(10))

