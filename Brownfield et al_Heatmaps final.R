#Load libraries
library(gdata)
library(gplots)
library(reshape2)
library(reshape)
library(ggplot2)
library(FactoMineR)
library(plyr)
library(class)
library(doBy)
library(raster)
library(reshape2)
library(reshape)
library(class)

##################################Figure 1 Brownfield et al###################################################################################

###### Create heatmap for E18 epithelial cells
genes <- c("Fgfr2", "Fzd8", "Cd36", "Cd74", "Ngfrap1", "Lyz2", "Sftpc", "Ager", "S100a14", "Gapdh", "Actb", "Ppia")
countassembly <- read.csv("C:/R/analysis/Counts/E18ep/E18ep_fpkm_all_01272022.csv", header=T, check.names=FALSE, row.names = 1);
#The AT1<--BP--> order was previously established (See Treutlein et al, Nature 2014)
cellorder <- read.csv("C:/R/analysis/Counts/Cellorder/fig3order.csv", header=T);
cellorder <- cellorder$Order

count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- log2(count.mat);
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
ep_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=FALSE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######


###### Create heatmap for E18 endothelial cells
genes <- c("Fgf1", "Fgf3", "Fgf7", "Fgf10", "Fgf21", "Fgf22", "Pecam1", "Cdh5", "Col1a1", "Col1a2", "Gapdh", "Actb", "Ppia")
countassembly <- read.csv("C:/R/analysis/E18endo_fpkm_all.csv", header=T, check.names=FALSE, row.names=1);
cellorder <- colnames(countassembly)

#Remove cells that low/no expression of 2 or more "housekeeping genes" Gapdh, Actb, and Ppia
cellorder <- cellorder[cellorder != cellorder[grep('C82', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C85', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C93', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C95', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C96', cellorder)] ]

count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- log2(count.mat);
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
endo_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=TRUE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######


###### Create heatmap for E18 Mesenchymal cells
genes <- c("Fgf1", "Fgf3", "Fgf7", "Fgf10", "Fgf21", "Fgf22", "Pecam1", "Cdh5", "Col1a1", "Col1a2", "Gapdh", "Actb", "Ppia")
countassembly <- read.csv("C:/R/analysis/E18mes_fpkm_all.csv", header=T, check.names=FALSE, row.names = 1);
cellorder <- colnames(countassembly)

cellorder <- cellorder[cellorder != cellorder[grep('C01', cellorder)] ]

count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- log2(count.mat);
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
mes_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=TRUE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######


############################################################Figure 5######################################################################

###### Create Adult AT2 heatmap
genes <- c("Fgfr2", "Etv5", "Spry1", "Spry2", "Dusp1", "Dusp6", "Stat3", "Sftpc", "Lyz2", "Ppia", "Actb","Ubc")
countassembly <- read.csv("C:/R/analysis/Counts/Adult/AT2_062214_fpkm.csv", header=T, check.names=FALSE, row.names = 1);
cellorder <- colnames(countassembly)

#Sort out cells not expressing Ppia, Actb, and Ubc
cellorder <- cellorder[cellorder != cellorder[grep('C50', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C94', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C09', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C38', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C63', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C92', cellorder)] ]
cellorder <- cellorder[cellorder != cellorder[grep('C68', cellorder)] ]

count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- count.mat;
#count.trans <- log2(count.mat);
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
AT2_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=TRUE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######

###### Create Adult Mesenchyme heatmap
genes <- c("Fgf7", "Fgf10", "Col1a1", "Col1a2", "Mgp", "Ppia", "Actb", "Ubc", "Pdgfra")
countassembly <- read.csv("C:/R/analysis/Counts/Adult/Adult_mes_all_fpkm.csv", header=T, check.names=FALSE, row.names = 1);
cellorder <- colnames(countassembly)

#Sort out cells not expressing Ppia, Actb, and Ubc
cellorder <- cellorder[cellorder != cellorder[grep('C34', cellorder)] ]


count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- count.mat;
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
Adult-mes_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=TRUE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######


#############################################################Supplemental Figure 1 Brownfield et al########################################

###### Create Fgr2 isoform heatmap for E18 epithelial cells
genes <- c("NM_201601", "NM_010207")
countassembly <- read.csv("C:/R/analysis/Counts/E18ep/Excel/E18ep_fpkm_isoforms_05052015.csv", header=T, check.names=FALSE, row.names = 1);
cellorder <- read.csv("C:/R/analysis/Counts/Cellorder/fig3order.csv", header=T);
cellorder <- cellorder$Order

count.mat = matrix(, nrow=length(genes), ncol=length(cellorder));
row.names(count.mat) <- genes;
colnames(count.mat) <- cellorder;

j<- 1;
for (j in 1:length(cellorder))
{
  cell.col <- grep(cellorder[j], colnames(countassembly), fixed = TRUE);
  i<- 1;
  for (i in 1:length(genes))
  {
    gene.row <- grep(genes[i], row.names(countassembly), fixed = TRUE);
    count.mat[i,j] <- countassembly[gene.row[1], cell.col[1]];
  }
}

count.trans <- log2(count.mat);
count.trans[count.trans < 0] <- 0;
count.trans[is.na(count.trans)] <- 0;

Sys.sleep(0);
palette.breaks <- seq(0, 18.85423, 0.01);
color.palette = colorRampPalette(c("midnightblue","dodgerblue3", "white", "goldenrod1", "darkorange2"), space="Lab");
episo_counts_heatmap <- heatmap.2(count.trans, Rowv=FALSE, Colv=FALSE, key=TRUE, col=color.palette, breaks=palette.breaks, scale= "none", margins=c(6,15), trace= "none", density= "none", cexRow=0.8, cexCol=0.8);
Sys.sleep(0);
######
