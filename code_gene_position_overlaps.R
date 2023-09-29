BiocManager::install("Homo.sapiens")
library(Homo.sapiens)
library(GenomicRanges)
library(annotatr)
library(rtracklayer)
library(pheatmap)
library(tidyverse)
library(RCircos)
library(dplyr)
library(corrplot)
# corplot between values
corrplot.mixed(cor(m1),main="ChIP peak distribution")
# heatmap of the common bins
pheatmap(common_chip_bins,angle_col = 45,cluster_rows = T,cluster_cols = F,display_numbers = T,oma=c(6,6,6,6))
# calculate p-values for the ratio in chip_counts
chip_counts$ratio <- chip_counts$scc104/chip_counts$scc152
# convert bins to factor
chip_counts$bins <- 1:nrow(chip_counts)
# to handle NaN values which corrupt calculation
chip_counts[is.na(chip_counts) | chip_counts=="Inf"] = NA
# linear model for ratio vs bins  
lmod <- lm(ratio~bins,chip_counts,na.action = na.omit)
# test for differences between ratios
anova(lm(ratio~bins,chip_counts))
pvals <- summary(lmod)$coef[,4]
sigp <- coef( lmod )[ pvals < 0.05 ]
# ******Annotate the positions***********
# separate chr and position
m1<- m1%>%separate(chr,into=c("chr","start"),sep="_")
# 1MB bin so add 1MB to create stop
m1$start <- as.integer(m1$start)
m1 <- m1%>%mutate(stop=m1$start+1000000)
# re-order columns and provide proper names
names(m1) <- c("chrom","chromStart","chromEnd","scc104","scc152","ratio")
# create a GRanges object from the binned dataframe
bins.GRanges <- makeGRangesFromDataFrame(m1)
# update the TxDb for both strands -not needed here
TxDb <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene, single.strand.genes.only = FALSE)
# merge the binned dataframe with TxDb annotation
annotated.bins <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), bins.GRanges)
# above is not exactly in the format needed- we shall use annotatr
# list the set of annotations needed
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries','hg19_genes_promoters',
           'hg19_genes_5UTRs','hg19_lncrna_gencode')
annots_genes="hg19_genes_cds"
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
  regions = bins.GRanges,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# Coerce to a data.frame
df_dm_annotated = data.frame(dm_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))

annotations_genes = build_annotations(genome = 'hg19', annotations = annots_genes)

dm_annotated_genes = annotate_regions(
  regions = bins.GRanges,
  annotations = annotations_genes,
  ignore.strand = TRUE,
  quiet = FALSE)

df_dm_annotated_genes = data.frame(dm_annotated_genes)
annot_genes$track <- rep(80,nrow(annot_genes))
# plot the significant bins by RCircos
# RCircos ideogram
data("UCSC.HG19.Human.CytoBandIdeogram")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info,chr.exclude = NULL,
                            tracks.inside = 10,tracks.outside = 0)
rcircos_params <- RCircos.Get.Plot.Parameters()  
rcircos_params$base.per.unit <- 200000
rcircos_params$text.size <- 0.7
rcircos_params$heatmap.width <- 100
rcircos_params$hist.width <- 100
RCircos.Reset.Plot.Parameters(rcircos_params)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
# RCircos tile
data("annot_genes")
RCircos.Heatmap.Plot(annot_genes,track.num=8,side="in",data.col=5)
RCircos.Gene.Name.Plot(annot_genes,                 
                         name.col=4,track.num=3, side="in")

# skin epidermis data with score>500
skin_epidermis_ENCFF003TGC<- skin_epidermis_ENCFF003TGC[skin_epidermis_ENCFF003TGC$score>500,]

skin_epidermis_ENCFF423CQD<- skin_epidermis_ENCFF423CQD[skin_epidermis_ENCFF423CQD$score>500,]

# merge skin epidermis with the binned_peaks df
binned_peaks%>%inner_join(binned_skin_epiderm_ENCFF003TGC,by=)
# dataframe plotted
corrplot.mixed(cor(corr_binned_enriched))
corrplot.mixed(cor(corr_binned_peaks),main="ChIP peak distribution")
# ******Annotate the positions***********
# separate chr and position
valid_chip_count <- chip_counts[(chip_counts$scc104>=1) & (chip_counts$scc152>=1),]
valid_chip_count$pos <- row.names(valid_chip_count)
# partition the pos  column
valid_chip_count<- valid_chip_count%>%separate(pos,into=c("chr","start"),sep="_")
# 1MB bin so add 1MB to create stop
valid_chip_count$start <- as.integer(valid_chip_count$start)
valid_chip_count <- valid_chip_count%>%mutate(stop=valid_chip_count$start+1000000)
# counts same as those for cell lines with HPV
valid_chip_count <- merge(valid_chip_count,binned_peaks,by=c("chrom","chromStart","chromEnd"))
valid_chip_count <- valid_chip_count[,c("chrom","chromStart","chromEnd","scc104.x","scc152.x","003tgc","423cqd")]
# corplot between values
dummy <- dummy%>%filter(scc104.x>1,scc152.x>1,`003tgc`>1,`423cqd`>1)
# corrplot with bigger text and numbers
corrplot.mixed(cor(m2),tl.cex=2,number.cex=2)