library(pheatmap)
library(dplyr)
library(corrplot)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(annotatr)
library(rtracklayer)
library(biomaRt)
library(enrichplot)
library(clusterProfiler)
library(ChIPseeker)
library(ggplot2)

# analysis
squamous_085 <- ENCFF085QXI
names(squamous_085) <- c("chrom","chromStart","chromEnd","name","score","strand",
                         "signalValue","pValue","qValue","peak")
squamous_466 <- ENCFF466DPJ
names(squamous_466) <- c("chrom","chromStart","chromEnd","name","score","strand",
                         "signalValue","pValue","qValue","peak")
# select bins which are significant by score
squamous_085 <- squamous_085[squamous_085$pValue>50,]
squamous_466 <- squamous_466[squamous_466$pValue>50,]

squamous_085 <- arrange(squamous_085,chrom,chromStart)
squamous_466 <- arrange(squamous_466,chrom,chromStart)
write.table(squamous_466,"squamous_466.tab",quote = F,row.names = F,sep="\t")
write.table(squamous_085,"squamous_085.tab",quote = F,row.names = F,sep="\t")
# binning at BEAR portal with bedtools
# load bedtools
#module load BEDTools/2.25.0-foss-2018b
#tail -n+2 squamous_466>squamous_466_no_head.tab
#bedtools intersect -wa -wb -a binned_genome.bed -b squamous_466_no_head.tab -c>binned_squamous_466.tab

merged_bins <- H3K27ac_common_peaks%>%merge(binned_squamous_085,by=c("chrom","chromStart","chromEnd"))
merged_bins <- merged_bins%>%merge(binned_squamous_466,by=c("chrom","chromStart","chromEnd"))
corr_merged_bins <- merged_bins[,c("scc104","scc152","sqm_085","sqm_466")]
# row.names
row.names(merged_bins) <- paste(merged_bins$chrom,merged_bins$chromStart,sep = "_")
# heatmap
pheatmap(corr_merged_bins,angle_col = 45,display_numbers = T,
         method="ward",oma=c(6,6,6,6),fontsize = 14)
# plot cor(merged_bins) and corrplot of this correlation
corrplot.mixed(cor(corr_merged_bins),tl.cex=1.3,number.cex=1.3)

# annotate the significant bins
# create a GRanges object from the binned dataframe
key_locations <- merged_bins[c("chr1_209000000","chr3_189000000","chr9_127000000"),]
key_locations<- key_locations[,c("chrom","chromStart","chromEnd")]
# make a GRanges object
bins.GRanges <- makeGRangesFromDataFrame(key_locations)

# list the set of annotations needed
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries','hg19_genes_promoters',
           'hg19_genes_5UTRs','hg19_lncrna_gencode')
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we read in with the annotations
key_locations_annotated = annotate_regions(
  regions = bins.GRanges,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# Coerce to a data.frame
df_key_locations_annotated = data.frame(key_locations_annotated)
df_key_locations_annotated <- df_key_locations_annotated[complete.cases(df_key_locations_annotated),]
unique_genes <- unique(df_key_locations_annotated$annot.symbol)
write.table(unique_genes,"unique_genes.tab",row.names = F,quote = F,sep="\t")
# list of all genes from the all bins g.ranges object
all_genes <- all_bins_annotated@elementMetadata$annot$symbol
all_genes <- unique(all_genes)
# pathway analysis using biomart
# call biomart function to interconvert all gene names to entrez ids
gene_2_id_universe<-getBM(attributes = c('external_gene_name','entrezgene_id'),
                          filters = 'external_gene_name',
                          values = all_genes, # here all_data_gene_name is list of gene names
                          mart = ensembl)
gene_2_id_universe<-gene_2_id_universe[complete.cases(gene_2_id_universe),]
# extract gene_ids for GO enrichment
all_gene_ids<-gene_2_id_universe$entrezgene_id # here all_genes is list of entrez ids

# call the function to interconvert genes to entrez ids of unique genes
gene_2_id_unique<-getBM(attributes = c('external_gene_name','entrezgene_id'),
                        filters = 'external_gene_name',
                        values = unique_genes, 
                        mart = ensembl)
gene_2_id_unique<-gene_2_id_unique[complete.cases(gene_2_id_unique),]
# extract sample gene ids for GO enrichment
unique_gene_ids<-gene_2_id_unique$entrezgene_id

# create enrichGO object
ego_unique <- enrichGO(gene = unique_gene_ids,
                     universe      = all_gene_ids,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE)
# view the enrichGO dataframe
dim(ego_unique)
head(ego_unique[,c(1,2)],20)
# plot the data
barplot(ego_unique,font.size=30,showCategory = 5,x="p.adjust",ylab="pathways")
dotplot(ego_unique,showCategory=16,x="p.adjust") + ggtitle("pathways in bins")

# plot the chip seq data using chipseeker package
# create GRanges object from the narrow peak dataframe
# cell lines
scc104.GRanges <- makeGRangesFromDataFrame(SCC104_narrow)
scc152.GRanges <- makeGRangesFromDataFrame(SCC152_narrow)
# squamous epithelium
sq466.GRanges <- makeGRangesFromDataFrame(squamous_466)
sq085.GRanges <- makeGRangesFromDataFrame(squamous_085)
# create a coverage plot
covplot(scc104.GRanges)
covplot(scc152.GRanges)
# annotate the peak
peakAnno_scc104 <- annotatePeak(scc104.GRanges, tssRegion=c(-3000, 3000),
TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb="org.Hs.eg.db")

peakAnno_scc152 <- annotatePeak(scc152.GRanges, tssRegion=c(-3000, 3000),
                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb="org.Hs.eg.db")

peakAnno_sq466 <- annotatePeak(sq466.GRanges, tssRegion=c(-3000, 3000),
                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb="org.Hs.eg.db")

peakAnno_sq085 <- annotatePeak(sq085.GRanges, tssRegion=c(-3000, 3000),
                               TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb="org.Hs.eg.db")
# creating named list for joint analysis
peaklist <- list("scc104"=peakAnno_scc104,"scc152"=peakAnno_scc152,
                 "sqms.ep.466"=peakAnno_sq466,"sqms.ep.085"=peakAnno_sq085)

plotAnnoBar(peaklist)+ggplot2::theme(axis.text = element_text(size = 20))

