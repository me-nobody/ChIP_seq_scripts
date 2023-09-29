library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(GenomicRanges)
library(annotatr)
library(Homo.sapiens)


# create the mart object
ensembl <-useEnsembl(biomart="genes", dataset = "hsapiens_gene_ensembl")

# create a GRanges object from the binned dataframe
bins.GRanges_hpv_zero <- makeGRangesFromDataFrame(hpv_zero_H3K27ac)
# merge the binned dataframe with TxDb annotation
#annotated.bins <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), bins.GRanges_all)
# list the set of annotations needed
annots_genes="hg19_genes_cds"
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we read in with the annotations
hpv_zero_bins_annotated = annotate_regions(
  regions = bins.GRanges_hpv_zero,
  annotations = annotations_genes, # only genes annotated
  ignore.strand = TRUE,
  quiet = FALSE)
# Coerce to a data.frame
hpv_zero_bins_annotated_df <-  data.frame(hpv_zero_bins_annotated)
hpv_zero_bins_annotated_df <- hpv_zero_bins_annotated_df[complete.cases(hpv_zero_bins_annotated_df),]


# get the list of important genes in ChIP peak
hpv_zero_genes <- hpv_zero_bins_annotated_df$annot.symbol
hpv_zero_genes <- unique(hpv_zero_genes)

# create a merged dataframe from annotated bins and chip-seq bins
hpv_zero_bins_annotated_df <- hpv_zero_bins_annotated_df%>%left_join(hpv_zero_H3K27ac,by=c("chrom","chromStart","chromEnd"))
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
filtered_gene_ids<-gene_2_id_unique$entrezgene_id
# create enrichGO object
ego_chip <- enrichGO(gene = filtered_gene_ids,
                           universe      = all_gene_ids,
                           OrgDb         = org.Hs.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
# view the enrichGO dataframe
dim(ego_chip)
head(ego_chip[,c(1,2)],20)
# plot the data
barplot(ego_chip,showCategory = 16,x="p.adjust",ylab="pathways")
dotplot(ego_chip,showCategory=16,x="p.adjust") + ggtitle("overrepresented pathways in HPV infection")
