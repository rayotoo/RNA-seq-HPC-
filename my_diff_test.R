library(edgeR)

base.dir = path.expand("~/rnaseq")
base.dir

# meta data
meta_file = paste(base.dir, "data/ibs_class/meta_data.txt", sep="/")
sample_list <- read.table(meta_file, header=TRUE, stringsAsFactors=FALSE)
samples <- sample_list[order(sample_list$Sample_ID),]

# count data
counts_file = paste(base.dir, "out/HTSeq/CountMerged.txt", sep="/")
counts <- read.table(counts_file, header=TRUE, sep ="\t", stringsAsFactors = FALSE)
counts <- counts[,c("id",samples$Sample_ID)]

row.names(counts) <- counts$id
counts$id <- NULL
colnames(counts) <- samples$Sample_ID

# gene filtering
keep <- rowSums(counts >= 5) >= 3
counts <- counts[keep,]

# TMM normalization in edgeR
# edgeR normalizes by total count
# edgeR is concerned with differential expression analysis rather than with the quantification of expression levels. 

# Function definition
tmm <- function(counts, groups=NA){
  d <- DGEList(counts=counts, group=factor(groups))
  d <- calcNormFactors(d, method="TMM") 
  return(d)
}

# Functional call
dgList <- tmm(counts, samples$Group)

# After filtering, it is a good idea to reset the library sizes
# the lib.size and norm.factors are multiplied together to act as the effective library size
dgList$samples$lib.size <- colSums(dgList$counts)

#===============================================#
# Differential Expression Analysis using t-test #
#===============================================#

# CPM, Count Per Million
# log2 - transform
d1 <- cpm(dgList, normalized.lib.sizes=TRUE,log=TRUE)

# Control sample ids
ctrl <- samples[samples$Group == "control","Sample_ID"]

# case sample ids
case <- samples[samples$Group == "case","Sample_ID"]

# Function
test_fun <- function(x) {
  t.test(x[ctrl], x[case], alternative = c("two.sided"), paired = FALSE, var.equal=TRUE)$p.value
}

# t-test
all.p.genes = apply(counts, 1, test_fun)

# FDR, 
fdr.value = p.adjust(all.p.genes, method="fdr")

# Fold Change calculation
ctrl_mean <- rowMeans(d1[,ctrl])
grp_mean <- rowMeans(d1[,case])
fc.value <-  grp_mean-ctrl_mean

# Combine results 
fc.plus.pvals = cbind(fc.value,all.p.genes, fdr.value)

# Save results
write.table(fc.plus.pvals,file=paste(base.dir, "/star_htseq_ttest_results.txt", sep=""), row.names=T, col.names=NA, quote=F, sep="\t")

#==============================================#
# Differential Expression Analysis using edgeR #
#==============================================#
# print sample details
samples

samples[samples$Group == "control","dm"] <- "1"
samples[samples$Group == "case","dm"] <- "2"
samples$dm

design <- model.matrix(~0+samples$dm)
colnames(design) <- levels(samples$dm)
design

dg <- estimateDisp(dgList,design)
et <- exactTest(dg, pair=2:1)
results_edgeR <- topTags(et, n = nrow(dgList$counts), sort.by = "none")

# Save results
write.table(results_edgeR,file=paste(base.dir, "/star_htseq_edgeR_results.txt", sep=""), row.names=T, col.names=NA, quote=F, sep="\t")
