#!/usr/bin/env Rscript
# Human and viral mature-miRNA counts in the same prespecified count models.
# Viral reads were counted with a BK-only reference; the B1-3p mature sequence
# is also JCPyV J1-3p, so this model does not assign its viral species of origin.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Usage: Rscript run_joint_mirna_models.R ROOT OUT RLIB")
root <- normalizePath(args[1])
out <- args[2]
if (file.exists(out)) stop("Output directory already exists: ", out)
dir.create(out, recursive=TRUE)
.libPaths(c(normalizePath(args[3]), .libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))

human_all <- as.matrix(read.csv(file.path(root, "inputs/recount/human_unique_counts.csv"),
                                row.names=1, check.names=FALSE))
viral_all <- as.matrix(read.csv(file.path(root, "inputs/recount/viral_historical_counts.csv"),
                                row.names=1, check.names=FALSE))
stopifnot(nrow(viral_all)==2, all(human_all==round(human_all)),
          all(viral_all==round(viral_all)),
          identical(colnames(human_all), colnames(viral_all)))

for (strategy in c("main_ge2000", "all22", "ge1000", "old_rule",
                   "main_without_S3", "main_without_S24")) {
  target <- file.path(out, strategy)
  dir.create(target)
  sf <- read.csv(file.path(root, "results/count_models", strategy, "size_factors.csv"))
  ef <- read.csv(file.path(root, "results/count_models", strategy,
                           "edger_library_factors.csv"))
  ss <- sf$sample
  stopifnot(!anyDuplicated(ss), setequal(ss,ef$sample),
            all(ss %in% colnames(human_all)), all(ss %in% colnames(viral_all)),
            all(is.finite(sf$size_factor)), all(sf$size_factor>0))
  ef <- ef[match(ss, ef$sample),,drop=FALSE]
  human <- human_all[,ss,drop=FALSE]
  human <- human[rowSums(human >= 10) >= 3,,drop=FALSE]
  viral <- viral_all[,ss,drop=FALSE]
  stopifnot(nrow(human)==244)
  all_counts <- rbind(human, viral)
  storage.mode(all_counts) <- "integer"
  cd <- data.frame(group=factor(sf$group,levels=c("control","treatment")),
                   row.names=ss)
  dds <- DESeqDataSetFromMatrix(all_counts, cd, ~group)
  # Estimate exposure from host features alone, before adding the viral rows.
  sizeFactors(dds) <- setNames(sf$size_factor, ss)
  dds <- DESeq(dds, fitType="local", minReplicatesForReplace=Inf, quiet=TRUE)
  res <- as.data.frame(results(dds, contrast=c("group","treatment","control"),
                               independentFiltering=FALSE, alpha=.05))
  res$miRNA <- rownames(res)
  res$feature_type <- ifelse(res$miRNA %in% rownames(viral), "viral", "human")
  res$wald_lower95 <- res$log2FoldChange - 1.96*res$lfcSE
  res$wald_upper95 <- res$log2FoldChange + 1.96*res$lfcSE
  cooks <- assays(dds)[["cooks"]]
  res$max_Cooks <- apply(cooks,1,max,na.rm=TRUE)
  res <- res[,c("miRNA","feature_type",setdiff(names(res),c("miRNA","feature_type")))]
  write.csv(res,file.path(target,"deseq2_results.csv"),row.names=FALSE,na="")
  write.csv(as.data.frame(cooks),file.path(target,"Cooks_distances.csv"))
  write.csv(data.frame(sample=ss,group=sf$group,host_size_factor=sf$size_factor,
                       host_library_total=ef$lib.size,host_TMMwsp_factor=ef$norm.factors),
            file.path(target,"fixed_host_normalization.csv"),row.names=FALSE)
  write.csv(data.frame(miRNA=rownames(viral),counts(dds,normalized=TRUE)[rownames(viral),,
                                                                   drop=FALSE],check.names=FALSE),
            file.path(target,"viral_normalized_counts.csv"),row.names=FALSE)

  # Diagnostic only: this intentionally disables DESeq2's outlier safeguard.
  # Never use these p-values as the primary result when Cook's filtering acts.
  unfiltered <- as.data.frame(results(dds, contrast=c("group","treatment","control"),
                                      independentFiltering=FALSE,cooksCutoff=FALSE))
  unfiltered$miRNA <- rownames(unfiltered)
  write.csv(unfiltered[rownames(viral),c("miRNA","log2FoldChange","lfcSE","pvalue","padj")],
            file.path(target,"viral_no_Cooks_filter_DIAGNOSTIC_ONLY.csv"),row.names=FALSE)

  # Same host-only library scaling as the existing edgeR human sensitivity.
  y <- DGEList(counts=all_counts,group=cd$group)
  y$samples$lib.size <- ef$lib.size
  y$samples$norm.factors <- ef$norm.factors
  design <- model.matrix(~group,cd)
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  test <- glmQLFTest(fit,coef=2)
  tab <- topTags(test,n=Inf,sort.by="none")$table
  tab$miRNA <- rownames(tab)
  tab$feature_type <- ifelse(tab$miRNA %in% rownames(viral),"viral","human")
  write.csv(tab[,c("miRNA","feature_type",setdiff(names(tab),c("miRNA","feature_type")))],
            file.path(target,"edger_robust_QL_results.csv"),row.names=FALSE)
  saveRDS(dds,file.path(target,"deseq2_fit.rds"))
  message(strategy, ": ", sum(cd$group=="treatment"), " cases, ",
          sum(cd$group=="control"), " controls")
  print(res[res$feature_type=="viral",
            c("miRNA","log2FoldChange","pvalue","padj","max_Cooks")])
  print(tab[rownames(viral),c("logFC","PValue","FDR")])
}
# Influence diagnostic for the viral 5p Cook's-distance outlier. This is not a
# newly selected primary subset and must not be used to discard S7 post hoc.
main <- readRDS(file.path(out,"main_ge2000","deseq2_fit.rds"))
ss <- setdiff(colnames(main),"S7")
influence <- DESeqDataSetFromMatrix(counts(main)[,ss,drop=FALSE],
                                    as.data.frame(colData(main))[ss,,drop=FALSE],~group)
sizeFactors(influence) <- sizeFactors(main)[ss]
influence <- DESeq(influence,fitType="local",minReplicatesForReplace=Inf,quiet=TRUE)
ir <- as.data.frame(results(influence,contrast=c("group","treatment","control"),
                            independentFiltering=FALSE))
ir$miRNA <- rownames(ir)
write.csv(ir[rownames(viral_all),c("miRNA","log2FoldChange","lfcSE","pvalue","padj")],
          file.path(out,"main_without_S7_DIAGNOSTIC_ONLY.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
writeLines("Completed joint models; see per-subset outputs.",file.path(out,"SUCCESS"))
