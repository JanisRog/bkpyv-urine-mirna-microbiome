args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4) stop('Usage: Rscript run_count_models.R RECOUNT METADATA OUTPUT RLIB')
.libPaths(c(normalizePath(args[4]),.libPaths()))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
dir.create(args[3],recursive=TRUE,showWarnings=FALSE)
cts <- as.matrix(read.csv(file.path(args[1],'human_unique_counts.csv'),row.names=1,check.names=FALSE))
meta <- read.csv(args[2],stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(meta$sample),setequal(meta$sample,colnames(cts)),all(is.finite(cts)),all(cts>=0),all(cts==round(cts)))
rownames(meta)<-meta$sample;meta$group<-factor(meta$group,levels=c('control','treatment'))
cts<-cts[,meta$sample,drop=FALSE];storage.mode(cts)<-'integer'
yield<-colSums(cts);detected<-colSums(cts>=10)
sets<-list(main_ge2000=names(yield)[yield>=2000],all22=names(yield),ge1000=names(yield)[yield>=1000],old_rule=names(yield)[yield>=2000 & detected>=30])
sets$main_without_S3<-setdiff(sets$main_ge2000,'S3')
sets$main_without_S24<-setdiff(sets$main_ge2000,'S24')
summaries<-list();set.seed(20260922)
for(nm in names(sets)) {
 message('ANALYSIS ',nm)
 ss<-sets[[nm]];out<-file.path(args[3],nm);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 dat<-cts[,ss,drop=FALSE];keep<-rowSums(dat>=10)>=3;dat<-dat[keep,,drop=FALSE]
 cd<-meta[ss,,drop=FALSE]
 dds<-DESeqDataSetFromMatrix(dat,cd,~group)
 dds<-DESeq(dds,sfType='poscounts',fitType='local',minReplicatesForReplace=Inf,quiet=TRUE)
 res<-as.data.frame(results(dds,contrast=c('group','treatment','control'),independentFiltering=FALSE,alpha=.05))
 res$miRNA<-rownames(res);res$wald_lower95<-res$log2FoldChange-1.96*res$lfcSE;res$wald_upper95<-res$log2FoldChange+1.96*res$lfcSE
 res$max_Cooks<-apply(assays(dds)[['cooks']],1,max,na.rm=TRUE)
 res$dispersion<-dispersions(dds)
 res<-res[,c('miRNA',setdiff(names(res),'miRNA'))]
 write.csv(res,file.path(out,'deseq2_results.csv'),row.names=FALSE,na='')
 write.csv(data.frame(miRNA=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),file.path(out,'deseq2_normalized_counts.csv'),row.names=FALSE)
 sf<-data.frame(sample=ss,group=cd$group,raw_total=yield[ss],size_factor=sizeFactors(dds))
 write.csv(sf,file.path(out,'size_factors.csv'),row.names=FALSE)
 write.csv(as.data.frame(assays(dds)[['cooks']]),file.path(out,'Cooks_distances.csv'))
 vsd<-varianceStabilizingTransformation(dds,blind=FALSE)
 write.csv(data.frame(miRNA=rownames(vsd),assay(vsd),check.names=FALSE),file.path(out,'vst_counts.csv'),row.names=FALSE)
 pdf(file.path(out,'dispersion_diagnostic.pdf'));plotDispEsts(dds);dev.off()
 y<-DGEList(counts=dat,group=cd$group)
 # Full unique-mature library totals retained when restricting tested features.
 y$samples$lib.size<-yield[ss]
 y<-calcNormFactors(y,method='TMMwsp')
 design<-model.matrix(~group,cd)
 y<-estimateDisp(y,design)
 fit<-glmQLFit(y,design)
 test<-glmQLFTest(fit,coef=2)
 tab<-topTags(test,n=Inf,sort.by='none')$table
 tab$miRNA<-rownames(tab)
 write.csv(tab[,c('miRNA',setdiff(names(tab),'miRNA'))],file.path(out,'edger_TMMwsp_results.csv'),row.names=FALSE,na='')
 write.csv(data.frame(sample=ss,y$samples),file.path(out,'edger_library_factors.csv'),row.names=FALSE)
 summaries[[nm]]<-data.frame(strategy=nm,n=length(ss),cases=sum(cd$group=='treatment'),controls=sum(cd$group=='control'),features=nrow(dat),deseq2_tested=sum(is.finite(res$pvalue)),deseq2_FDR05=sum(res$padj<.05,na.rm=TRUE),edger_FDR05=sum(tab$FDR<.05,na.rm=TRUE),deseq2_min_FDR=min(res$padj,na.rm=TRUE),edger_min_FDR=min(tab$FDR,na.rm=TRUE))
 saveRDS(dds,file.path(out,'deseq2_fit.rds'))
}
write.csv(do.call(rbind,summaries),file.path(args[3],'model_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[3],'R_sessionInfo.txt'))
writeLines('Complete: all specified count-model analyses finished.',file.path(args[3],'SUCCESS'))
