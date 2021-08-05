setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/XinPancreas")
library(cola)
library(scRNAseq)
data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/XinPancreas/XinPancreas_data.rds')
mat = as.matrix(assays(data)$rpkm)
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)[, c("age", "cell.type")]
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "XinPancreas_cola_rh.rds")

cola_report(rh, output = "XinPancreas_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'XinPancreas'")
