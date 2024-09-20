
## data

- Treatments (12): EGFP, E2D1_aCS3, E2D1_aCS3_V33R, E2D1_aCS3_C85A, E2D1_aCS3_C85A/V33R, E2D1_aCS3_F62A, E2D1_aNSa5, E2B_aCS3, VHL_aCS3, VHL_aCS3_V33R, siRNA (Control), siRNA (SHP2)
- Replicates each (4)

[1] PE203_2021_GP_bioPRO_sample1_reinj.raw.PG.Quantity
[2] PE203_2021_GP_bioPRO_sample2.raw.PG.Quantity
[3] PE203_2021_GP_bioPRO_sample3.raw.PG.Quantity
[4] PE203_2021_GP_bioPRO_sample4.raw.PG.Quantity
[5] PE203_2021_GP_bioPRO_sample5.raw.PG.Quantity
[6] PE203_2021_GP_bioPRO_sample6.raw.PG.Quantity
[7] PE203_2021_GP_bioPRO_sample7.raw.PG.Quantity
[8] PE203_2021_GP_bioPRO_sample8.raw.PG.Quantity
[9] PE203_2021_GP_bioPRO_sample9.raw.PG.Quantity
[10] PE203_2021_GP_bioPRO_sample10.raw.PG.Quantity
[11] PE203_2021_GP_bioPRO_sample11_reinj.raw.PG.Quantity
[12] PE203_2021_GP_bioPRO_sample12.raw.PG.Quantity
[13] PE203_2021_GP_bioPRO_sample13_reinj.raw.PG.Quantity
[14] PE203_2021_GP_bioPRO_sample14.raw.PG.Quantity
[15] PE203_2021_GP_bioPRO_sample15.raw.PG.Quantity
[16] PE203_2021_GP_bioPRO_sample16.raw.PG.Quantity
[17] PE203_2021_GP_bioPRO_sample17.raw.PG.Quantity
[18] PE203_2021_GP_bioPRO_sample18_reinj.raw.PG.Quantity
[19] PE203_2021_GP_bioPRO_sample19.raw.PG.Quantity
[20] PE203_2021_GP_bioPRO_sample20.raw.PG.Quantity
[21] PE203_2021_GP_bioPRO_sample21.raw.PG.Quantity
[22] PE203_2021_GP_bioPRO_sample22.raw.PG.Quantity
[23] PE203_2021_GP_bioPRO_sample23.raw.PG.Quantity
[24] PE203_2021_GP_bioPRO_sample24.raw.PG.Quantity
[25] PE203_2021_GP_bioPRO_sample25.raw.PG.Quantity
[26] PE203_2021_GP_bioPRO_sample26.raw.PG.Quantity
[27] PE203_2021_GP_bioPRO_sample27.raw.PG.Quantity
[28] PE203_2021_GP_bioPRO_sample28.raw.PG.Quantity
[29] PE203_2021_GP_bioPRO_sample29.raw.PG.Quantity
[30] PE203_2021_GP_bioPRO_sample30.raw.PG.Quantity
[31] PE203_2021_GP_bioPRO_sample31.raw.PG.Quantity
[32] PE203_2021_GP_bioPRO_sample32.raw.PG.Quantity
[33] PE203_2021_GP_bioPRO_sample33.raw.PG.Quantity
[34] PE203_2021_GP_bioPRO_sample34.raw.PG.Quantity
[35] PE203_2021_GP_bioPRO_sample35.raw.PG.Quantity
[36] PE203_2021_GP_bioPRO_sample36.raw.PG.Quantity
[37] PE203_2021_GP_bioPRO_sample37.raw.PG.Quantity
[38] PE203_2021_GP_bioPRO_sample38.raw.PG.Quantity
[39] PE203_2021_GP_bioPRO_sample39.raw.PG.Quantity
[40] PE203_2021_GP_bioPRO_sample40.raw.PG.Quantity
[41] PE203_2021_GP_bioPRO_sample41.raw.PG.Quantity
[42] PE203_2021_GP_bioPRO_sample42.raw.PG.Quantity
[43] PE203_2021_GP_bioPRO_sample43.raw.PG.Quantity
[44] PE203_2021_GP_bioPRO_sample44.raw.PG.Quantity
[45] PE203_2021_GP_bioPRO_sample45.raw.PG.Quantity
[46] PE203_2021_GP_bioPRO_sample46.raw.PG.Quantity
[47] PE203_2021_GP_bioPRO_sample47.raw.PG.Quantity
[48] PE203_2021_GP_bioPRO_sample48.raw.PG.Quantity

- Comparisons:
  - E2D1_aCS3 vs EGFP
  - E2D1_aCS3_V33R vs EGFP
  - E2D1_aCS3_C85A vs EGFP
  - E2D1_aCS3_C85A/V33R vs EGFP
  - E2D1_aCS3_F62A vs EGFP
  - E2D1_aNSa5 vs EGFP
  - E2B_aCS3 vs EGFP
  - VHL_aCS3 vs EGFP
  - VHL_aCS3_V33R vs EGFP
  - EGFP vs siRNA (Control)
  - E2D1_aCS3 vs siRNA (Control)
  - E2D1_aNSa5 vs siRNA (Control)
  - E2B_aCS3 vs siRNA (Control)
  - VHL_aCS3 vs siRNA (Control)
  - siRNA (SHP2) vs siRNA (Control)

48 samples were all processed on one plate and randomized for injection on the mass spec to avoid any bias from changing instrument performance during acquisition.

The DIA runs were analyzed in Spectronaut which has a local normalization option based on a RT dependent local regression model described by Callister et al 2006.



## pxanalytics and limma analysis

```{r}
cd ~/data
R

library(devtools)
load_all("~/R/3.6.0/pxanalytics")
library(tidyverse)
#library(magrittr)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)

# Load and select data
data <- fread("20211207_134123_PE203_2021_bioPRO_normalized_proteins.txt")
data <- data[,c(2,51:98)]
colnames(data) <- c("genes", "EGFP_1", "EGFP_2", "EGFP_3", "EGFP_4", "E2D1_aCS3_1", "E2D1_aCS3_2", "E2D1_aCS3_3", "E2D1_aCS3_4", "E2D1_aNSa5_1", "E2D1_aNSa5_2", "E2D1_aNSa5_3", "E2D1_aNSa5_4", "E2D1_aCS3_V33R_1", "E2D1_aCS3_V33R_2", "E2D1_aCS3_V33R_3", "E2D1_aCS3_V33R_4", "E2D1_aCS3_C85A_1", "E2D1_aCS3_C85A_2", "E2D1_aCS3_C85A_3", "E2D1_aCS3_C85A_4", "E2D1_aCS3_C85A_V33R_1", "E2D1_aCS3_C85A_V33R_2", "E2D1_aCS3_C85A_V33R_3", "E2D1_aCS3_C85A_V33R_4", "VHL_aCS3_V33R_1", "VHL_aCS3_V33R_2", "VHL_aCS3_V33R_3", "VHL_aCS3_V33R_4", "siRNA_control_1", "siRNA_control_2", "siRNA_control_3", "siRNA_control_4", "siRNA_shp2_1", "siRNA_shp2_2", "siRNA_shp2_3", "siRNA_shp2_4", "E2D1_aCS3_F62A_1", "E2D1_aCS3_F62A_2", "E2D1_aCS3_F62A_3", "E2D1_aCS3_F62A_4", "E2B_aCS3_1", "E2B_aCS3_2", "E2B_aCS3_3", "E2B_aCS3_4", "VHL_aCS3_1", "VHL_aCS3_2", "VHL_aCS3_3", "VHL_aCS3_4")
data <- unique(data, by = "genes") # remove duplicated rows # distinct function
nrow(data) # 7657
data <- data.frame(data)
rownames(data) <- data$genes

# Assay data
assay <- na.omit(data[,2:49])
nrow(data) # 7657
nrow(assay) # 7657
length(rownames(assay)) # 7657
length(unique(rownames(assay))) # 7657


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_")), "Replicate" = sapply(colnames(assay), function (x) paste(tail(unlist(strsplit(x, "_")), 1), collapse="_")), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
px_norm <- px_raw %>% pxnormalise(method = "total")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "E2D1_aCS3_vs_EGFP" = E2D1_aCS3-EGFP,
  "E2D1_aCS3_V33R_vs_EGFP" = E2D1_aCS3_V33R-EGFP,
  "E2D1_aCS3_C85A_vs_EGFP" = E2D1_aCS3_C85A-EGFP,
  "E2D1_aCS3_C85A_V33R_vs_EGFP" = E2D1_aCS3_C85A_V33R-EGFP,
  "E2D1_aCS3_F62A_vs_EGFP" = E2D1_aCS3_F62A-EGFP,
  "E2D1_aNSa5_vs_EGFP" = E2D1_aNSa5-EGFP,
  "E2B_aCS3_vs_EGFP" = E2B_aCS3-EGFP,
  "VHL_aCS3_vs_EGFP" = VHL_aCS3-EGFP,
  "VHL_aCS3_V33R_vs_EGFP" = VHL_aCS3_V33R-EGFP,
  "EGFP_vs_siRNA_control" = EGFP-siRNA_control,
  "E2D1_aCS3_vs_siRNA_control" = E2D1_aCS3-siRNA_control,
  "E2D1_aNSa5_vs_siRNA_control" = E2D1_aNSa5-siRNA_control,
  "E2B_aCS3_vs_siRNA_control" = E2B_aCS3-siRNA_control,
  "VHL_aCS3_vs_siRNA_control" = VHL_aCS3-siRNA_control,
  "siRNA_shp2_vs_siRNA_control" = siRNA_shp2-siRNA_control,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# stats, PTPN11, tables, volcanos

for (c in colnames(contrast.matrix)){
  detable <- data.table(topTable(fit2, coef=c, adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)
  cat(c, "\t", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05]), "\t", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05]), "\n")
  print(detable[rn == "PTPN11"])
  # write tables
  ## down
  write.table(detable[logFC < -0.5 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], paste("../tables/20220112/", c, ".down.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  ## up
  write.table(detable[logFC > 0.5 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], paste("../tables/20220112/", c, ".up.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  cat("****************************************\n")
  # volcano
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = head(detable[logFC < -0.5 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
  geom_text_repel(data = head(detable[logFC > 0.5 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)
  ggsave(paste('../figures/20220112/volcano_', c, '.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
}

E2D1_aCS3_vs_EGFP 	 93 	 143
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.165073 17.13661 -8.166509 6.350491e-10 6.096292e-07 12.47788
****************************************
E2D1_aCS3_V33R_vs_EGFP 	 124 	 407
       rn     logFC  AveExpr         t     P.Value  adj.P.Val         B
1: PTPN11 -0.433324 17.13661 -3.037359 0.004274164 0.04357826 -2.335119
****************************************
E2D1_aCS3_C85A_vs_EGFP 	 144 	 241
       rn       logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 -0.08704418 17.13661 -0.610131 0.545367 0.7893251 -6.395549
****************************************
E2D1_aCS3_C85A_V33R_vs_EGFP 	 243 	 284
       rn        logFC  AveExpr           t   P.Value adj.P.Val         B
1: PTPN11 -0.001972656 17.13661 -0.01382722 0.9890394 0.9958144 -6.590937
****************************************
E2D1_aCS3_F62A_vs_EGFP 	 228 	 88
       rn      logFC  AveExpr          t  P.Value adj.P.Val         B
1: PTPN11 -0.1318922 17.13661 -0.9244906 0.361001 0.6506443 -5.781113
****************************************
E2D1_aNSa5_vs_EGFP 	 110 	 57
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.8227283 17.13661 -5.766865 1.146419e-06 0.0005486334 5.425834
****************************************
E2B_aCS3_vs_EGFP 	 242 	 147
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.890644 17.13661 -13.25236 6.565184e-16 2.513481e-12 24.69009
****************************************
VHL_aCS3_vs_EGFP 	 372 	 45
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.947983 17.13661 -13.65427 2.534462e-16 9.703189e-13 25.46349
****************************************
VHL_aCS3_V33R_vs_EGFP 	 48 	 45
       rn       logFC  AveExpr          t   P.Value adj.P.Val        B
1: PTPN11 -0.04845112 17.13661 -0.3396153 0.7359958 0.9208365 -6.15588
****************************************
EGFP_vs_siRNA_control 	 265 	 221
       rn     logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 0.0641864 17.13661 0.4499107 0.655303 0.8252947 -6.213415
****************************************
E2D1_aCS3_vs_siRNA_control 	 654 	 708
       rn     logFC  AveExpr         t     P.Value    adj.P.Val        B
1: PTPN11 -1.100886 17.13661 -7.716598 2.50521e-09 5.480684e-07 11.30448
****************************************
E2D1_aNSa5_vs_siRNA_control 	 988 	 522
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.7585419 17.13661 -5.316954 4.771653e-06 0.0001574851 3.965503
****************************************
E2B_aCS3_vs_siRNA_control 	 87 	 119
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.826458 17.13661 -12.80245 1.948807e-15 7.461008e-12 24.14485
****************************************
VHL_aCS3_vs_siRNA_control 	 174 	 65
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.883797 17.13661 -13.20436 7.364813e-16 1.879746e-12 25.15552
****************************************
siRNA_shp2_vs_siRNA_control 	 36 	 34
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.426508 17.13661 -9.999022 3.027922e-12 7.843096e-09 16.91789
****************************************


# jitter plot SHP2
data_jitter <- data.table(ID = colnames(data.table(data)[genes=="PTPN11"][,2:49]), Intensity = as.numeric(data.table(data)[genes=="PTPN11"][,2:49]))
data_jitter[, Sample := sapply(ID, function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_"))]
data_jitter[, Sample := factor(Sample, levels = c("EGFP", "E2D1_aCS3", "E2D1_aNSa5", "E2D1_aCS3_C85A", "E2D1_aCS3_V33R", "E2D1_aCS3_C85A_V33R", "E2D1_aCS3_F62A", "E2B_aCS3", "VHL_aCS3", "VHL_aCS3_V33R", "siRNA_control", "siRNA_shp2"))]

gg <- ggplot(data_jitter, aes(x=Sample, y=Intensity)) +
geom_point() +
ylab("SHP2 normalised intensity") +
xlab("") +
theme_classic() +
theme(axis.title = element_text(size=14), axis.text = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('../figures/20220112/jitter_PTPN11.pdf', width = 6, height = 6, useDingbats = FALSE)

gg <- ggplot(data_jitter, aes(x=Sample, y=log2(Intensity))) +
geom_point() +
ylab("log2(SHP2 normalised intensity)") +
xlab("") +
theme_classic() +
theme(axis.title = element_text(size=14), axis.text = element_text(size=12, color = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('../figures/20220112/jitter_PTPN11_log2.pdf', width = 6, height = 6, useDingbats = FALSE)



# pca
pca_res <- prcomp(data[2:49], scale. = TRUE)
pca_res_rot <- data.table(pca_res$rotation, keep.rownames=TRUE)
pca_res_rot[, Sample := sapply(rn, function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_"))]
pca_res_rot[, Sample := factor(Sample, levels = c("EGFP", "E2D1_aCS3", "E2D1_aNSa5", "E2D1_aCS3_C85A", "E2D1_aCS3_V33R", "E2D1_aCS3_C85A_V33R", "E2D1_aCS3_F62A", "E2B_aCS3", "VHL_aCS3", "VHL_aCS3_V33R", "siRNA_control", "siRNA_shp2"))]

gg <- ggplot(pca_res_rot, aes(x = PC1, y = PC2, col = Sample)) +
geom_point(size = 1) +
theme_bw()

ggsave('../figures/20220112/pca.pdf', width = 8, height = 6, useDingbats = FALSE)

round(as.numeric(100*pca_res$sdev/sum(pca_res$sdev))[1],2) # 49% PC1
round(as.numeric(100*pca_res$sdev/sum(pca_res$sdev))[2],2) # 9% PC2
```



### venn diagrams and intersecting lists

```r
#cd ~/tables/20220112
#R

library(data.table)
library(VennDiagram)

# function for NA padding when combining tables with different number of rows
cbindPad <- function(...){
args <- list(...)
n <- sapply(args,nrow)
mx <- max(n)
pad <- function(x, mx){
    if (nrow(x) < mx){
        nms <- colnames(x)
        padTemp <- matrix(NA, mx - nrow(x), ncol(x))
        colnames(padTemp) <- nms
        if (ncol(x)==0) {
          return(padTemp)
        } else {
        return(rbind(x,padTemp))
          }
    }
    else{
        return(x)
    }
}
rs <- lapply(args,pad,mx)
return(do.call(cbind,rs))
}


# Enlarge the view width when printing tables
options(width = 300)


comparisons <- list(
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_C85A_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_F62A_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "E2D1_aCS3_C85A_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_C85A_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aNSa5_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2B_aCS3_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("E2B_aCS3_vs_EGFP", "E2D1_aNSa5_vs_EGFP"),
  c("E2D1_aNSa5_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("E2B_aCS3_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("VHL_aCS3_vs_EGFP", "VHL_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "VHL_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("E2B_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("E2D1_aNSa5_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("VHL_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control")
)


########
# down #
########

for (p in comparisons) {
  print(p)
  p1 <- fread(sprintf("%s.down.txt", p[1]))
  p2 <- fread(sprintf("%s.down.txt", p[2]))
  # venn
  venn.plot <- draw.pairwise.venn(
    area1 = length(p1$rn),
    area2 = length(p2$rn),
    cross.area = length(intersect(p1$rn, p2$rn)),
    category = c(sprintf("%s\n (%s)", p[1], length(p1$rn)), sprintf("%s\n (%s)", p[2], length(p2$rn))),
    euler.d = FALSE,
    scaled = FALSE,
    inverted = FALSE,
    fill = c("white", "white"),
    cex = 1.5,
    fontfamily = "sans",
    cat.pos = c(45, -45),
    cat.dist = c(0.15, 0.15),
    cat.cex = 1.5,
    cat.fontfamily = "sans",
    print.mode = "raw",
    margin = 0.25)
  pdf(sprintf("../../figures/20220112/venn_%s_cf_%s.down.pdf", p[1], p[2]), width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
  g <- grid.draw(venn.plot)
  dev.off()
  i <- p1[rn %in% intersect(p1$rn,p2$rn)][,.(rn, logFC, adj.P.Val)]
  l <- p1[rn %in% setdiff(p1$rn,p2$rn)][,.(rn, logFC, adj.P.Val)]
  r <- p2[rn %in% setdiff(p2$rn,p1$rn)][,.(rn, logFC, adj.P.Val)]
  lir <- data.frame(cbindPad(l,i,r))
  lir <- setNames(lir, c(sprintf("gene.%s", p[1]), sprintf("logFC.%s", p[1]), sprintf("FDR.%s", p[1]), sprintf("gene.%s", "common"), sprintf("logFC.%s", "common"), sprintf("FDR.%s", "common"), sprintf("gene.%s", p[2]), sprintf("logFC.%s", p[2]), sprintf("FDR.%s", p[2])))
  write.table(lir, sprintf("%s_cf_%s.down.txt", p[1], p[2]), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
}
```



## pxanalytics and limma analysis with pathways

```{r}
#cd ~/data
#R

library(devtools)
load_all("~/R/3.6.0/pxanalytics")
library(tidyverse)
#library(magrittr)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 300)


# Load and select data
data <- fread("20211207_134123_PE203_2021_bioPRO_normalized_proteins.txt")
data <- data[,c(2,51:98)]
colnames(data) <- c("genes", "EGFP_1", "EGFP_2", "EGFP_3", "EGFP_4", "E2D1_aCS3_1", "E2D1_aCS3_2", "E2D1_aCS3_3", "E2D1_aCS3_4", "E2D1_aNSa5_1", "E2D1_aNSa5_2", "E2D1_aNSa5_3", "E2D1_aNSa5_4", "E2D1_aCS3_V33R_1", "E2D1_aCS3_V33R_2", "E2D1_aCS3_V33R_3", "E2D1_aCS3_V33R_4", "E2D1_aCS3_C85A_1", "E2D1_aCS3_C85A_2", "E2D1_aCS3_C85A_3", "E2D1_aCS3_C85A_4", "E2D1_aCS3_C85A_V33R_1", "E2D1_aCS3_C85A_V33R_2", "E2D1_aCS3_C85A_V33R_3", "E2D1_aCS3_C85A_V33R_4", "VHL_aCS3_V33R_1", "VHL_aCS3_V33R_2", "VHL_aCS3_V33R_3", "VHL_aCS3_V33R_4", "siRNA_control_1", "siRNA_control_2", "siRNA_control_3", "siRNA_control_4", "siRNA_shp2_1", "siRNA_shp2_2", "siRNA_shp2_3", "siRNA_shp2_4", "E2D1_aCS3_F62A_1", "E2D1_aCS3_F62A_2", "E2D1_aCS3_F62A_3", "E2D1_aCS3_F62A_4", "E2B_aCS3_1", "E2B_aCS3_2", "E2B_aCS3_3", "E2B_aCS3_4", "VHL_aCS3_1", "VHL_aCS3_2", "VHL_aCS3_3", "VHL_aCS3_4")
data <- unique(data, by = "genes") # remove duplicated rows # distinct function
nrow(data) # 7657
data <- data.frame(data)
rownames(data) <- data$genes

# Assay data
## Assay data is a data frame of protein intensities with sample names given as
## column names and protein identifiers reported as row names.
assay <- na.omit(data[,2:49])
nrow(data) # 7657
nrow(assay) # 7657
length(rownames(assay)) # 7657
length(unique(rownames(assay))) # 7657


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_")), "Replicate" = sapply(colnames(assay), function (x) paste(tail(unlist(strsplit(x, "_")), 1), collapse="_")), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
px_norm <- px_raw %>% pxnormalise(method = "total")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "E2D1_aCS3_vs_EGFP" = E2D1_aCS3-EGFP,
  "E2D1_aCS3_V33R_vs_EGFP" = E2D1_aCS3_V33R-EGFP,
  "E2D1_aCS3_C85A_vs_EGFP" = E2D1_aCS3_C85A-EGFP,
  "E2D1_aCS3_C85A_V33R_vs_EGFP" = E2D1_aCS3_C85A_V33R-EGFP,
  "E2D1_aCS3_F62A_vs_EGFP" = E2D1_aCS3_F62A-EGFP,
  "E2D1_aNSa5_vs_EGFP" = E2D1_aNSa5-EGFP,
  "E2B_aCS3_vs_EGFP" = E2B_aCS3-EGFP,
  "VHL_aCS3_vs_EGFP" = VHL_aCS3-EGFP,
  "VHL_aCS3_V33R_vs_EGFP" = VHL_aCS3_V33R-EGFP,
  "EGFP_vs_siRNA_control" = EGFP-siRNA_control,
  "E2D1_aCS3_vs_siRNA_control" = E2D1_aCS3-siRNA_control,
  "E2D1_aNSa5_vs_siRNA_control" = E2D1_aNSa5-siRNA_control,
  "E2B_aCS3_vs_siRNA_control" = E2B_aCS3-siRNA_control,
  "VHL_aCS3_vs_siRNA_control" = VHL_aCS3-siRNA_control,
  "siRNA_shp2_vs_siRNA_control" = siRNA_shp2-siRNA_control,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# pathway data
pathway <- fread("pathway_lists.csv", header = TRUE)
pathway[, c("V2","V4","V6","V8"):=NULL]

pathway[Ub_pathway!=""]$Ub_pathway


# stats, PTPN11, tables, volcanos

for (c in colnames(contrast.matrix)){
  detable <- data.table(topTable(fit2, coef=c, adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)
  cat(c, "\t", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05]), "\t", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05]), "\n")
  print(detable[rn == "PTPN11"])
  detable[, Ub_pathway := ifelse(rn %in% pathway[Ub_pathway!=""]$Ub_pathway, "y", "n")]
  detable[, Cell_cycle_pathway := ifelse(rn %in% pathway[Cell_cycle_pathway!=""]$Cell_cycle_pathway, "y", "n")]
  detable[, E2D1_interactors := ifelse(rn %in% pathway[E2D1_interactors!=""]$E2D1_interactors, "y", "n")]
  detable[, MAPK_pathway := ifelse(rn %in% pathway[MAPK_pathway!=""]$MAPK_pathway, "y", "n")]
  detable[, p53_pathway := ifelse(rn %in% pathway[p53_pathway!=""]$p53_pathway, "y", "n")]
  # write table
  write.table(detable[,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220131/", c, ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  cat("****************************************\n")
  # volcano Ub_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220131/volcano_', c, '.Ub_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano Cell_cycle_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220131/volcano_', c, '.Cell_cycle_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano E2D1_interactors
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220131/volcano_', c, '.E2D1_interactors.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano MAPK_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220131/volcano_', c, '.MAPK_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano p53_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-0.5, linetype="longdash", color = "black") +
  geom_vline(xintercept=0.5, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -0.5 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 0.5 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -0.5 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 0.5 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220131/volcano_', c, '.p53_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
}

E2D1_aCS3_vs_EGFP 	 93 	 143
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.165073 17.13661 -8.166509 6.350491e-10 6.096292e-07 12.47788
****************************************
E2D1_aCS3_V33R_vs_EGFP 	 124 	 407
       rn     logFC  AveExpr         t     P.Value  adj.P.Val         B
1: PTPN11 -0.433324 17.13661 -3.037359 0.004274164 0.04357826 -2.335119
****************************************
E2D1_aCS3_C85A_vs_EGFP 	 144 	 241
       rn       logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 -0.08704418 17.13661 -0.610131 0.545367 0.7893251 -6.395549
****************************************
E2D1_aCS3_C85A_V33R_vs_EGFP 	 243 	 284
       rn        logFC  AveExpr           t   P.Value adj.P.Val         B
1: PTPN11 -0.001972656 17.13661 -0.01382722 0.9890394 0.9958144 -6.590937
****************************************
E2D1_aCS3_F62A_vs_EGFP 	 228 	 88
       rn      logFC  AveExpr          t  P.Value adj.P.Val         B
1: PTPN11 -0.1318922 17.13661 -0.9244906 0.361001 0.6506443 -5.781113
****************************************
E2D1_aNSa5_vs_EGFP 	 110 	 57
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.8227283 17.13661 -5.766865 1.146419e-06 0.0005486334 5.425834
****************************************
E2B_aCS3_vs_EGFP 	 242 	 147
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.890644 17.13661 -13.25236 6.565184e-16 2.513481e-12 24.69009
****************************************
VHL_aCS3_vs_EGFP 	 372 	 45
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.947983 17.13661 -13.65427 2.534462e-16 9.703189e-13 25.46349
****************************************
VHL_aCS3_V33R_vs_EGFP 	 48 	 45
       rn       logFC  AveExpr          t   P.Value adj.P.Val        B
1: PTPN11 -0.04845112 17.13661 -0.3396153 0.7359958 0.9208365 -6.15588
****************************************
EGFP_vs_siRNA_control 	 265 	 221
       rn     logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 0.0641864 17.13661 0.4499107 0.655303 0.8252947 -6.213415
****************************************
E2D1_aCS3_vs_siRNA_control 	 654 	 708
       rn     logFC  AveExpr         t     P.Value    adj.P.Val        B
1: PTPN11 -1.100886 17.13661 -7.716598 2.50521e-09 5.480684e-07 11.30448
****************************************
E2D1_aNSa5_vs_siRNA_control 	 988 	 522
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.7585419 17.13661 -5.316954 4.771653e-06 0.0001574851 3.965503
****************************************
E2B_aCS3_vs_siRNA_control 	 87 	 119
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.826458 17.13661 -12.80245 1.948807e-15 7.461008e-12 24.14485
****************************************
VHL_aCS3_vs_siRNA_control 	 174 	 65
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.883797 17.13661 -13.20436 7.364813e-16 1.879746e-12 25.15552
****************************************
siRNA_shp2_vs_siRNA_control 	 36 	 34
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.426508 17.13661 -9.999022 3.027922e-12 7.843096e-09 16.91789
****************************************
```



## pxanalytics and limma analysis logFC 1/-1 w/ and w/o pathways

```{r}
#cd ~/data
#R

#library(devtools)
#load_all("~/R/3.6.0/pxanalytics")
library(pxanalytics)
library(tidyverse)
#library(magrittr)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 300)


# Load and select data
data <- fread("20211207_134123_PE203_2021_bioPRO_normalized_proteins.txt")
data <- data[,c(2,51:98)]
colnames(data) <- c("genes", "EGFP_1", "EGFP_2", "EGFP_3", "EGFP_4", "E2D1_aCS3_1", "E2D1_aCS3_2", "E2D1_aCS3_3", "E2D1_aCS3_4", "E2D1_aNSa5_1", "E2D1_aNSa5_2", "E2D1_aNSa5_3", "E2D1_aNSa5_4", "E2D1_aCS3_V33R_1", "E2D1_aCS3_V33R_2", "E2D1_aCS3_V33R_3", "E2D1_aCS3_V33R_4", "E2D1_aCS3_C85A_1", "E2D1_aCS3_C85A_2", "E2D1_aCS3_C85A_3", "E2D1_aCS3_C85A_4", "E2D1_aCS3_C85A_V33R_1", "E2D1_aCS3_C85A_V33R_2", "E2D1_aCS3_C85A_V33R_3", "E2D1_aCS3_C85A_V33R_4", "VHL_aCS3_V33R_1", "VHL_aCS3_V33R_2", "VHL_aCS3_V33R_3", "VHL_aCS3_V33R_4", "siRNA_control_1", "siRNA_control_2", "siRNA_control_3", "siRNA_control_4", "siRNA_shp2_1", "siRNA_shp2_2", "siRNA_shp2_3", "siRNA_shp2_4", "E2D1_aCS3_F62A_1", "E2D1_aCS3_F62A_2", "E2D1_aCS3_F62A_3", "E2D1_aCS3_F62A_4", "E2B_aCS3_1", "E2B_aCS3_2", "E2B_aCS3_3", "E2B_aCS3_4", "VHL_aCS3_1", "VHL_aCS3_2", "VHL_aCS3_3", "VHL_aCS3_4")
data <- unique(data, by = "genes") # remove duplicated rows # distinct function
nrow(data) # 7657
data <- data.frame(data)
rownames(data) <- data$genes

# Assay data
assay <- na.omit(data[,2:49])
nrow(data) # 7657
nrow(assay) # 7657
length(rownames(assay)) # 7657
length(unique(rownames(assay))) # 7657


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_")), "Replicate" = sapply(colnames(assay), function (x) paste(tail(unlist(strsplit(x, "_")), 1), collapse="_")), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
px_norm <- px_raw %>% pxnormalise(method = "total")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "E2D1_aCS3_vs_EGFP" = E2D1_aCS3-EGFP,
  "E2D1_aCS3_V33R_vs_EGFP" = E2D1_aCS3_V33R-EGFP,
  "E2D1_aCS3_C85A_vs_EGFP" = E2D1_aCS3_C85A-EGFP,
  "E2D1_aCS3_C85A_V33R_vs_EGFP" = E2D1_aCS3_C85A_V33R-EGFP,
  "E2D1_aCS3_F62A_vs_EGFP" = E2D1_aCS3_F62A-EGFP,
  "E2D1_aNSa5_vs_EGFP" = E2D1_aNSa5-EGFP,
  "E2B_aCS3_vs_EGFP" = E2B_aCS3-EGFP,
  "VHL_aCS3_vs_EGFP" = VHL_aCS3-EGFP,
  "VHL_aCS3_V33R_vs_EGFP" = VHL_aCS3_V33R-EGFP,
  "VHL_aCS3_vs_VHL_aCS3_V33R" = VHL_aCS3-VHL_aCS3_V33R,
  "EGFP_vs_siRNA_control" = EGFP-siRNA_control,
  "E2D1_aCS3_vs_siRNA_control" = E2D1_aCS3-siRNA_control,
  "E2D1_aNSa5_vs_siRNA_control" = E2D1_aNSa5-siRNA_control,
  "E2B_aCS3_vs_siRNA_control" = E2B_aCS3-siRNA_control,
  "VHL_aCS3_vs_siRNA_control" = VHL_aCS3-siRNA_control,
  "siRNA_shp2_vs_siRNA_control" = siRNA_shp2-siRNA_control,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# pathway data
pathway <- fread("pathway_lists.csv", header = TRUE)
pathway[, c("V2","V4","V6","V8"):=NULL]

pathway[Ub_pathway!=""]$Ub_pathway


# stats, PTPN11, tables, volcanos

for (c in colnames(contrast.matrix)){
  # load data
  detable <- data.table(topTable(fit2, coef=c, adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)
  # print data
  cat(c, "\t", nrow(detable[logFC < -1 & adj.P.Val < 0.05]), "\t", nrow(detable[logFC > 1 & adj.P.Val < 0.05]), "\n")
  print(detable[rn == "PTPN11"])
  # manage pathways
  detable[, Ub_pathway := ifelse(rn %in% pathway[Ub_pathway!=""]$Ub_pathway, "y", "n")]
  detable[, Cell_cycle_pathway := ifelse(rn %in% pathway[Cell_cycle_pathway!=""]$Cell_cycle_pathway, "y", "n")]
  detable[, E2D1_interactors := ifelse(rn %in% pathway[E2D1_interactors!=""]$E2D1_interactors, "y", "n")]
  detable[, MAPK_pathway := ifelse(rn %in% pathway[MAPK_pathway!=""]$MAPK_pathway, "y", "n")]
  detable[, p53_pathway := ifelse(rn %in% pathway[p53_pathway!=""]$p53_pathway, "y", "n")]
  # write tables
  ## down
  write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220215/", c, ".down.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  ## up
  write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220215/", c, ".up.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  ## all
  write.table(detable[,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220215/", c, ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  cat("****************************************\n")
  # volcano
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
  geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)
  ggsave(paste('../figures/20220215/volcano_', c, '.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano Ub_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -1 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 1 & adj.P.Val < 0.05 & Ub_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220215/volcano_', c, '.Ub_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano Cell_cycle_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -1 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 1 & adj.P.Val < 0.05 & Cell_cycle_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220215/volcano_', c, '.Cell_cycle_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano E2D1_interactors
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -1 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 1 & adj.P.Val < 0.05 & E2D1_interactors=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220215/volcano_', c, '.E2D1_interactors.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano MAPK_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -1 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 1 & adj.P.Val < 0.05 & MAPK_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220215/volcano_', c, '.MAPK_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
  # volcano p53_pathway
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-1, linetype="longdash", color = "black") +
  geom_vline(xintercept=1, linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = detable[logFC < -1 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1) +
  geom_text_repel(data = detable[logFC > 1 & adj.P.Val < 0.05 & p53_pathway=="y"], aes(label = rn), force = 10, size = 3, box.padding = 1)
  ggsave(paste('../figures/20220215/volcano_', c, '.p53_pathway.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
}

E2D1_aCS3_vs_EGFP 	 34 	 28
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.165073 17.13661 -8.166509 6.350491e-10 6.096292e-07 12.47788
****************************************
E2D1_aCS3_V33R_vs_EGFP 	 45 	 68
       rn     logFC  AveExpr         t     P.Value  adj.P.Val         B
1: PTPN11 -0.433324 17.13661 -3.037359 0.004274164 0.04357826 -2.335119
****************************************
E2D1_aCS3_C85A_vs_EGFP 	 64 	 92
       rn       logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 -0.08704418 17.13661 -0.610131 0.545367 0.7893251 -6.395549
****************************************
E2D1_aCS3_C85A_V33R_vs_EGFP 	 129 	 94
       rn        logFC  AveExpr           t   P.Value adj.P.Val         B
1: PTPN11 -0.001972656 17.13661 -0.01382722 0.9890394 0.9958144 -6.590937
****************************************
E2D1_aCS3_F62A_vs_EGFP 	 140 	 19
       rn      logFC  AveExpr          t  P.Value adj.P.Val         B
1: PTPN11 -0.1318922 17.13661 -0.9244906 0.361001 0.6506443 -5.781113
****************************************
E2D1_aNSa5_vs_EGFP 	 51 	 16
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.8227283 17.13661 -5.766865 1.146419e-06 0.0005486334 5.425834
****************************************
E2B_aCS3_vs_EGFP 	 172 	 29
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.890644 17.13661 -13.25236 6.565184e-16 2.513481e-12 24.69009
****************************************
VHL_aCS3_vs_EGFP 	 251 	 16
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.947983 17.13661 -13.65427 2.534462e-16 9.703189e-13 25.46349
****************************************
VHL_aCS3_V33R_vs_EGFP 	 32 	 18
       rn       logFC  AveExpr          t   P.Value adj.P.Val        B
1: PTPN11 -0.04845112 17.13661 -0.3396153 0.7359958 0.9208365 -6.15588
****************************************
VHL_aCS3_vs_VHL_aCS3_V33R 	 4 	 0
       rn     logFC  AveExpr         t      P.Value    adj.P.Val       B
1: PTPN11 -1.899532 17.13661 -13.31466 5.657717e-16 4.332114e-12 20.6478
****************************************
EGFP_vs_siRNA_control 	 64 	 109
       rn     logFC  AveExpr         t  P.Value adj.P.Val         B
1: PTPN11 0.0641864 17.13661 0.4499107 0.655303 0.8252947 -6.213415
****************************************
E2D1_aCS3_vs_siRNA_control 	 192 	 261
       rn     logFC  AveExpr         t     P.Value    adj.P.Val        B
1: PTPN11 -1.100886 17.13661 -7.716598 2.50521e-09 5.480684e-07 11.30448
****************************************
E2D1_aNSa5_vs_siRNA_control 	 269 	 179
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.7585419 17.13661 -5.316954 4.771653e-06 0.0001574851 3.965503
****************************************
E2B_aCS3_vs_siRNA_control 	 35 	 23
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.826458 17.13661 -12.80245 1.948807e-15 7.461008e-12 24.14485
****************************************
VHL_aCS3_vs_siRNA_control 	 73 	 17
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.883797 17.13661 -13.20436 7.364813e-16 1.879746e-12 25.15552
****************************************
siRNA_shp2_vs_siRNA_control 	 12 	 6
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.426508 17.13661 -9.999022 3.027922e-12 7.843096e-09 16.91789
****************************************
```



### intersecting tables and logFC scatterplots

```r
cd ~/tables/20220215
R

library(data.table)
library(ggplot2)

# function for NA padding when combining tables with different number of rows
cbindPad <- function(...){
args <- list(...)
n <- sapply(args,nrow)
mx <- max(n)
pad <- function(x, mx){
    if (nrow(x) < mx){
        nms <- colnames(x)
        padTemp <- matrix(NA, mx - nrow(x), ncol(x))
        colnames(padTemp) <- nms
        if (ncol(x)==0) {
          return(padTemp)
        } else {
        return(rbind(x,padTemp))
          }
    }
    else{
        return(x)
    }
}
rs <- lapply(args,pad,mx)
return(do.call(cbind,rs))
}


# Enlarge the view width when printing tables
options(width = 300)


comparisons <- list(
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_C85A_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aCS3_F62A_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "E2D1_aCS3_C85A_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_C85A_vs_EGFP", "E2D1_aCS3_C85A_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2D1_aNSa5_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "E2B_aCS3_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("E2B_aCS3_vs_EGFP", "E2D1_aNSa5_vs_EGFP"),
  c("E2D1_aNSa5_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("E2B_aCS3_vs_EGFP", "VHL_aCS3_vs_EGFP"),
  c("VHL_aCS3_vs_EGFP", "VHL_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_V33R_vs_EGFP", "VHL_aCS3_V33R_vs_EGFP"),
  c("E2D1_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("E2B_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("E2D1_aNSa5_vs_EGFP", "siRNA_shp2_vs_siRNA_control"),
  c("VHL_aCS3_vs_EGFP", "siRNA_shp2_vs_siRNA_control")
)

###############
# down and up #
###############

for (p in comparisons) {
  print(p)
  # load data an filter
  p1 <- fread(sprintf("%s.txt", p[1]))
  p1_filter <- p1[adj.P.Val<0.05 & ((logFC < -1) | (logFC > 1))]
  p2 <- fread(sprintf("%s.txt", p[2]))
  p2_filter <- p2[adj.P.Val<0.05 & ((logFC < -1) | (logFC > 1))]
  # intersecting tables
  i1 <- p1_filter[rn %in% intersect(p1_filter$rn,p2_filter$rn)][,.(rn, logFC, adj.P.Val)]
  i2 <- p2_filter[rn %in% intersect(p1_filter$rn,p2_filter$rn)][,.(rn, logFC, adj.P.Val)]
  i <- merge(i1, i2, by = "rn", suffixes = c(".1", ".2"))[order(adj.P.Val.1)]
  l <- p1_filter[rn %in% setdiff(p1_filter$rn,p2_filter$rn)][,.(rn, logFC, adj.P.Val)]
  r <- p2_filter[rn %in% setdiff(p2_filter$rn,p1_filter$rn)][,.(rn, logFC, adj.P.Val)]
  lir <- data.frame(cbindPad(l,i,r))
  lir <- setNames(lir, c(sprintf("gene.%s.only", p[1]), sprintf("logFC.%s.only", p[1]), sprintf("FDR.%s.only", p[1]), sprintf("gene.%s", "common"), sprintf("logFC.%s.common", p[1]), sprintf("FDR.%s.common", p[1]), sprintf("logFC.%s.common", p[2]), sprintf("FDR.%s.common", p[2]), sprintf("gene.%s.only", p[2]), sprintf("logFC.%s.only", p[2]), sprintf("FDR.%s.only", p[2])))
  write.table(lir, sprintf("%s_cf_%s.txt", p[1], p[2]), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  # logFC scatterplots
  p1_p2 <- merge(p1, p2, by = "rn")[rn != ""]
  gg <- ggplot(p1_p2, aes(x = logFC.x, y = logFC.y)) +
  geom_point(size = 1, col = "black") +
  geom_vline(xintercept=0, linetype="longdash", color = "maroon") +
  geom_hline(yintercept=0, linetype="longdash", color = "maroon") +
  xlab(p[1]) +
  ylab(p[2]) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 7), clip = "off") +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"))
  ggsave(sprintf('../../figures/20220215/scatterplot_%s_cf_%s.pdf', p[1], p[2]), width = 6, height = 6, useDingbats = FALSE)
}
```



## pxanalytics and limma analysis logFC log2(1.5)/-log2(1.5) w/ and w/o pathways

New set of comparisons requested by Jonathan in the chat

```{r}
#cd ~/data
#R

#library(devtools)
#load_all("~/R/3.6.0/pxanalytics")
library(pxanalytics)
library(tidyverse)
#library(magrittr)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 300)


# Load and select data
data <- fread("20211207_134123_PE203_2021_bioPRO_normalized_proteins.txt")
data <- data[,c(2,51:98)]
colnames(data) <- c("genes", "EGFP_1", "EGFP_2", "EGFP_3", "EGFP_4", "E2D1_aCS3_1", "E2D1_aCS3_2", "E2D1_aCS3_3", "E2D1_aCS3_4", "E2D1_aNSa5_1", "E2D1_aNSa5_2", "E2D1_aNSa5_3", "E2D1_aNSa5_4", "E2D1_aCS3_V33R_1", "E2D1_aCS3_V33R_2", "E2D1_aCS3_V33R_3", "E2D1_aCS3_V33R_4", "E2D1_aCS3_C85A_1", "E2D1_aCS3_C85A_2", "E2D1_aCS3_C85A_3", "E2D1_aCS3_C85A_4", "E2D1_aCS3_C85A_V33R_1", "E2D1_aCS3_C85A_V33R_2", "E2D1_aCS3_C85A_V33R_3", "E2D1_aCS3_C85A_V33R_4", "VHL_aCS3_V33R_1", "VHL_aCS3_V33R_2", "VHL_aCS3_V33R_3", "VHL_aCS3_V33R_4", "siRNA_control_1", "siRNA_control_2", "siRNA_control_3", "siRNA_control_4", "siRNA_shp2_1", "siRNA_shp2_2", "siRNA_shp2_3", "siRNA_shp2_4", "E2D1_aCS3_F62A_1", "E2D1_aCS3_F62A_2", "E2D1_aCS3_F62A_3", "E2D1_aCS3_F62A_4", "E2B_aCS3_1", "E2B_aCS3_2", "E2B_aCS3_3", "E2B_aCS3_4", "VHL_aCS3_1", "VHL_aCS3_2", "VHL_aCS3_3", "VHL_aCS3_4")
data <- unique(data, by = "genes") # remove duplicated rows # distinct function
nrow(data) # 7657
data <- data.frame(data)
rownames(data) <- data$genes

# Assay data
assay <- na.omit(data[,2:49])
nrow(data) # 7657
nrow(assay) # 7657
length(rownames(assay)) # 7657
length(unique(rownames(assay))) # 7657


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(head(unlist(strsplit(x, "_")),-1), collapse="_")), "Replicate" = sapply(colnames(assay), function (x) paste(tail(unlist(strsplit(x, "_")), 1), collapse="_")), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
px_norm <- px_raw %>% pxnormalise(method = "total")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "E2D1_aCS3_vs_E2D1_aCS3_C85A" = E2D1_aCS3-E2D1_aCS3_C85A,
  "E2D1_aCS3_vs_E2D1_aCS3_V33R" = E2D1_aCS3-E2D1_aCS3_V33R,
  "E2D1_aCS3_vs_E2D1_aCS3_C85A_V33R" = E2D1_aCS3-E2D1_aCS3_C85A_V33R,
  "E2D1_aCS3_vs_E2D1_aCS3_F62A" = E2D1_aCS3-E2D1_aCS3_F62A,
  "E2D1_aCS3_vs_E2D1_aNSa5" = E2D1_aCS3-E2D1_aNSa5,
  "E2D1_aCS3_vs_E2B_aCS3" = E2D1_aCS3-E2B_aCS3,
  "E2D1_aCS3_vs_VHL_aCS3" = E2D1_aCS3-VHL_aCS3,
  "E2B_aCS3_vs_E2D1_aNSa5" = E2B_aCS3-E2D1_aNSa5,
  "E2B_aCS3_vs_VHL_aCS3" = E2B_aCS3-VHL_aCS3,
  "E2B_aCS3_vs_E2D1_aCS3_C85A_V33R" = E2B_aCS3-E2D1_aCS3_C85A_V33R,
  "E2D1_aNSa5_vs_VHL_aCS3" = E2D1_aNSa5-VHL_aCS3,
  "E2D1_aNSa5_vs_E2D1_aCS3_C85A_V33R" = E2D1_aNSa5-E2D1_aCS3_C85A_V33R,
  "E2D1_aCS3_V33R_vs_VHL_aCS3_V33R" = E2D1_aCS3_V33R-VHL_aCS3_V33R,
  "E2D1_aCS3_V33R_vs_E2D1_aCS3_C85A_V33R" = E2D1_aCS3_V33R-E2D1_aCS3_C85A_V33R,
  "E2D1_aCS3_vs_siRNA_control" = E2D1_aCS3-siRNA_control,
  "E2D1_aCS3_C85A_vs_siRNA_control" = E2D1_aCS3_C85A-siRNA_control,
  "E2D1_aCS3_V33R_vs_siRNA_control" = E2D1_aCS3_V33R-siRNA_control,
  "E2D1_aCS3_C85A_V33R_vs_siRNA_control" = E2D1_aCS3_C85A_V33R-siRNA_control,
  "E2D1_aCS3_F62A_vs_siRNA_control" = E2D1_aCS3_F62A-siRNA_control,
  "E2D1_aNSa5_vs_siRNA_control" = E2D1_aNSa5-siRNA_control,
  "E2B_aCS3_vs_siRNA_control" = E2B_aCS3-siRNA_control,
  "VHL_aCS3_vs_siRNA_control" = VHL_aCS3-siRNA_control,
  "VHL_aCS3_V33R_vs_siRNA_control" = VHL_aCS3_V33R-siRNA_control,
  "VHL_aCS3_V33R_vs_siRNA_control" = VHL_aCS3_V33R-siRNA_control,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# pathway data
pathway <- fread("pathway_lists.csv", header = TRUE)
pathway[, c("V2","V4","V6","V8"):=NULL]

pathway[Ub_pathway!=""]$Ub_pathway


# stats, PTPN11, tables, volcanos

for (c in colnames(contrast.matrix)){
  # load data
  detable <- data.table(topTable(fit2, coef=c, adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)
  # print data
  cat(c, "\t", nrow(detable[logFC < -log2(1.5) & adj.P.Val < 0.05]), "\t", nrow(detable[logFC > log2(1.5) & adj.P.Val < 0.05]), "\n")
  print(detable[rn == "PTPN11"])
  # manage pathways
  detable[, Ub_pathway := ifelse(rn %in% pathway[Ub_pathway!=""]$Ub_pathway, "y", "n")]
  detable[, Cell_cycle_pathway := ifelse(rn %in% pathway[Cell_cycle_pathway!=""]$Cell_cycle_pathway, "y", "n")]
  detable[, E2D1_interactors := ifelse(rn %in% pathway[E2D1_interactors!=""]$E2D1_interactors, "y", "n")]
  detable[, MAPK_pathway := ifelse(rn %in% pathway[MAPK_pathway!=""]$MAPK_pathway, "y", "n")]
  detable[, p53_pathway := ifelse(rn %in% pathway[p53_pathway!=""]$p53_pathway, "y", "n")]
  # write tables
  ## down
  write.table(detable[logFC < -log2(1.5) & adj.P.Val < 0.05][,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220606/", c, ".down.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  ## up
  write.table(detable[logFC > log2(1.5) & adj.P.Val < 0.05][,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220606/", c, ".up.txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  ## all
  write.table(detable[,.(rn, logFC, adj.P.Val, Ub_pathway, Cell_cycle_pathway, E2D1_interactors, MAPK_pathway, p53_pathway)], paste("../tables/20220606/", c, ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
  cat("****************************************\n")
  # volcano
  gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 1, col = "lightgray") +
  geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
  geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
  geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
  geom_point(data = detable[logFC < -log2(1.5) & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[logFC > log2(1.5) & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
  geom_point(data = detable[rn=="PTPN11"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
  ylab(expression("-log"[10]*"FDR")) +
  xlab(expression("log"[2]*"FC")) +
  ggtitle(paste(unlist(strsplit(c, "_vs_"))[1], unlist(strsplit(c, "_vs_"))[2], sep=" vs ")) +
  theme_bw() +
  coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20), clip = "off") +
  annotate("text", x = -5, y = 20, label = sprintf("n = %s", nrow(detable[logFC < -log2(1.5) & adj.P.Val < 0.05])), size = 5) +
  annotate("text", x = 5, y = 20, label = sprintf("n = %s", nrow(detable[logFC > log2(1.5) & adj.P.Val < 0.05])), size = 5) +
  theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
  geom_text_repel(data = head(detable[logFC < -log2(1.5) & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
  geom_text_repel(data = head(detable[logFC > log2(1.5) & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)
  ggsave(paste('../figures/20220606/volcano_', c, '.pdf', sep=""), width = 6, height = 6, useDingbats = FALSE)
}

E2D1_aCS3_vs_E2D1_aCS3_C85A 	 202 	 164
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.078029 17.13661 -7.556378 4.104252e-09 2.856932e-06 10.75968
****************************************
E2D1_aCS3_vs_E2D1_aCS3_V33R 	 3 	 0
       rn      logFC  AveExpr        t      P.Value  adj.P.Val        B
1: PTPN11 -0.7317489 17.13661 -5.12915 8.629205e-06 0.03303691 2.964064
****************************************
E2D1_aCS3_vs_E2D1_aCS3_C85A_V33R 	 148 	 133
       rn   logFC  AveExpr         t      P.Value    adj.P.Val       B
1: PTPN11 -1.1631 17.13661 -8.152682 6.621965e-10 4.225366e-07 12.5032
****************************************
E2D1_aCS3_vs_E2D1_aCS3_F62A 	 164 	 389
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.033181 17.13661 -7.242018 1.088465e-08 1.666876e-05 9.781887
****************************************
E2D1_aCS3_vs_E2D1_aNSa5 	 1 	 3
       rn      logFC  AveExpr         t    P.Value adj.P.Val         B
1: PTPN11 -0.3423446 17.13661 -2.399644 0.02136177 0.3755698 -3.312481
****************************************
E2D1_aCS3_vs_E2B_aCS3 	 306 	 419
       rn     logFC  AveExpr        t      P.Value    adj.P.Val       B
1: PTPN11 0.7255716 17.13661 5.085851 9.888689e-06 0.0006595082 3.36899
****************************************
E2D1_aCS3_vs_VHL_aCS3 	 204 	 733
       rn     logFC  AveExpr        t     P.Value    adj.P.Val        B
1: PTPN11 0.7829101 17.13661 5.487762 2.77909e-06 0.0001911148 4.566244
****************************************
E2B_aCS3_vs_E2D1_aNSa5 	 325 	 578
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.067916 17.13661 -7.485494 5.109934e-09 2.445423e-06 10.59813
****************************************
E2B_aCS3_vs_VHL_aCS3 	 12 	 17
       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
1: PTPN11 0.05733852 17.13661 0.4019109 0.6899766  0.916232 -6.066214
****************************************
E2B_aCS3_vs_E2D1_aCS3_C85A_V33R 	 270 	 293
       rn     logFC  AveExpr         t      P.Value   adj.P.Val        B
1: PTPN11 -1.888672 17.13661 -13.23853 6.785997e-16 1.29901e-12 25.66144
****************************************
E2D1_aNSa5_vs_VHL_aCS3 	 336 	 524
       rn    logFC  AveExpr        t      P.Value   adj.P.Val        B
1: PTPN11 1.125255 17.13661 7.887405 1.484187e-09 8.74186e-07 11.77387
****************************************
E2D1_aNSa5_vs_E2D1_aCS3_C85A_V33R 	 428 	 184
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.8207557 17.13661 -5.753038 1.197876e-06 0.0001069218 5.343792
****************************************
E2D1_aCS3_V33R_vs_VHL_aCS3_V33R 	 143 	 299
       rn      logFC  AveExpr         t    P.Value  adj.P.Val         B
1: PTPN11 -0.3848729 17.13661 -2.697743 0.01031341 0.08041729 -3.182731
****************************************
E2D1_aCS3_V33R_vs_E2D1_aCS3_C85A_V33R 	 101 	 199
       rn      logFC  AveExpr         t     P.Value  adj.P.Val         B
1: PTPN11 -0.4313513 17.13661 -3.023531 0.004434564 0.06107096 -2.357558
****************************************
E2D1_aCS3_vs_siRNA_control 	 535 	 600
       rn     logFC  AveExpr         t     P.Value    adj.P.Val        B
1: PTPN11 -1.100886 17.13661 -7.716598 2.50521e-09 5.480684e-07 11.30448
****************************************
E2D1_aCS3_C85A_vs_siRNA_control 	 321 	 357
       rn       logFC  AveExpr          t   P.Value adj.P.Val         B
1: PTPN11 -0.02285779 17.13661 -0.1602203 0.8735469 0.9440356 -6.716446
****************************************
E2D1_aCS3_V33R_vs_siRNA_control 	 430 	 716
       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
1: PTPN11 -0.3691376 17.13661 -2.587448 0.0135727 0.0539036 -3.624738
****************************************
E2D1_aCS3_C85A_V33R_vs_siRNA_control 	 456 	 363
       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
1: PTPN11 0.06221374 17.13661 0.4360835 0.6652166 0.8166689 -6.640516
****************************************
E2D1_aCS3_F62A_vs_siRNA_control 	 135 	 83
       rn       logFC  AveExpr          t   P.Value adj.P.Val         B
1: PTPN11 -0.06770582 17.13661 -0.4745799 0.6377728 0.8773142 -6.308144
****************************************
E2D1_aNSa5_vs_siRNA_control 	 815 	 427
       rn      logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -0.7585419 17.13661 -5.316954 4.771653e-06 0.0001574851 3.965503
****************************************
E2B_aCS3_vs_siRNA_control 	 77 	 85
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.826458 17.13661 -12.80245 1.948807e-15 7.461008e-12 24.14485
****************************************
VHL_aCS3_vs_siRNA_control 	 155 	 51
       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
1: PTPN11 -1.883797 17.13661 -13.20436 7.364813e-16 1.879746e-12 25.15552
****************************************
VHL_aCS3_V33R_vs_siRNA_control 	 95 	 107
       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
1: PTPN11 0.01573527 17.13661 0.1102955 0.9127491 0.9823852 -6.473276
****************************************
VHL_aCS3_V33R_vs_siRNA_control 	 95 	 107
       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
1: PTPN11 0.01573527 17.13661 0.1102955 0.9127491 0.9823852 -6.473276
****************************************
```


## Additional volcano plots

```{r}
cd ~/data/20230207
R

library(data.table)
library(ggplot2)
library(ggrepel)

# Enlarge the view width when printing tables
options(width = 300)

# Enable ggsave of png files
# https://stackoverflow.com/questions/24999983/r-unable-to-start-device-png-capabilities-has-true-for-png
options(bitmapType='cairo')


#############
# E2D1_aCS3 #
#############

# Load and select data
data <- fread("E2D1_aCS3.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5)) +
geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001][ID != c("UBE2D1")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_aCS3.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_aCS3.png', width = 4, height = 4, dpi = 300)


############
# E2B_aCS3 #
############

# Load and select data
data <- fread("E2B_aCS3.csv")
data <- data[ID!=""]


# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.000015][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-8, -1.5)) +
geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001][ID != c("UBE2B")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "UBE2B"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2B_aCS3.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2B_aCS3.png', width = 4, height = 4, dpi = 300)



############
# VHL_aCS3 #
############

# Load and select data
data <- fread("VHL_aCS3.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5)) +
geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001][ID != c("VHL")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "VHL"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkorchid3", color = "white", segment.color = "black")
ggsave('../../figures/20230207/VHL_aCS3.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/VHL_aCS3.png', width = 4, height = 4, dpi = 300)


#############
# E2D1_NSa5 #
#############

# Load and select data
data <- fread("E2D1_NSa5.csv")
data <- data[ID!=""]

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00004][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(-8, -1.5)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00004][ID != c("UBE2D1")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_NSa5.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_NSa5.png', width = 4, height = 4, dpi = 300)

# Volcano including SHP2
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00004][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(-8, -1.5), ylim = c(4, Inf)) +
geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-8, -5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00004][ID != c("UBE2D1")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_NSa5_2.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_NSa5_2.png', width = 4, height = 4, dpi = 300)


##################
# E2D1_aCS3_V33R #
##################

# Load and select data
data <- fread("E2D1_aCS3_V33R.csv")
data <- data[ID!=""]

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.000001][ID != c("UBE2D1")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(0, 8)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_aCS3_V33R.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_aCS3_V33R.png', width = 4, height = 4, dpi = 300)


#################
# VHL_aCS3_V33R #
#################

# Load and select data
data <- fread("VHL_aCS3_V33R.csv")
data <- data[ID!=""]

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001][ID != c("VHL")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1.5, Inf)) +
geom_label_repel(data = data[ID == "VHL"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkorchid3", color = "white", segment.color = "black")
ggsave('../../figures/20230207/VHL_aCS3_V33R.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/VHL_aCS3_V33R.png', width = 4, height = 4, dpi = 300)


##################
# E2D1_aCS3_C85A #
##################

# Load and select data
data <- fread("E2D1_aCS3_C85A.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, -1.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00000001][ID != c("UBE2D1")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0, xlim = c(-8, 8)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(10, 15), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_aCS3_C85A.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_aCS3_C85A.png', width = 4, height = 4, dpi = 300)


##############
# siRNA_SHP2 #
##############

# Load and select data
data <- fread("siRNA_SHP2.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0)) +
geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0), ylim = c(5, 7.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(0, Inf)) +
#geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/siRNA_SHP2.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/siRNA_SHP2.png', width = 4, height = 4, dpi = 300)


##################
# E2D1_aCS3_F62A #
##################

# Load and select data
data <- fread("E2D1_aCS3_F62A.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0), ylim = c(5, 7.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00001][ID != "UBE2D1"], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(1, Inf)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(1.5, 8), ylim = c(13, Inf), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_aCS3_F62A.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_aCS3_F62A.png', width = 4, height = 4, dpi = 300)


#######################
# E2D1_aCS3_C85A_V33R #
#######################

# Load and select data
data <- fread("E2D1_aCS3_C85A_V33R.csv")

# Volcano
gg <- ggplot(data, aes(x = log2FC, y = -log10(adj.p))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-log2(1.5), linetype="longdash", color = "black") +
geom_vline(xintercept=log2(1.5), linetype="longdash", color = "black") +
geom_point(data = data[log2FC < -log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkblue') +
geom_point(data = data[log2FC > log2(1.5) & adj.p < 0.05], aes(x = log2FC, y = -log10(adj.p)), color='darkred') +
ylab("False Discovery Rate") +
xlab(expression("Fold Change (log"[2]*"FC)")) +
theme_bw() +
coord_cartesian(xlim = c(-8, 8), ylim = c(0, 20), clip = "off") +
annotate("text", x = -2, y = 20, label = nrow(data[log2FC < -log2(1.5) & adj.p < 0.05]), size = 5, color='darkblue') +
annotate("text", x = 2, y = 20, label = nrow(data[log2FC > log2(1.5) & adj.p < 0.05]), size = 5, color='darkred') +
theme(axis.title = element_blank(), axis.text = element_text(size=14, color = "black")) +
geom_text_repel(data = data[log2FC < 0 & adj.p < 0.00001][ID != c("SHP2")], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0)) +
#geom_label_repel(data = data[ID == "SHP2"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(-Inf, 0), ylim = c(5, 7.5), fill = "darkblue", color = "white", segment.color = "black") +
geom_text_repel(data = data[log2FC > 0 & adj.p < 0.00000001][ID != "UBE2D1"], aes(label = ID), force = 1, size = 3.5, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, xlim = c(-Inf, 9)) +
geom_label_repel(data = data[ID == "UBE2D1"], aes(label = ID), force = 1, size = 4.5, max.overlaps = Inf, box.padding = 1, min.segment.length = 0, xlim = c(3, 9), ylim = c(15, 20), fill = "darkred", color = "white", segment.color = "black")
ggsave('../../figures/20230207/E2D1_aCS3_C85A_V33R.pdf', width = 4, height = 4, useDingbats = FALSE)
ggsave('../../figures/20230207/E2D1_aCS3_C85A_V33R.png', width = 4, height = 4, dpi = 300)
```
