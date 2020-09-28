"~/Nakano_RNAseq/network_analysis/script/omake/baratuki_list_fisher.R"
CY_annotation <- read.csv("~/Nakano_RNAseq/network_analysis/base/CY_annotation.csv")
CY_annotation$AGI <- toupper(CY_annotation$AGI)
baratuki <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/negishi_baratuki/group_normCV3.txt", sep = "\t", header = T, stringsAsFactors = F)
CY15_AGI <- rownames(allRNASeq[allRNASeq$CY15_1h_q_value < 0.05 | allRNASeq$CY15_3h_q_value < 0.05 | allRNASeq$CY15_12h_q_value < 0.05 | allRNASeq$CY15_24h_q_value < 0.05, ])
CY16_AGI <- rownames(allRNASeq[allRNASeq$CY16_1h_q_value < 0.05 | allRNASeq$CY16_3h_q_value < 0.05 | allRNASeq$CY16_12h_q_value < 0.05 | allRNASeq$CY16_24h_q_value < 0.05, ])
CY20_AGI <- rownames(allRNASeq[allRNASeq$CY20_1h_q_value < 0.05 | allRNASeq$CY20_3h_q_value < 0.05 | allRNASeq$CY20_12h_q_value < 0.05 | allRNASeq$CY20_24h_q_value < 0.05, ])
Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/B.cinerea.txt", sep = "\t", header = T, stringsAsFactors = F)
PstDC3000 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
PstDC3000_hrp <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/PstDC3000_hrp-_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
MeJA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/MeJA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/TF/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)

T_AGI <- list(CY15_AGI, CY16_AGI, CY20_AGI, Botrytis_cinerea$B.cenerera, PstDC3000$pseudomonas.syringae, PstDC3000_hrp$AGI, MeJA_DEGs$AGI, BTH$AGI)
obname <- c("CY15", "CY16", "CY20", "B.cinerea", "PstDC3000", "PstDC3000_hrp-", "MeJA", "BTH")
T_pvalue <- c()
test <- c()
data <- list()
i <- 1
for(i in i:length(T_AGI)){
  a <- intersect(baratuki$AGI, T_AGI[[i]])
  b <- setdiff(baratuki$AGI, T_AGI[[i]])
  c <- setdiff(T_AGI[[i]], baratuki$AGI)
  d <- nrow(CY15_Table_20180801)-length(a)-length(b)-length(c)
  mx <- matrix(c(length(a), length(b), length(c), d), nrow=2, byrow=T)
  T_pvalue <- c(T_pvalue, fisher.test(mx)$p.value)
  names(T_pvalue)[i] <- obname[i]
  colnames(mx) <- c(paste0(obname[i], "_Yes"), paste0(obname[i], "_No"))
  rownames(mx) <- c(paste0("negishi_list", "_Yes"), paste0("negishi_list", "_No"))
  test <- cbind(test, mx)
  data <- c(data, list(a))
  i <- i+1
}
write.table(test, file = "~/Nakano_RNAseq/network_analysis/result_Table/negishi_list/negishi_fisher.txt", append=F, quote = F, sep = "\t", row.names = T, col.names = T)

temp <- matrix(rep("No"), nrow = length(unique(unlist(data))), ncol = c(length(obname)+2))
colnames(temp) <- c("AGI", "TF", obname)
temp[, "AGI"] <- unique(unlist(data))
a <- intersect(TF_family$AGI, unique(unlist(data)))
temp[match(a, unique(unlist(data))), "TF"] <- TF_family$TF[match(a, TF_family$AGI)]
n <- 3
for(n in n:ncol(temp)){
  temp[match(data[[n-2]], temp[, "AGI"]), colnames(temp)[n]] <- "Yes"
  print(n)
  n <- n+1
}
write.table(temp, file = "~/Nakano_RNAseq/network_analysis/result_Table/negishi_list/intersection.txt", append=F, quote = F, sep = "\t", row.names = F, col.names = T)
