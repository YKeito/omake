#"~/Nakano_RNAseq/network_analysis/script/omake/GL3GR_list_fisher.R"
Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/B.cinerea_newlist.txt", sep = "\t", header = T, stringsAsFactors = F)
Botrytis_cinerea$Gene.Locus <- toupper(Botrytis_cinerea$Gene.Locus)
PstDC3000 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
PstDC3000_hrp <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/PstDC3000_hrp-_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
                            
GL3GR_DEGs_4h <- araiGL3GR_FC[araiGL3GR_FC$X4h_q_value < 0.05, c("test_id", "D4M4_log2")]
GL3GR_DEGs_8h <- araiGL3GR_FC[araiGL3GR_FC$X8h_q_value < 0.05, c("test_id", "D8M8_log2")]
GL3GR_DEGs_12h <- araiGL3GR_FC[araiGL3GR_FC$X12h_q_value < 0.05, c("test_id", "D12M12_log2")]
GL3GR_DEGs_16h <- araiGL3GR_FC[araiGL3GR_FC$X16h_q_value < 0.05, c("test_id", "D16M16_log2")]
GL3GR_DEGs_20h <- araiGL3GR_FC[araiGL3GR_FC$X20h_q_value < 0.05, c("test_id", "D20M20_log2")]
GL3GR_DEGs_24h <- araiGL3GR_FC[araiGL3GR_FC$X24h_q_value < 0.05, c("test_id", "D24M24_log2")]
GL3GR_DEGs_48h <- araiGL3GR_FC[araiGL3GR_FC$X48h_q_value < 0.05, c("test_id", "D48M48_log2")]
GL3GR_DEGs <- list(GL3GR_DEGs_4h$test_id, GL3GR_DEGs_8h$test_id, GL3GR_DEGs_12h$test_id, GL3GR_DEGs_16h$test_id, GL3GR_DEGs_20h$test_id, GL3GR_DEGs_24h$test_id, GL3GR_DEGs_48h$test_id)
GL3GR_DEGs <- unique(unlist(GL3GR_DEGs))

T_AGI <- list(CY15_AGI, CY16_AGI, CY20_AGI, Botrytis_cinerea$Gene.Locus, PstDC3000$pseudomonas.syringae, PstDC3000_hrp$AGI, MeJA_DEGs$AGI, BTH$AGI)
obname <- c("CY15", "CY16", "CY20", "B.cinerea", "PstDC3000", "PstDC3000_hrp-", "MeJA", "BTH")
T_pvalue <- c()
test <- c()
data <- list()
i <- 1
for(i in i:length(T_AGI)){
  a <- intersect(GL3GR_DEGs, T_AGI[[i]])
  b <- setdiff(GL3GR_DEGs, T_AGI[[i]])
  c <- setdiff(T_AGI[[i]], GL3GR_DEGs)
  d <- nrow(CY15_Table_20180801)-length(a)-length(b)-length(c)
  mx <- matrix(c(length(a), length(b), length(c), d), nrow=2, byrow=T)
  T_pvalue <- c(T_pvalue, fisher.test(mx)$p.value)
  names(T_pvalue)[i] <- obname[i]
  colnames(mx) <- c(paste0(obname[i], "_Yes"), paste0(obname[i], "_No"))
  rownames(mx) <- c(paste0("GL3GR_list", "_Yes"), paste0("negishi_list", "_No"))
  test <- cbind(test, mx)
  data <- c(data, list(a))
  i <- i+1
}
write.table(test, file = "~/Nakano_RNAseq/network_analysis/result_Table/arai_list/GL3GR_fisher.txt", append=F, quote = F, sep = "\t", row.names = T, col.names = T)


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
write.table(temp, file = "~/Nakano_RNAseq/network_analysis/result_Table/arai_list/GL3GR_intersection.txt", append=F, quote = F, sep = "\t", row.names = F, col.names = T)
