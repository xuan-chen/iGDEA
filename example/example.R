rm(list= ls())

setwd("G:/Shared drives/BIOSTAT615_Group_Project/iGDEA")
source("R/main.R")

load("./data/GSE5583.RData")
obj = Create_DE_object(count)
obj@sample_info$genotype = c(rep("WT", 3), rep("KO", 3))
# str(obj)
obj = Filter_Data(obj, 0, 0)

gse_nb = Normalize_Data(obj, "none")
gse_nb = DEA_test(gse_nb, group = "genotype", option = "NB")
gse_nb = Get_DEG(gse_nb, "logfc", 0.5)
gse_nb = Get_DEG(gse_nb, "p", 0.05)

gse_t = Normalize_Data(obj)
gse_t = DEA_test(gse_t, group = "genotype", option = "T")
gse_t = Get_DEG(gse_t, "logfc", 0.5)
gse_t = Get_DEG(gse_t, "p", 0.05)

gse_w = Normalize_Data(obj, "none")
gse_w = DEA_test(gse_w, group = "genotype", option = "Wilcox")
gse_w = Get_DEG(gse_w, "logfc", 0.5)
gse_w = Get_DEG(gse_w, "p", 0.05)

save(gse_nb, gse_t, gse_w, file = "gse_result.RData")


######
rm(list=ls())
source("R/main.R")
load("./data/simulate_data.RData")
obj = Create_DE_object(countdata)
obj@sample_info$treatment = c(rep("A", 6), rep("B", 6))
str(obj)
obj = Filter_Data(obj, 0, 0)

sim_nb = Normalize_Data(obj, "none")
sim_nb = DEA_test(sim_nb, group = "treatment", option = "NB")
sim_nb = Get_DEG(sim_nb, "logfc", 0.5)
sim_nb = Get_DEG(sim_nb, "p", 0.05)

sim_t = Normalize_Data(obj, "none")
sim_t = DEA_test(sim_t, group = "treatment", option = "T")
sim_t = Get_DEG(sim_t, "logfc", 0.5)
sim_t = Get_DEG(sim_t, "p", 0.05)

sim_w = Normalize_Data(obj, "none")
sim_w = DEA_test(sim_w, group = "treatment", option = "Wilcox")
sim_w = Get_DEG(sim_w, "logfc", 0.5)
sim_w = Get_DEG(sim_w, "p", 0.05)

save(sim_nb, sim_t, sim_w, file = "sim_result.RData")

######
rm(list=ls())
load("./data/PBMC.RData")
count = as.matrix(ifnb@assays$RNA@counts)
source("R/main.R")

obj = Create_DE_object(count)
obj@sample_info$treatment = ifnb$stim
str(obj)
obj = Filter_Data(obj, 5, 5)

# sim_nb = Normalize_Data(obj, "none")
# sim_nb = DEA_test(sim_nb, group = "treatment", option = "NB")
# sim_nb = Get_DEG(sim_nb, "logfc", 0.5)
# sim_nb = Get_DEG(sim_nb, "p", 0.05)

pbmc_t = Normalize_Data(obj, "log")
pbmc_t = DEA_test(pbmc_t, group = "treatment", option = "T")
pbmc_t = Get_DEG(pbmc_t, "logfc", 1)
pbmc_t = Get_DEG(pbmc_t, "p", 0.05)
tmp_t = p.adjust(pbmc_t@DE_info$p, "bonferroni")
table(tmp_t < 0.05)
str(pbmc_t)

pbmc_w = Normalize_Data(obj, "log")
pbmc_w = DEA_test(pbmc_w, group = "treatment", option = "Wilcox")
pbmc_w = Get_DEG(pbmc_w, "logfc", 1)
pbmc_w = Get_DEG(pbmc_w, "p", 0.05)
tmp_w = p.adjust(pbmc_w@DE_info$p, "bonferroni")
table(tmp_w < 0.05)
str(pbmc_w)

pbmc_t_deinfo = pbmc_t@DE_info
pbmc_w_deinfo = pbmc_w@DE_info
save(pbmc_t_deinfo, pbmc_w_deinfo, tmp_t, tmp_w, file = "pbmc_result.RData")

