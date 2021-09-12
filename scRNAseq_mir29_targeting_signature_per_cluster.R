
adult_ctrl = subset(mir29ko, subset = state == "Cre-")
adult_flox = subset(mir29ko, subset = state == "Cre+")
neo = subset(mir29ko, subset = state == "Neo")

actrl_norm = as.data.frame(as.matrix(GetAssayData(adult_ctrl, 'data')))
aflox_norm = as.data.frame(as.matrix(GetAssayData(adult_flox, 'data')))
neo_norm = as.data.frame(as.matrix(GetAssayData(neo, 'data')))

# for all miRNAs
# find target genes (names here) - context score < -0.2
# find the target genes for miRNAs ~50MB file
targets_mouse <- read.table("Targetscan_mouse_filtered_genes.txt")
colnames(targets_mouse) <- c("transcript_id_ver", "gene_name","seed", "mirna_name")
targets_mouse$transcript_id <- substr(targets_mouse$transcript_id_ver, 1, 18)


# find weak targets of miRNAs, used to exclude in background genes - context score >= -0.2
weak_targets_all <- read.table("Targetscan_mouse_weak_targets_genes_names.txt")
colnames(weak_targets_all) <- c("transcript_id_ver", "gene_name","seed", "mirna_name")
weak_targets_all$transcript_id <- substr(weak_targets_all$transcript_id_ver, 1, 18)


mir_family_info <- read.table("miR_Family_Info.txt", fill = T, header = T, sep = "\t")
mmu_mir_family_info <- mir_family_info[substr(mir_family_info[,'MiRBase.ID'] ,1 ,3) == "mmu", ]
mmu_con_mir_family_info <- mmu_mir_family_info[mmu_mir_family_info[,'Family.Conservation.'] >= 2, ]
seed_con <- mmu_con_mir_family_info[,c('Seed.m8', 'MiRBase.ID')]

library(plyr)
mirna_seeds_con <- ddply(seed_con, .(Seed.m8), summarize,
                         miRNA=paste(MiRBase.ID,collapse=";"))


colnames(mirna_seeds_con) <- c("seed", "miRNA")

targets_mouse_con_sites <- read.table("/home/hz543/data/targetscan/mouse/Targetscan_mouse_filtered_genes.txt")
colnames(targets_mouse_con_sites) <- c("transcript_id_ver", "gene_name","seed", "mirna_name")
targets_mouse_con_sites$transcript_id <- substr(targets_mouse_con_sites$transcript_id_ver, 1, 18)
targets_mouse_seed <- targets_mouse[targets_mouse_con_sites$'seed' %in% seed_con$Seed.m8,]
gene_w_con_mirna_sites <- targets_mouse_seed$gene_name
 


targeting_cluster_bg_adult <- function(cluster, seed, low_counts_cutoff){
    seed <- toupper(seed)
    targets_name <- targets_mouse[targets_mouse[,'seed'] == toString(seed),]$gene_name
    
    actrl_norm_subset = actrl_norm[, names(adult_ctrl$seurat_clusters)[adult_ctrl$seurat_clusters == cluster]]
    aflox_norm_subset = aflox_norm[, names(adult_flox$seurat_clusters)[adult_flox$seurat_clusters == cluster]]
    
    actrl_norm_filter <- actrl_norm_subset[rowMeans(actrl_norm_subset) > low_counts_cutoff,] # 10367 genes
    aflox_norm_filter <- aflox_norm_subset[rowMeans(aflox_norm_subset) > low_counts_cutoff,] # 10336
    
    actrl_norm_both <- actrl_norm_filter[rownames(actrl_norm_filter) %in% rownames(aflox_norm_filter),] # 10116 genes
    aflox_norm_both <- aflox_norm_filter[rownames(aflox_norm_filter) %in% rownames(actrl_norm_filter),]
    
    fold_change = log2(rowMeans(actrl_norm_both) / rowMeans(aflox_norm_both))
    fold_change_sites = fold_change[names(fold_change) %in% gene_w_con_mirna_sites]
        
    print(nrow(actrl_norm_both))
    if (length(targets_name) < 5){
        return("No enough targets")
    } else {
        targets_weak <- weak_targets_all[weak_targets_all[,'seed'] == toString(seed),]$gene_name
        exp <- fold_change_sites[names(fold_change_sites) %in% targets_name]
        exp_background <- fold_change_sites[!names(fold_change_sites) %in% targets_name]
        exp_background_no_weak <- exp_background[!names(exp_background) %in% targets_weak]
    
        # wilcox test
        wilcox_test <- wilcox.test(exp, exp_background_no_weak)
        p_value <- wilcox_test$p.value
        p_value_sign = sign(median(exp) - median(exp_background_no_weak)) * p_value
    
        return (p_value_sign)
    }
}

#targeting_bg(0, 'AGCACCA')
cluster_target_pvals_bg = c()
for (i in c(0:8)){
    print(i)
    cluster_target_pvals_bg = c(cluster_target_pvals_bg, targeting_cluster_bg_adult(i, 'AGCACCA', 0.01))
}

