# === Lab 7: Statistical Analysis of Gene Expression Data in R ===
# written by: Nilesh Domah    domah001    5216096

# === Part 1: Understanding the Breast Cancer Study ===
# please refer to the text file in the zip folder


# === Part 2: Loading and Manipulating Expression Data in R ===

# a) 
gene_expression <- as.matrix(read.table("BC_MicroArray.txt"), header = T, row.names = 1, sep = "\t", quote = "")
gene_status <- as.matrix(read.table("BC_MicroArray_status.txt"), header = T, row.names = 1, sep = "\t", quote = "")

# b) 
str(gene_expression)
str(gene_status)

# c) 
BRCA1 <- gene_expression["BRCA1", ]
print(BRCA1[2])

# d) 
GSM <- gene_expression[, "GSM519791"]
print(GSM[1])

# e) 
BRCA1_mean <- mean(BRCA1)
BRCA1_sd <- sd(BRCA1)

print(BRCA1_mean)
print(BRCA1_sd)

# f) 
GSM_mean <- mean(GSM)
GSM_sd <- sd(GSM)

# g) 
print(sort(GSM, decreasing = TRUE)[1:10])

# h) 
ERpos_samples <- gene_status[gene_status[ ,3] == "1", 2]
print(ERpos_samples)

# i) 
ERneg_samples <- gene_status[gene_status[ ,3] == "0", 2]
print(ERneg_samples)

# j) 
ER_pos <- (gene_expression[, ERpos_samples])
ER_neg <- (gene_expression[, ERneg_samples])

ER_pos_mean <- rowMeans(ER_pos)
ER_neg_mean <- rowMeans(ER_neg)

expr_difference <- ER_pos_mean - ER_neg_mean
print(expr_difference)

# k) 
expr_difference <- sort(expr_difference, decreasing = TRUE)
gene_names <- names(expr_difference)

large_pos_diff <- gene_names[1]
print(large_pos_diff)

small_pos_diff <- gene_names[length(gene_names)]
print(small_pos_diff)



# === Part 3: Statistical Analysis of Microarray Gene Expresson Data ===

# a)
numerator <- ER_pos_mean - ER_neg_mean

s2x1 <- apply(ER_pos, 1, var)
s2x2 <- apply(ER_neg, 1, var)

denominator <- sqrt((s2x1/43) + (s2x2/45))

t_stat <- numerator / denominator
print(t_stat)

# b) 
hist(t_stat, breaks = 100)

# c)
pval <- 2 * pt(-abs(t_stat) ,86)
print(pval)

# d)
# please refer to the text file in the zip folder

# e-i)
pval_sort <- sort(t_stat, index.return = TRUE)
sig_pval_sort <- pval_sort$x[pval_sort$x < 0.05]
pval_sort_stat <- t_stat[pval_sort$ix[pval_sort$x < 0.05]]
pval_sort_gene <- names(sig_pval_sort)
pval_sort_list <- list(genes = pval_sort_gene, t_statistics = pval_sort_stat, p_value = sig_pval_sort)
print(pval_sort_list)

#e-ii)
print(length(pval_sort_list$p_value))


bonf_pval <- 2 * pt(-abs(t_stat) ,32862)
bonf_pval_sort <- sort(t_stat, index.return = TRUE)
bonf_sig_pval_sort <- bonf_pval_sort$x[bonf_pval_sort$x < 0.05]
bonf_pval_sort_stat <- t_stat[bonf_pval_sort$ix[bonf_pval_sort$x < 0.05]]
bonf_pval_sort_gene <- names(bonf_sig_pval_sort)
bonf_pval_sort_list <- list(genes = bonf_pval_sort_gene, t_statistics = bonf_pval_sort_stat, p_value = bonf_sig_pval_sort)

print(length(bonf_pval_sort_list$p_value))

#e-iii)
print(length(pval_sort_list$p_value > mean(bonf_sig_pval_sort)))

#e-iv)
print(length(pval_sort_list$p_value < mean(bonf_sig_pval_sort)))

#e-v)
# please refer to the text file in the zip folder



# please refer to the rest of this question in text file in the zip folder

# === Part 4: Using R Packages to Process RNA-Seq Data ===
# please refer to Part4_updated_version2.R





