packages_to_use = c(
 "poppr", "mmod", "adegenet", "ape", "LEA",
 "treemap", "tidyverse",
 "vcfR", "pegas", "hierfstat",
 "ggpubr")
for (i in packages_to_use) { 
library(i, character.only = T)
}
setwd("D:\\postdoc_rprojects\\Genotype karite/")

rm(list=ls())
load("input_data_RefSNPs.RData")
groups <- with(population_information, 
  list(countries = unique(Country),
   regions = unique(Region),
   populations = unique(Locality),
   use = unique(Land_use)))
  
 # genetic_summary <- lapply(levels(pop(genind)), function(x) {
 # inds <- pop(genind) == x
 
 # basic.stats(genind[inds], diploid = T)
 # })
 # names(genetic_summary) <- levels(pop(genind))

# gen_summ_pops <- 
 # sapply(levels(pop(genind)), function(x) {
 # genetic_summary[[x]]$overall[c("Ho", "Hs", "Ht", "Fis")]
 # }) %>% t() %>% data.frame()
 
 # exp_het <- adegenet::Hs(genind, pop = pop(genind))
# exp_het_countries_mean <- 
 # sapply(groups$countries, function(x) {
 # a = exp_het[grep(x, levels(pop(genind)))] 
 # c(mean = mean(a), std = sd(a))
 # })
 
 # fis_countries_mean <- sapply(groups$countries, function(x) {
 # a = gen_summ_pops$Fis[grep(x, levels(pop(genind)))] 
 # c(mean = mean(a), std = sd(a))
 # })
 
# write.table(gen_summ_pops, file = "reference_snps_summary_genetics.txt", sep = "\t", dec = ".", quote = F, row.names = F, col.names = T)

# write.table(gen_summ_pops, file = "reference_snps_summary_genetics.txt", sep = "\t", dec = ".", quote = F, row.names = T, col.names = T)

# allecic_richness <- lapply(levels(pop(genind)), function(x) {
  # inds <- pop(genind) == x
 # data = data.frame(pop = x, tab(genind[inds]))
 # allelic.richness(data = data, diploid = T)
 # })

# for (i in names(allelic_richness)) {
	# if (!exists("boxplot_ar")) {
 		# boxplot_ar = data.frame(allelic_richness[[i]]$Ar)
		# names(boxplot_ar) = i
	# } else {
 		# a = ncol(boxplot_ar) + 1
	 	# boxplot_ar[, a] = allelic_richness[[i]]$Ar
	 	# names(boxplot_ar)[a] = i
# }}

# boxplot_ar$loci = rownames(boxplot_ar)

# pivot_longer(boxplot_ar, cols = !loci, 
	# names_to = "Populations", 
	# values_to = "Allelic_richness") %>%
 # ggplot(aes(Populations, Allelic_richness)) +
 # geom_boxplot() + 
 # theme(axis.text.x = element_text(angle = 90)) +
 # scale_y_continuous(limits = c(0.5, 3.5))


gff <- read.table("input data - DO NOT MODIFY/Vitpa_all.gff", 
									sep = "\t", quote = "", header = F)

correlations <- read.delim("input data - DO NOT MODIFY/correlations_vcf_assem.txt", header = F)
names(correlations) = c("assembly", "vcf")

for (i in 1:nrow(correlations)) {
	gff$V1 = gsub(correlations$assembly[i], correlations$vcf[i], gff$V1)
	gff$V9 = gsub(correlations$assembly[i], correlations$vcf[i], gff$V9)
}

vcf_separated <- lapply(correlations$vcf, function(x) {
file = paste0("input data - DO NOT MODIFY/separated_vcf_noimp/", x, ".vcf")
read.vcfR(file)
})

dna <- lapply(correlations$vcf, function(x) {
file = paste0("input data - DO NOT MODIFY/separated_fasta/", x, ".fasta")
a = read.dna(file,format = "fasta")
dimnames(a)[[1]][1] = x
a
})

names(dna) = names(vcf_separated) = correlations$vcf

for (x in correlations$vcf) {
print(dimnames(dna[[x]])[[1]][1])
}

chroms <- lapply(names(dna), function(x) {
create.chromR(name=x, vcf=vcf_separated[[x]], seq=dna[[x]], ann=gff[gff$V1 == x,])
})

plot(chrom)

chromoqc(chrom, dp.alpha=20)

head(chrom@var.info)

plot(chrom)


