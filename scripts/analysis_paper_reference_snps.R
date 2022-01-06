########## Load packages 
packages_to_use = c(
	"poppr", "mmod", "adegenet", "ape", "LEA",
	"treemap", "tidyverse", "gridExtra",
	"vcfR", "pegas", "hierfstat", "ggdendro",
	"ggpubr", "RColorBrewer")

# install.packages(packages_to_use,
# 								 dependencies = TRUE)
devtools::install_github("kwstat/pals")
for (i in packages_to_use) { 
library(i, character.only = T)}

library(pals)


######### Load pop info

population_information <- read.csv(file = "input data - DO NOT MODIFY/Population.csv")
names(population_information) <- c("Taxa", "Taxa.Name", "Country",
																	 "Region", "Location.Population",
																	 "Type.Land.Use")

population_information = 
	population_information[- which(is.na(population_information$Taxa)),]



########## Load data
vcf <- read.vcfR("input data - DO NOT MODIFY/third_filters_m3_p50_x0_S2.singleton.unlinked_100k.vcf.imputed_Q14.vcf")

genind <- vcfR2genind(vcf, sep = "/")
ind_new_data <- data.frame(individuals = rownames(tab(genind)))
pop_info <- 
	unique(population_information[
		, c("Country", "Region", "Location.Population", "Type.Land.Use")])


localities <- gsub("^[CMT][a-z]{2}[A-Z][0-9]-([a-z]+)_([a-z]{2})-[a-z0-9A-Z\\-]+", 
									 "\\1", ind_new_data$individuals) %>% unique %>% sort()

locality_true_name <- sapply(localities, function(x) {
	pop_info[grep(x, pop_info$Location.Population, 
								ignore.case = T, value = F), c(2,3)] %>% 
		unique
}) %>% t() %>% as.data.frame()

locality_true_name["ktaba",] <- c("Ouest", "Koutaba")
locality_true_name["mayo",] <- c("Logone Oriental", "Mayongo")

land_use <- gsub("^[CMT][a-z]{2}[A-Z][0-9]-([a-z]+)_([a-z]{2})-[a-z0-9A-Z\\-]+", 
								 "\\2", ind_new_data$individuals) %>% unique %>% sort()


land_use_true_name <- data.frame(
	truename = c("Champ agricole", "Champ", 
							 "Foret", "Jachare", "Jardin de case", 
							 "Savane Naturelle", "Verger"),
	acronym = land_use)

population_information <- data.frame(matrix(nrow = 0, ncol = 5))
names(population_information) <- 
	c("inds", "Country", "Region", "Locality", "Land_use")

for (i in ind_new_data$individuals) {
	tmp <- data.frame(matrix(nrow = 1, ncol = 5))
	names(tmp) <- 
		c("inds", "Country", "Region", "Locality", "Land_use")
	tmp$inds[1] = i
	tmp$Country[1] = 
		ifelse(
			grepl("^Tch", i),
			"Tchad", ifelse(
				grepl("Cmr", i), 
				"Cameroun", ifelse(
					grepl("Mal", i), "Mali", "")
			))
	ind_loc <- gsub("^[CMT][a-z]{2}[A-Z][0-9]-([a-z]+)_([a-z]{2})-[a-z0-9A-Z\\-]+", 
									"\\1", i)
	
	ind_lu <- gsub("^[CMT][a-z]{2}[A-Z][0-9]-([a-z]+)_([a-z]{2})-[a-z0-9A-Z\\-]+", 
								 "\\2", i)
	tmp[1, c(3, 4)] = locality_true_name[ind_loc,]
	tmp[1,5] = land_use_true_name$truename[
		land_use_true_name$acronym == ind_lu]
	population_information = data.frame(rbind(population_information, tmp))
}

strata(genind) <- population_information
setPop(genind) <-
	~Country/Region/Locality


############ Pairwise Fst

groups <- with(population_information, 
							 list(countries = unique(Country),
							 		 regions = unique(Region),
							 		 populations = unique(Locality),
							 		 use = unique(Land_use)))

pwfst <- genet.dist(genind, diploid = T, method = "Fst")
nj_fst <- pwfst %>% nj()
nj_fst$edge.length <- abs(nj_fst$edge.length)
nj_fst$tip.label <- 
	levels(pop(genind)) %>%
	gsub(pattern = "_", replacement = " - ")

pwfst_df <- as.matrix(pwfst) %>% 
	as.data.frame()
names(pwfst_df) <- rownames(pwfst_df) <- 
	levels(pop(genind))

dend_pwfst <- hclust(pwfst_df %>% as.matrix %>% as.dist)
pops_countries <- dend_pwfst$labels %>% 
	str_split(pattern = "_") %>% 
	lapply(function(x) x[1]) %>%
	unlist %>% as.factor

############# Pairwise Fst dend plot

ggdend_pwfst <- as.dendrogram(dend_pwfst) %>%
	dendro_data()

ggdend_pwfst$labels$label <- 
	gsub("_", " - ", ggdend_pwfst$labels$label)

ggdend_pwfst$labels$country <- 
	ggdend_pwfst$labels$label %>% str_split(pattern = " ") %>%
	lapply(function(x) x[1]) %>% unlist %>% as.factor 

ggplot() +
	geom_segment(data = ggdend_pwfst$segments, 
							 aes(x = y, y = x, 
							 		xend = yend, yend = xend))+
	geom_text(data = ggdend_pwfst$labels, 
						aes(y, x, label = label, color = country),
						hjust = 1.1, size = 3) +
	scale_x_continuous(
		limits = c(-0.25, 0.18), 
		breaks = seq(0, 0.2, 0.05)) + 
	theme_minimal() + 
	theme(axis.text.y = element_blank(), 
				panel.grid = element_blank()) + 
	labs(x = "Pairwise Fst", y = "")

ggsave(paste0("Reference_snps_PairwiseFST_njtree.pdf"))
knitr::kable(pwfst_df, digits = 3)


########## Summary per country
country_genind <- genind
setPop(country_genind) <-
	~Country

genetic_summary_countries <- 
	lapply(groups$countries, function(x) {
		inds <- pop(country_genind) == x
		
		basic.stats(country_genind[inds])
	})

names(genetic_summary_countries) <- groups$countries


exp_het_countries <- adegenet::Hs(country_genind, 
																	pop = pop(country_genind))

sapply(groups$countries, function(x) {
	genetic_summary_countries[[x]]$overall[c("Ho", "Hs", "Ht", "Fis")]
}) %>% t() %>% data.frame()

exp_het_countries

############### summary per pop

genetic_summary <- lapply(levels(pop(genind)), function(x) {
	inds <- pop(genind) == x
	
	basic.stats(genind[inds], diploid = T)
})
names(genetic_summary) <- levels(pop(genind))


gen_summ_pops <- 
	sapply(levels(pop(genind)), function(x) {
		genetic_summary[[x]]$overall[c("Ho", "Hs", "Ht", "Fis")]
	}) %>% t() %>% data.frame()

exp_het <- adegenet::Hs(genind, pop = pop(genind))

exp_het_countries_mean <- 
	sapply(groups$countries, function(x) {
		a = exp_het[grep(x, levels(pop(genind)))] 
		c(mean = mean(a), std = sd(a))
	})

fis_countries_mean <- sapply(groups$countries, function(x) {
	a = gen_summ_pops$Fis[grep(x, levels(pop(genind)))] 
	c(mean = mean(a), std = sd(a))
})
exp_het_countries_mean
fis_countries_mean

################# Amova with Mali
amova_mali <- poppr.amova(x = genind, within = F,
													hier = ~Country/Region/Locality,
													threads = parallel::detectCores(), cutoff = 0.95
)

(rand_amova_mali <- randtest(amova_mali, nrepet = 999))

plot(rand_amova_mali)

knitr::kable(amova_mali$results, caption = 
						 	"Result AMOVA including the Mali population.\
						 Hierarchy: _Country - Region - Locality_.\
						 Sample = Location.Population")

knitr::kable(amova_mali$componentsofcovariance, caption = 
						 	"Components of covariance from AMOVA including the Mali population.\
						 Hierarchy: _Country - Region - Locality_.\
						 Sample = Location.Population")


##################### Amova without Mali
exclude_mali <- population_information$inds[
	population_information$Country != "Mali"]

genind_nomali <- genind[exclude_mali]
strata(genind_nomali) <- population_information[
	population_information$Country != "Mali",] %>%
	apply(MARGIN = 2, FUN = as.character) %>% as.data.frame
setPop(genind_nomali) <-
	~Country/Region/Locality

amova_nomali <- poppr.amova(x = genind_nomali, within = F,
														hier = ~Country/Region/Locality,
														threads = parallel::detectCores(), cutoff = 0.95
)

(rand_amova_nomali <- randtest(amova_nomali, nrepet = 999))

plot(rand_amova_nomali)

knitr::kable(amova_nomali$results, caption = 
						 	"Result AMOVA excluding the Mali population.\
						 Hierarchy: _Country - Region - Locality_.\
						 Sample = Location.Population")


################## glPCA

genlight <- vcfR2genlight(vcf, n.cores = parallel::detectCores())

glPCA_result <- glPca(
	genlight, parallel = T, nf = 20,
	n.cores = parallel::detectCores()
)

glPCA_df <- data.frame(
	Taxa.Name = rownames(glPCA_result$scores),
	glPCA_result$scores
) %>% merge(population_information[, -1], by = "inds") 

#devtools::install_github("jaredhuling/jcolors")
library("jcolors")

ggplot(glPCA_df, 
			 aes(PC1, PC2,shape = Country, color = Region)) +
	geom_point() + scale_color_jcolors(palette = "pal8")

# DAPC

## DAPC using clusters obtained with find.clusters function

############### DAPC with find.clusters
theme_set(theme_bw())
strata(genind) <- population_information
setPop(genind) <-
	~Country/Region/Locality

maxK <- 20
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
	grp <- find.clusters(genind, n.pca = 490, 
											 choose.n.clust = FALSE,  max.n.clust = maxK)
	myMat[i,] <- grp$Kstat
}

myMat <- as.data.frame(myMat)
myMat$Group = 1:nrow(myMat)
my_df <- pivot_longer(myMat, cols = !Group, names_to = "K", values_to = "BIC")
my_df$Group <- factor(my_df$Group, levels = 1:10)
my_df$K <- factor(my_df$K, levels = 1:20)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot() + 
	p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

mean_BIC <- apply(myMat[, -21], 2, mean)
nK = (mean_BIC == min(mean_BIC)) %>% which


############ DAPC with find.clusters

clusters_karite <- find.clusters(
	genind, n.pca = 490, stat = "BIC",
	method = "ward", n.clust = nK)

table(pop(genind), clusters_karite$grp)

table.value(
	table(population_information$Country, clusters_karite$grp),
	col.labels = paste("K ", levels(clusters_karite$grp)))

dapc_countries <- dapc(genind, 
											 clusters_karite$grp, 
											 n.pca = 90, n.da = 40)


dapc_da <- data.frame(
	population_information[, 2:5],
	clusters = dapc_countries$grp,
	dapc_countries$ind.coord
)



########## DAPC scatter plot

theme_set(theme_bw())
scatter_dapc <- list(
	eigen = (data.frame(
		DAs = factor(paste0("Discriminant function ", 1:dapc_countries$n.da),
								 levels = paste0("Discriminant function ", 1:dapc_countries$n.da)),
		eigenvalues = dapc_countries$eig) %>%
			ggplot(aes(DAs, eigenvalues, fill = eigenvalues)) +
			geom_col() + 
			labs(x = "", y = "Eigenvalues", 
					 title = "Eigenvalues of retained discriminant functions") +
			#scale_y_continuous(breaks = seq(0, 6000, 1000)) +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
						legend.position = "none")
	), 
	ld1_ld2 =
		(ggplot(dapc_da, aes(LD1, LD2, color = clusters)) + 
		 	geom_point(aes(shape = Country)) + 
		 	stat_ellipse(aes(color = clusters)) +
		 	labs(x = "DAPC LD1", y = "DAPC LD2") +
		 	scale_color_manual(values = pals::glasbey(nK))),
	ld1_ld3 =
		(ggplot(dapc_da, aes(LD1, LD3, color = clusters)) + 
		 	geom_point(aes(shape = Country)) + 
		 	stat_ellipse(aes(color = clusters)) +
		 	labs(x = "DAPC LD1", y = "DAPC LD3") +
		 	scale_color_manual(values = pals::glasbey(nK))),
	ld2_ld3 =
		(ggplot(dapc_da, aes(LD2, LD3, color = clusters)) + 
		 	geom_point(aes(shape = Country)) + 
		 	stat_ellipse(aes(color = clusters)) +
		 	labs(x = "DAPC LD2", y = "DAPC LD3") +
		 	scale_color_manual(values = pals::glasbey(nK)))
)

ggsave(
	filename = paste0("ReferenceSNP_Scatters_DAPC_findclusters", Sys.Date(), ".pdf"), 
	plot = marrangeGrob(scatter_dapc, nrow=1, ncol=1)
)

save(scatter_dapc, "Scatters_DAPC_findclusters.RData")

########## DAPC PCA plot ###########

dapc_pca <- 
	data.frame(inds = population_information$inds, 
					 grp = dapc_countries$grp, 
					 Country = population_information$Country,
					 PC1 = dapc_countries$tab$`PCA-pc.1`, 
					 PC2 = dapc_countries$tab$`PCA-pc.2`) %>% 
	ggplot(aes(PC1, PC2, color = grp, shape = Country)) + geom_point() +
	scale_color_manual(values = pals::glasbey(nK))


ggsave(dapc_pca,
	filename = paste0("ReferenceSNP_PCA_DAPC_findclusters", Sys.Date(), ".pdf")
)

save(dapc_pca, "PCA_DAPC_findclusters.RData")
########### STRUCTURE-like plot ###############

setPop(genind) <- ~Country/Region/Locality
membership_probability = data.frame(
	individuals = rownames(dapc_countries[["posterior"]]),
	dapc_countries[["posterior"]]
) %>% merge(
	data.frame(individuals = population_information$inds,
						 pops = as.character(pop(genind))),
	by = "individuals") %>%
	pivot_longer(
		-c(individuals, pops), names_to = "Clusters", 
		values_to = "probability", names_prefix = "X")

membership_probability$individuals <- 
	factor(membership_probability$individuals,
				 levels = population_information$inds[
				 	order(population_information$inds)])

mem_prob_pops <- lapply(levels(pop(genind)), function(i) {
	title = paste0(
		i, "\nN=", which(pop(genind) == i) %>% length
	)
	filter(membership_probability, pops == i) %>%
		ggplot(
			aes(individuals, probability, fill = Clusters)) +
		geom_col() + labs(x = "", y = "",
											title = title) +
		scale_fill_manual(values = glasbey(nK), breaks = 1:nK) +
		theme(axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					axis.title.y = element_text(vjust = 0.5),
					panel.grid = element_blank(), 
					plot.title = element_text(hjust = 0.5)
		)
})

plot_groups <- list(
	ggarrange(plotlist = lapply(1:6, function(i) {mem_prob_pops[[i]]}),
						ncol = 2, nrow = 3, labels = "AUTO", 
						common.legend = T, legend = "bottom"),
	ggarrange(plotlist = lapply(7:12, function(i) {mem_prob_pops[[i]]}),
						ncol = 2, nrow = 3, labels = "AUTO", 
						common.legend = T, legend = "bottom"),
	ggarrange(plotlist = lapply(13:18, function(i) {mem_prob_pops[[i]]}),
						ncol = 2, nrow = 3, labels = "AUTO", 
						common.legend = T, legend = "bottom"),
	ggarrange(plotlist = lapply(19:21, function(i) {mem_prob_pops[[i]]}),
						ncol = 1, nrow = 3, labels = "AUTO", 
						common.legend = T, legend = "bottom"))
library(gridExtra)

ggsave(
	filename = paste0("ReferenceSNP_Cluster_membership_DAPC_DA10_PC90_k", nK, Sys.Date(), ".pdf"), 
	plot = marrangeGrob(plot_groups, nrow=1, ncol=1)
)

save(plot_groups,  file = "plots_clusterMembership.RData")

############### poppr analysis #############
poppr_idx <- 
		poppr(genind)

poppr_idx_df <- sapply(1:length(poppr_idx), function(x) {
	a = as.data.frame(poppr_idx[[x]]) 
	names(a) <- names(poppr_idx[[x]])
	a
})  %>% t() 


save(poppr_idx, poppr_idx_df, file = "poppr_results.RData")