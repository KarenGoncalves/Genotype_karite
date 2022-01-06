library(tidyverse)

inds <- read.delim("input data - DO NOT MODIFY/inds", header = F)


population_information <- read.csv(file = "input data - DO NOT MODIFY/Population.csv")
names(population_information) <- c("Taxa", "Taxa.Name", "Country",
																	 "Region", "Location.Population",
																	 "Type.Land.Use")

population_information = 
	population_information[- which(is.na(population_information$Taxa)),]

pop_info <- 
	unique(population_information[
		, c("Country", "Region", "Location.Population", "Type.Land.Use")])


localities <- gsub("^[0-9]+_[CMT][a-z]{2}[A-Z][0-9]_([a-z]+)_([a-z]{2})_pl[0-9]+_[A-Z0-9]+", 
									 "\\1", inds$V1) %>% unique %>% sort()

locality_true_name <- sapply(localities, function(x) {
	pop_info[grep(x, pop_info$Location.Population, 
								ignore.case = T, value = F), c(2,3)] %>% 
		unique
}) %>% t() %>% as.data.frame()

locality_true_name["ktaba",] <- c("Ouest", "Koutaba")
locality_true_name["mayo",] <- c("Logone Oriental", "Mayongo")

land_use <- gsub("^[0-9]+_[CMT][a-z]{2}[A-Z][0-9]_([a-z]+)_([a-z]{2})_pl[0-9]+_[A-Z0-9]+", 
								 "\\2", inds$V1) %>% unique %>% sort()


land_use_true_name <- data.frame(
	truename = c("Champ agricole", "Champ", 
							 "Foret", "Jachare", "Jardin de case", 
							 "Savane", "Savane Naturelle", "Verger"),
	acronym = land_use)

population_information <- data.frame(matrix(nrow = 0, ncol = 5))
names(population_information) <- 
	c("inds", "Country", "Region", "Locality", "Land_use")

for (i in inds$V1) {
	tmp <- data.frame(matrix(nrow = 1, ncol = 5))
	names(tmp) <- 
		c("inds", "Country", "Region", "Locality", "Land_use")
	tmp$inds[1] = i
	tmp$Country[1] = 
		ifelse(
			grepl("^[0-9]+_Tch", i),
			"Tchad", ifelse(
				grepl("^[0-9]+_Cmr", i), 
				"Cameroun", ifelse(
					grepl("^[0-9]+_Mal", i), "Mali", "")
			))
	ind_loc <- gsub("^[0-9]+_[CMT][a-z]{2}[A-Z][0-9]_([a-z]+)_([a-z]{2})_pl[0-9]+_[A-Z0-9]+", 
									"\\1", i)
	
	ind_lu <- gsub("^[0-9]+_[CMT][a-z]{2}[A-Z][0-9]_([a-z]+)_([a-z]{2})_pl[0-9]+_[A-Z0-9]+", 
								 "\\2", i)
	tmp[1, c(3, 4)] = locality_true_name[ind_loc,]
	tmp[1,5] = land_use_true_name$truename[
		land_use_true_name$acronym == ind_lu]
	population_information = data.frame(rbind(population_information, tmp))
}

population_information <- apply(population_information, 2, as.character) %>%
	as.data.frame()
write_delim(file = "input data - DO NOT MODIFY/Whole_experiment_popinfo.txt", 
						x = population_information, 
						append = F, delim = "\t", col_names = T)

pop_map <- 
	data.frame(inds = inds$V1, 
						 pops = with(population_information, 
						 						paste(Country, Region, Locality, Land_use, sep = "_")))

write_delim(file = "input data - DO NOT MODIFY/popmap.txt", 
						x = pop_map, 
						append = F, delim = "\t", col_names = F)
