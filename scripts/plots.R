#!/bin/env Rscript

#################### SETUP ####################

### Environment
library(dplyr)
require(data.table)
library(argparse)

### Ex. round_any(132, 10, f = ceiling) should output 140 (can also use f = round)
round_ceiling = function(x, accuracy, f = ceiling){
	f(x / accuracy) * accuracy
}

### Capitalize first letter of a given string
caps <- function(x) {
	paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

### Syndromic flag accepted entries
usage <- function() {
	help.text <- "Syndromes accepted include: 
	all (default)
	respir (or respiratory)
	blood
	gastro (or gastrointestinal)
	human (for all human-relevant virus families)
	csf (meaning cerebral spine fluid)
	urine
	sjf (meaning skin or joint fluid)
	special (referring to special pathogens, zoonotics, vectorborne)
	
	"
	return(help.text)
}

### All syndromic lists (from syndromic_filter.py)
### Consider: examining the "Unclassified"?
#142
all <- list("Ackermannviridae","Adenoviridae","Adomaviridae","Alloherpesviridae","Alphaflexiviridae","Alphasatellitidae","Alphatetraviridae","Alvernaviridae","Amalgaviridae","Amnoonviridae","Ampullaviridae","Anelloviridae","Arenaviridae",
		"Arteriviridae","Ascoviridae","Asfarviridae","Aspiviridae","Astroviridae","Autolykiviridae","Bacilladnaviridae","Baculoviridae","Barnaviridae","Benyviridae","Betaflexiviridae",
		"Bicaudaviridae","Bidnaviridae","Birnaviridae","Bornaviridae","Bromoviridae","Caliciviridae","Carmotetraviridae","Caulimoviridae","Chrysoviridae","Chuviridae","Circoviridae","Clavaviridae","Closteroviridae","Coronaviridae",
		"Corticoviridae","Cruciviridae","Cruliviridae","Cystoviridae","Deltaflexiviridae","Dicistroviridae","Endornaviridae","Euroniviridae","Filoviridae","Fimoviridae","Flaviviridae","Flexiviridae","Fusariviridae","Fuselloviridae","Gammaflexiviridae",
		"Geminiviridae","Genomoviridae","Globuloviridae","Guttaviridae","Hantaviridae","Hepadnaviridae","Hepeviridae","Herpesviridae","Hypoviridae","Hytrosaviridae","Iflaviridae","Inoviridae","Iridoviridae","Lavidaviridae",
		"Leviviridae","Lipothrixviridae","Luteoviridae","Malacoherpesviridae","Marnaviridae","Marseilleviridae","Medioniviridae","Megabirnaviridae","Mesoniviridae","Metaviridae","Microviridae","Mimiviridae","Mononiviridae","Mymonaviridae","Myoviridae",
		"Mypoviridae","Nairoviridae","Nanoviridae","Narnaviridae","Nimaviridae","Nodaviridae","Nudiviridae","Nyamiviridae","Orthomyxoviridae","Papillomaviridae","Paramyxoviridae","Partitiviridae",
		"Parvoviridae","Peribunyaviridae","Permutotetraviridae","Phasmaviridae","Phenuiviridae","Phycodnaviridae","Picobirnaviridae","Picornaviridae","Pithoviridae","Plasmaviridae",
		"Pleolipoviridae","Pneumoviridae","Podoviridae","Polycipiviridae","Polydnaviridae","Polyomaviridae","Portogloboviridae","Potyviridae","Poxviridae","Qinviridae","Quadriviridae","Reoviridae","Retroviridae","Rhabdoviridae","Roniviridae",
		"Rudiviridae","Sarthroviridae","Secoviridae","Siphoviridae","Smacoviridae","Solemoviridae","Solinviviridae","Sphaerolipoviridae","Sphaeromadae","Spiraviridae","Sunviridae",
		"Tectiviridae","Tobaniviridae","Togaviridae","Tolecusatellitidae","Tombusviridae","Totiviridae","Tristromaviridae","Turriviridae","Tymoviridae","Virgaviridae","Wupedeviridae","Yueviridae")
#20
respir <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Reoviridae", "Togaviridae")

#20
gastro <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Caliciviridae", "Circoviridae", "Coronaviridae", "Filoviridae", "Hantaviridae",
		"Hepadnaviridae", "Hepeviridae", "Nairoviridae", "Orthomyxoviridae", "Partitiviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae", "Polyomaviridae",
		"Reoviridae", "Tobaniviridae")

#26
blood <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae",
		"Hepeviridae", "Herpesviridae", "Nairoviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae",
		"Pneumoviridae", "Polyomaviridae", "Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Togaviridae")

#33
human <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Astroviridae", "Bornaviridae", "Caliciviridae", "Circoviridae", "Coronaviridae",
		"Filoviridae", "Flaviviridae", "Hantaviridae", "Hepadnaviridae", "Hepeviridae", "Herpesviridae", "Nairoviridae", "Orthomyxoviridae", "Papillomaviridae",
		"Paramyxoviridae", "Partitiviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picobirnaviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae",
		"Poxviridae", "Reoviridae", "Retroviridae", "Rhabdoviridae", "Smacoviridae", "Tobaniviridae", "Togaviridae")

#21
csf <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Orthomyxoviridae", "Paramyxoviridae", "Parvoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae")

#11
urine <- list("Adenoviridae", "Anelloviridae", "Arenaviridae", "Flaviviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Poxviridae", "Rhabdoviridae", "Togaviridae")

#9
sjf <- list("Flaviviridae", "Herpesviridae", "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Picornaviridae", "Polyomaviridae", "Rhabdoviridae", "Togaviridae")

#20
special <- list("Arenaviridae", "Arteriviridae", "Bornaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Hantaviridae", "Hepeviridae", "Herpesviridae",
		"Nairoviridae", "Paramyxoviridae", "Peribunyaviridae", "Phenuiviridae", "Picornaviridae", "Pneumoviridae", "Poxviridae", "Reoviridae", "Retroviridae",
		"Rhabdoviridae", "Togaviridae")
 
########## Breakdown of families ##########

family_proportion_plot <- function(input, output_dir) {

	data_subset <- data.frame(row.names = input$V1, input$V2)
	lbls <- rownames(data_subset)
	pct <- round(as.numeric(sub("%", "", data_subset$input.V2)), digits = 2)
	lbls <- paste(lbls, pct)
	lbls <- paste(lbls, "%", sep="")

	png(filename = file.path(output_dir, "breakdown_ref_families.png"), width=8, height=6, units="in", res=320)
	pie(pct, labels = lbls, cex = 0.7, radius = 1.05, col=rainbow(length(lbls)))
	dev.off()

}

########## Breakdown of balitmore classification ##########

baltimore_proportion_plot <- function(input, output_dir) {

	data_subset <- data.frame(row.names = input$V1, input$V2)
	lbls <- rownames(data_subset)
	pct <- round(as.numeric(sub("%", "", data_subset$input.V2)), digits = 2)
	lbls <- paste(lbls, pct)
	lbls <- paste(lbls, "%", sep="")

	png(filename = file.path(output_dir, "breakdown_ref_balitmore.png"), width=8, height=6, units="in", res=320)
	pie(pct, labels = lbls, cex = 0.7, radius = 1.05, col=rainbow(length(lbls)))
	dev.off()

}

########## Physical properties (GC Content and Melting temperature) ##########

phys_prop_plot <- function(input, output_dir) {

	png(filename = file.path(output_dir, "gc_plot.png"), width=6, height=6, units="in", res=320)
	gc <- hist(input$GCcontent, col = "blue", breaks = seq(20,60, by = 2), freq = TRUE, xlab ="GC Content of Probes (%)", main = NULL)
	dev.off()

	png(filename = file.path(output_dir, "tm_plot.png"), width=6, height=6, units="in", res=320)
	tm <- hist(input$`Melting Temperature (C)`, col = "red", breaks = seq(50,110, by = 2), freq = TRUE, xlab = "Melting Temperature (Tm) of Probes (\u00B0C)", main = NULL)
	dev.off()

	avg_gc <- mean(input$GCcontent)
	avg_tm <- mean(input$`Melting Temperature (C)`)

	cat(paste("Average GC Content (%) = ", avg_gc, "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
	cat(paste("Average Melting Temperature (\u00B0C) = ", avg_tm, "\n\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)

}


########## Distribution of coverage per virus family ##########

coverage_plot <- function(input, output_dir, syndrome, remove_zero, non_zero) {

### Reporting dataframe (not including # baits per virus family)
	
	headers <- c("Virus Family", "# References", "Average Depth of Coverage (fold)", "Median Depth of Coverage (fold)", "Minimum Depth of Coverage (fold)", "Maximum Depth of Coverage (fold)")
	report <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, headers)))

### Determine which dataset to analyze (all or non-zero)

	if (isTRUE(non_zero)) {
		chunk_avg <- subset(input$`Average Depth of Coverage Across Non-Zero`, (input$`Bait Origin` %in% syndrome))
	} else {
		chunk_avg <- subset(input$`Average Depth of Coverage`, (input$`Bait Origin` %in% syndrome))
	}

### Calculate average depth of coverage

	chunk_avg[!is.finite(as.numeric(sub("x", "", chunk_avg)))] <- 0
	avg_depth <- mean(as.numeric(sub("x", "", chunk_avg)))

	if (isTRUE(non_zero)) {
		cat(paste("Average fold coverage (non-zero regions only) = ", round(avg_depth, digits = 2), "\n\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
	} else {
		cat(paste("Average fold coverage = ", round(avg_depth, digits = 2), "\n\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
	}

### En masse coverage plot (all virus families in one plot)

	num_syn <- round_ceiling(sqrt(length(syndrome)),1)
	if (sqrt(length(syndrome))%%1 == 0) {
		num_syn_col <- num_syn
	} else if (num_syn > 1 && num_syn*(num_syn-1) >= length(syndrome)) {
		num_syn_col <- num_syn-1
	} else {
		num_syn_col <- num_syn
	}
	#png(filename = file.path(output_dir, "coverage_plots.png"), width=round_ceiling(sqrt(num_syn+round_ceiling(num_syn*0.75,1)-1),1), height=round_ceiling(sqrt(num_syn)+round_ceiling(num_syn*0.25,1)-1,1), units="in", res=320)
	png(filename = file.path(output_dir, "coverage_plots.png"), width=num_syn+round_ceiling(num_syn*0.25,1), height=num_syn+round_ceiling(num_syn*0.75,1), units="in", res=320)
	par(mfrow=c(num_syn,num_syn_col))

	for (family in syndrome) {

		if (isTRUE(non_zero)) {
			chunk <- subset(input$`Average Depth of Coverage Across Non-Zero`, (input$`Bait Origin` == family))
			chunk[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		} else {
			chunk <- subset(input$`Average Depth of Coverage`, (input$`Bait Origin` == family))
			chunk[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		}

		fold <- as.numeric(sub("x", "", chunk))
		if (length(fold) == 0) {
			fold = 0 # no coverage (and no values found)
		} else {
			fold[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		}
		lbls <- paste(fold, "x", sep="")
		
		max_family <- max(fold)
		if (is.finite(round_ceiling(max_family, 1)) && round_ceiling(max_family, 1) > 0) {
			max_break <- round_ceiling(max_family, 1) + 1
		} else { 
			max_break <- 10
		}

### PLOT (PNG)
		hist(x = fold, breaks = seq(0,max_break,0.5), col = "black", xlab = "Average Depth of Coverage", main = NULL, cex.lab = 0.7, cex.axis = 0.8) 
		title(family, cex.main = 0.9)
		
	}
	dev.off()

### individual virus family plots

	# Reporting parameters
	none_covered <- 0
	ref_count <- 0
	
	for (family in syndrome) {
		
		png(filename = file.path(output_dir, "coverage",  paste(family, ".png", sep="")), width = 8, height = 6, units = "in", res=320)

		if (isTRUE(non_zero)) {
			chunk <- subset(input$`Average Depth of Coverage Across Non-Zero`, (input$`Bait Origin` == family))
			chunk[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		} else {
			chunk <- subset(input$`Average Depth of Coverage`, (input$`Bait Origin` == family))
			chunk[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		}

		fold <- as.numeric(sub("x", "", chunk))
		if (length(fold) == 0) {
			fold = 0 # no coverage (and no values found)
		} else {
			fold[!is.finite(as.numeric(sub("x", "", chunk)))] <- 0
		}
		lbls <- paste(fold, "x", sep="")

### REPORTING
		len_family <- length(fold)
		cat(paste("Number of", family, "references: ", len_family, "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		avg_family <- mean(fold)
		cat(paste("Average fold coverage for ", family, "is ", round(avg_family, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		median_family <- median(fold)
		cat(paste("Median fold coverage for ", family, "is ", round(median_family, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		min_family <- min(fold)
		cat(paste("Minimum fold coverage for ", family, "is ", round(min_family, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		max_family <- max(fold)
		cat(paste("Maximum fold coverage for ", family, "is ", round(max_family, digits = 2), "\n-----\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		if (is.finite(round_ceiling(max_family, 1)) && round_ceiling(max_family, 1) > 0) {
			max_break <- round_ceiling(max_family, 1) + 1
		} else { 
			max_break <- 10
		}
		
		none_covered <- none_covered + length(subset(fold, (fold == 0)))
		ref_count <- ref_count + length(fold)
		
### Update reporter dataframe
		if (isTRUE(remove_zero)) {
			if (max_family > 0) {
				vf_update <- data.frame(family,len_family, round(avg_family, digits = 2), round(median_family, digits = 2), round(min_family, digits = 2), round(max_family, digits = 2))
				names(vf_update) <- headers
				report <- rbind(report, vf_update)
			}
		} else {
			vf_update <- data.frame(family,len_family, round(avg_family, digits = 2), round(median_family, digits = 2), round(min_family, digits = 2), round(max_family, digits = 2))
			names(vf_update) <- headers
			report <- rbind(report, vf_update)
		}

### PLOT (PNG)
		hist(x = fold, breaks = seq(0,max_break, 0.5), col = "black", xlab = "Average Depth of Coverage", main = NULL, cex.lab = 0.7, cex.axis = 0.8)
		title(family, cex.main = 0.9)
		dev.off()
		
	}
		
### REPORTING AND OUTPUT
		cat(paste("\nTotal number of virus genomes with 0 coverage: ", none_covered, "out of", ref_count, "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		cat("\n-----\n\n", file = file.path(output_dir, "summary.txt"), append = TRUE)
		write.table(report, file = file.path(output_dir, 'genome_coverage.txt'), quote = FALSE, sep = '\t', row.names = FALSE)
	

}

completeness_plot <- function(input, output_dir, syndrome, remove_zero) {

### Reporting dataframe (without # bait per virus family)
	
	headers <- c("Virus Family", "# References", "Average Genome Completeness (%)", "Median Genome Completeness (%)", "Minimum Genome Completeness (%)", "Maximum Genome Completeness (%)")
	report <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, headers)))

### Calculate average genome completeness
	
	chunk_cov <- subset(input$`Genome Completeness`, (input$`Bait Origin` %in% syndrome))
	chunk_cov[!is.finite(as.numeric(sub("%", "", chunk_cov)))] <- 0
	avg_cov <- mean(as.numeric(sub("%", "", chunk_cov)))

	cat(paste("Average % genome completeness = ", round(avg_cov, digits = 2), "\n\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
	
### En masse genome completeness

	num_syn <- round_ceiling(sqrt(length(syndrome)),1)
	if (sqrt(length(syndrome))%%1 == 0) {
		num_syn_col <- num_syn
	} else if (num_syn > 1 && (num_syn*(num_syn-1)) >= length(syndrome)) {
		num_syn_col <- num_syn-1
	} else {
		num_syn_col <- num_syn
	}
	#png(filename = file.path(output_dir, "completeness_plots.png"), width=round_ceiling(sqrt(num_syn+14),1), height=round_ceiling(sqrt(num_syn),1)+4, units="in", res=320)
	png(filename = file.path(output_dir, "completeness_plots.png"), width=num_syn+round_ceiling(num_syn*0.25,1), height=num_syn+round_ceiling(num_syn*0.75,1), units="in", res=320)
	par(mfrow=c(num_syn,num_syn_col))

	for (family in syndrome) {
		
		chunk <- subset(input$`Genome Completeness`, (input$`Bait Origin` == family))
		chunk[!is.finite(as.numeric(sub("%", "", chunk)))] <- 0

		complete <- as.numeric(sub("%", "", chunk))
		if (length(complete) == 0) {
			complete = 0 # no coverage (and no values found)
		} else {
			complete[!is.finite(as.numeric(sub("%", "", chunk)))] <- 0
		}
		lbls <- paste(complete, "%", sep="")
		
	
### PLOT (PNG)
		hist(x = complete, breaks = seq(0,100,10), xlim = c(0,100), col = "black", xlab = "Genome completeness (%)", main = NULL, cex.lab = 0.7, cex.axis = 0.8) 
		title(family, cex.main = 0.9)
	
	}
	dev.off()

### Individual plots 

	# Reporting parameters
	none_breadth <- 0
	ref_count_cov <- 0
	
	for (family in syndrome) {
		png(filename = file.path(output_dir, "completeness",  paste(family, ".png", sep="")), width = 8, height = 6, units = "in", res=320)
	
		chunk <- subset(input$`Genome Completeness`, (input$`Bait Origin` == family))
		chunk[!is.finite(as.numeric(sub("%", "", chunk)))] <- 0

		complete <- as.numeric(sub("%", "", chunk))
		if (length(complete) == 0) {
			complete = 0 # no coverage (and no values found)
		} else {
			complete[!is.finite(as.numeric(sub("%", "", chunk)))] <- 0
		}
		lbls <- paste(complete, "%", sep="")

### REPORTING
		len_complete <- length(complete)
		cat(paste("Number of", family, "references: ", len_complete, "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		avg_complete <- mean(complete)
		cat(paste("Average % genome completeness for ", family, "is ", round(avg_complete, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		median_complete <- median(complete)
		cat(paste("Median % genome completeness for ", family, "is ", round(median_complete, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		min_complete <- min(complete)
		cat(paste("Minimum % genome completeness for ", family, "is ", round(min_complete, digits = 2), "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		max_complete <- max(complete)
		cat(paste("Maximum % genome completeness for ", family, "is ", round(max_complete, digits = 2), "\n-----\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		
		none_breadth <- none_breadth + length(subset(complete, (complete == 0)))
		ref_count_cov <- ref_count_cov + length(complete)
		
### Update reporter dataframe
		if (isTRUE(remove_zero)) { 
			if (max_complete > 0) { 
				vf_update <- data.frame(family,len_complete, round(avg_complete, digits = 2), round(median_complete, digits = 2), round(min_complete, digits = 2), round(max_complete, digits = 2))
				names(vf_update) <- headers
				report <- rbind(report, vf_update)
			}
		} else {
			vf_update <- data.frame(family,len_complete, round(avg_complete, digits = 2), round(median_complete, digits = 2), round(min_complete, digits = 2), round(max_complete, digits = 2))
			names(vf_update) <- headers
			report <- rbind(report, vf_update)
		}

### PLOT (PNG)
		hist(x = complete, breaks = seq(0,100,10), xlim = c(0,100), col = "black", xlab = "Genome completeness (%)", main = NULL, cex.lab = 0.7, cex.axis = 0.8) 
		title(family, cex.main = 0.9)
		dev.off()
		
	}

### REPORTING AND OUTPUT
		cat(paste("\nTotal number of virus genomes with 0% genome completeness: ", none_breadth, "out of", ref_count_cov, "\n"), file = file.path(output_dir, "summary.txt"), append = TRUE)
		cat("\n-----\n\n", file = file.path(output_dir, "summary.txt"), append = TRUE)
		write.table(report, file = file.path(output_dir, 'genome_completeness.txt'), quote = FALSE, sep = '\t', row.names = FALSE)

}

########## PARSER, OUTPUT, AND FUNCTIONS ##########

parser <- ArgumentParser(description='Plot out bait coverage in a variety of forms')
parser$add_argument("--output", "-o", required=TRUE, help="Output directory name")
parser$add_argument("--phys-prop", "-p", default=NULL, help="File containing GC content and melting temperature (in that order) for a bait set in TSV format")
parser$add_argument("--proportion", "-f", default=NULL, help="File containing the number (and percentage) of baits per virus family in TSV format (family, number of baits, percentage of baits)")
parser$add_argument("--baltimore", "-b", default=NULL, help="File containing the number (and percentage) of baits per balitmore classifier (ex. +ssRNA, dsDNA, etc.) in TSV format (classifier, number of baits, percentage of baits)")
parser$add_argument("--coverage", "-c", default=NULL, help="File containing coverage information in TSV format (accession, genome completeness, average depth of coverage)")
parser$add_argument("--syndrome", "-s", default="all", help=paste("Examine coverage for specific groups of virus families (or syndromes). Input can be a comma-separated list of virus families or one of the accepted entries above", cat(usage(), sep = "\n")))
parser$add_argument("--remove-zero", action="store_true", help="Remove 0 coverage information from plots and tables (only applies if '--coverage' is provided)")
parser$add_argument("--non-zero", action="store_true", help="Assess (depth of) coverage across non-zero regions only. Useful for gapped coverages where breadth of coverage (completeness) is expected to be low")
args <- parser$parse_args()


### Output directory 
output_dir <- args$output 
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
	dir.create(file.path(output_dir, "coverage"))
	dir.create(file.path(output_dir, "completeness"))
} else {
	print("Output directory exists! Please rename and try again!")
	quit("no", 1, FALSE)
}

### Determine syndrome analysis (argument is string, match with preloaded list)
syndrome <- tolower(args$syndrome)
if (syndrome %in% c("gastro", "gastrointestinal")) {
	syndrome <- gastro
} else if (syndrome %in% c("blood")) {
	syndrome <- blood
} else if (syndrome %in% c("respir", "respiratory")) {
	syndrome <- respir
} else if (syndrome %in% c("human")) {
	syndrome <- human
} else if (syndrome %in% c("csf")) {
	syndrome <- csf
} else if (syndrome %in% c("urine")) {
	syndrome <- urine
} else if (syndrome %in% c("sjf")) {
	syndrome <- sjf
} else if (syndrome %in% c("special")) {
	syndrome <- special
} else if (syndrome %in% c("all")) {
	syndrome <- all
} else {
	syndrome <- lapply(as.list(strsplit(gsub(" ", "", syndrome, fixed = TRUE), ",")[[1]]), caps)
}

### Create summary file with all text output
cat("Summary of statistics:\n\n", file = file.path(output_dir, "summary.txt"), append = TRUE)

### All input files
if (!is.null(args$proportion)) { 
	if (file.exists(args$proportion)) {
		data_families <- as.data.frame(fread(args$proportion, header = FALSE))
		family_proportion_plot(data_families, output_dir)
	} else {
		print(paste("--proportion", args$proportion, "does not exist! Skipping"))
	}
}
if (!is.null(args$baltimore)) {
	if (file.exists(args$baltimore)) {
		data_balt <- as.data.frame(fread(args$baltimore, header = FALSE))
		baltimore_proportion_plot(data_balt, output_dir)
	} else {
		print(paste("--baltimore", args$baltimore, "does not exist! Skipping"))
	}
}

if (!is.null(args$phys_prop)) {
	if (file.exists(args$phys_prop)) {
		data_phys <- as.data.frame(fread(args$phys_prop, header = TRUE))
		phys_prop_plot(data_phys, output_dir)
	} else {
		print(paste("--phys-prop", args$phys_prop, "does not exist! Skipping"))
	}
}

if (!is.null(args$coverage)) {
	if (file.exists(args$coverage)) {
		data_coverage <- as.data.frame(fread(args$coverage, header = TRUE))
		coverage_plot(data_coverage, output_dir, syndrome, args$remove_zero, args$non_zero)
		
		### TODO: if non_zero, perhaps don't run completeness on total dataset?
		
		completeness_plot(data_coverage, output_dir, syndrome, args$remove_zero)
	} else {
		print(paste("--coverage", args$coverage, "does not exist! Skipping"))
	}

}



