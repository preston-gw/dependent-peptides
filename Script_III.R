# +-----------------------------------------------------------------------------------+
# | Supporting information file S3 Script                                             |
# +-----------------------------------------------------------------------------------+
# | Script III                                                                        |
# |                                                                                   |
# | To accompany the report 'Visualisation tools for dependent peptide searches to    |
# | support the exploration of in vitro protein modifications' by G. W. Preston, L.   |
# | Yang, D. H. Phillips and C. S. Maier                                              |
# +-----------------------------------------------------------------------------------+

#   Usage notes:
#   1. This script filters the dependent peptides from a given analysis (e.g. an 
#      analysis of NEM-treated BSA), enriches them on the basis of 'constant 
#      conjunction' (see note 2) and generates a dependent-peptide localisation plot 
#      and a mass-shift frequency histogram.
#   2. In order to be visualised, a dependent peptide must be 'constantly conjoined'
#      with the experimental condition (e.g. NEM treatment); that is, it must:
#      - have been detected in the given analysis;
#      - have been detected in a second analysis of the same sample; 
#      - have gone undetected in both of two analyses of a second sample (e.g. a 
#        sample of untreated BSA).
#   3. The script requires:
#      i.   R version 3.6.0 or later; 
#      ii.  Windows (a Windows-specific graphics device will be used);
#      iii. package 'seqinR';
#      iv.  four allPeptides.txt files (e.g. two for NEM-treated BSA and two for
#           untreated BSA);
#      v.   a *.fasta file containing the sequence of a protein of interest;
#      vi.  the ID (e.g. UniProt ID) of the protein of interest, exactly as it 
#           appears in the allPeptides.txt files.
#   4. Data from our BSA study, including the allPeptides.txt and *.fasta files 
#      specified in the script, are available via PRIDE/ProteomeXchange (accession 
#      PXD013040).
#   5. Further instructions are provided within the script. Comments beginning with a
#      greater-than symbol (>) are instructions for users.
#   6. Detailed explanatory notes can be found in Tables B and C of supporting 
#      information file S1 Text. 

# +-----------------------------------------------------------------------------------+

# Choose options
options("stringsAsFactors" = FALSE)

# Load seqinR
library(seqinr)

# Get protein sequence
# > Specify the path to a *.fasta file
protein.sequence <- read.fasta(paste("C:", "\\", "Users", "\\", "George", "\\",
	"MaxQuant_1.6.0.1", "\\", "Databases", "\\", "4f5s_A.fasta", sep = ""),
	seqtype = "AA",
	seqonly = TRUE)[[1]]

# > Specify protein of interest
protein.ID <- "4F5S"

# > Specify tolerance
tolerance <- 0.01

# Get search results (treated protein, analysis 1)
# > Specify the path to an allPeptides.txt file
treated_replicate1_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", 
	"George", "\\",	"Documents", "\\", "Model data", "\\", 
	"treated_replicate1_allPeptides.txt", sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Get search results (treated protein, analysis 2)
# > Specify the path to an allPeptides.txt file
treated_replicate2_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", 
	"George", "\\", "Documents", "\\", "Model data", "\\", 
	"treated_replicate2_allPeptides.txt", sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Get search results (untreated protein, analysis 1)
# > Specify the path to an allPeptides.txt file
control_replicate1_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", 
	"George", "\\", "Documents", "\\", "Model data", "\\", 
	"control_replicate1_allPeptides.txt", sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Get search results (untreated protein, analysis 2)
# > Specify the path to an allPeptides.txt file
control_replicate2_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", 
	"George", "\\", "Documents", "\\", "Model data", "\\", 
	"control_replicate2_allPeptides.txt", sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Filter search results (treated protein, analysis 1)
treated_replicate1_dependent.peptides <- treated_replicate1_all.peptides[
	!is.na(treated_replicate1_all.peptides$DP.Mass.Difference) &
	treated_replicate1_all.peptides$DP.Mass.Difference < 500.5 &
	treated_replicate1_all.peptides$DP.Mass.Difference >= -500.5 &
	treated_replicate1_all.peptides$DP.Score > 80 &	
	treated_replicate1_all.peptides$DP.PEP < 0.01 &
	treated_replicate1_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, treated_replicate1_all.peptides$DP.Proteins) &
	!(treated_replicate1_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(treated_replicate1_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	treated_replicate1_all.peptides$DP.AA == "C") &
	!grepl("Cation:", treated_replicate1_all.peptides$DP.Modification) &
	treated_replicate1_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Filter search results (treated protein, analysis 2)
treated_replicate2_dependent.peptides <- treated_replicate2_all.peptides[
	!is.na(treated_replicate2_all.peptides$DP.Mass.Difference) &
	treated_replicate2_all.peptides$DP.Mass.Difference < 500.5 &
	treated_replicate2_all.peptides$DP.Mass.Difference >= -500.5 &
	treated_replicate2_all.peptides$DP.Score > 80 &	
	treated_replicate2_all.peptides$DP.PEP < 0.01 &
	treated_replicate2_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, treated_replicate2_all.peptides$DP.Proteins) &
	!(treated_replicate2_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(treated_replicate2_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	treated_replicate2_all.peptides$DP.AA == "C") &
	!grepl("Cation:", treated_replicate2_all.peptides$DP.Modification) &
	treated_replicate2_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Filter search results (untreated protein, analysis 1)
control_replicate1_dependent.peptides <- control_replicate1_all.peptides[
	!is.na(control_replicate1_all.peptides$DP.Mass.Difference) &
	control_replicate1_all.peptides$DP.Mass.Difference < 500.5 &
	control_replicate1_all.peptides$DP.Mass.Difference >= -500.5 &
	control_replicate1_all.peptides$DP.Score > 80 &	
	control_replicate1_all.peptides$DP.PEP < 0.01 &
	control_replicate1_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, control_replicate1_all.peptides$DP.Proteins) &
	!(control_replicate1_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(control_replicate1_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	control_replicate1_all.peptides$DP.AA == "C") &
	!grepl("Cation:", control_replicate1_all.peptides$DP.Modification) &
	control_replicate1_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Filter search results (untreated protein, analysis 2)
control_replicate2_dependent.peptides <- control_replicate2_all.peptides[
	!is.na(control_replicate2_all.peptides$DP.Mass.Difference) &
	control_replicate2_all.peptides$DP.Mass.Difference < 500.5 &
	control_replicate2_all.peptides$DP.Mass.Difference >= -500.5 &
	control_replicate2_all.peptides$DP.Score > 80 &	
	control_replicate2_all.peptides$DP.PEP < 0.01 &
	control_replicate2_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, control_replicate2_all.peptides$DP.Proteins) &
	!(control_replicate2_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(control_replicate2_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	control_replicate2_all.peptides$DP.AA == "C") &
	!grepl("Cation:", control_replicate2_all.peptides$DP.Modification) &
	control_replicate2_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Make a non-redundant list of peptide sequences (untreated protein, analysis 1)
control_replicate1_peptide.sequences <- unique(c(
		control_replicate1_all.peptides$Sequence[
		control_replicate1_all.peptides$Score > 80 &
		grepl(protein.ID, control_replicate1_all.peptides$Proteins)],
	control_replicate1_all.peptides$DP.Base.Sequence[
		control_replicate1_all.peptides$DP.Score > 80 &
		control_replicate1_all.peptides$DP.PEP < 0.01 &
		control_replicate1_all.peptides$DP.Decoy != "+" &
		grepl(protein.ID, control_replicate1_all.peptides$DP.Proteins)]))

# Make a non-redundant list of peptide sequences (untreated protein, analysis 2)
control_replicate2_peptide.sequences <- unique(c(
	control_replicate2_all.peptides$Sequence[
		control_replicate2_all.peptides$Score > 80 &
		grepl(protein.ID, control_replicate2_all.peptides$Proteins)],
	control_replicate2_all.peptides$DP.Base.Sequence[
		control_replicate2_all.peptides$DP.Score > 80 &
		control_replicate2_all.peptides$DP.PEP < 0.01 &
		control_replicate2_all.peptides$DP.Decoy != "+" &
		grepl(protein.ID, control_replicate2_all.peptides$DP.Proteins)]))

# Compile sequences that are common to both lists
control_common.peptide.sequences <- control_replicate1_peptide.sequences[
	control_replicate1_peptide.sequences %in% 
	control_replicate2_peptide.sequences]

# Make empty data frame for sequences, mass shifts, retention times and test outcomes
treated_DP.constancy.results <- data.frame(
	Sequence = character(nrow(treated_replicate1_dependent.peptides)),
	Retention.time.1 = numeric(nrow(treated_replicate1_dependent.peptides)),
	Mass.difference.1 = numeric(nrow(treated_replicate1_dependent.peptides)),
	Retention.time.2 = numeric(nrow(treated_replicate1_dependent.peptides)),
	Mass.difference.2 = numeric(nrow(treated_replicate1_dependent.peptides)),
	Mean.mass.difference = numeric(nrow(treated_replicate1_dependent.peptides)),
	In.control_replicate1 = character(nrow(treated_replicate1_dependent.peptides)),
	In.control_replicate2 = character(nrow(treated_replicate1_dependent.peptides)))
for(i in 1:nrow(treated_DP.constancy.results))
	{for(j in 1:ncol(treated_DP.constancy.results))
	treated_DP.constancy.results[i, j] <- NA}

# Perform first round of pairwise comparisons (treated <-> treated), deposit data
for(i in 1:nrow(treated_replicate1_dependent.peptides))
	{for(j in 1:nrow(treated_replicate2_dependent.peptides))
	if(treated_replicate1_dependent.peptides$DP.Base.Sequence[i] ==
		treated_replicate2_dependent.peptides$DP.Base.Sequence[j] &
		sign(treated_replicate1_dependent.peptides$DP.Mass.Difference[i]) ==
		sign(treated_replicate2_dependent.peptides$DP.Mass.Difference[j]) &
		treated_replicate1_dependent.peptides$DP.Mass.Difference[i] + 
			tolerance >
		treated_replicate2_dependent.peptides$DP.Mass.Difference[j] - 
			tolerance &
		treated_replicate1_dependent.peptides$DP.Mass.Difference[i] - 
			tolerance <
		treated_replicate2_dependent.peptides$DP.Mass.Difference[j] + 
			tolerance)
	{treated_DP.constancy.results$Sequence[i] <- 
		treated_replicate1_dependent.peptides$DP.Base.Sequence[i]
	treated_DP.constancy.results$Retention.time.1[i] <- 
		treated_replicate1_dependent.peptides$Retention.time[i]
	treated_DP.constancy.results$Mass.difference.1[i] <- 
		treated_replicate1_dependent.peptides$DP.Mass.Difference[i]
	treated_DP.constancy.results$Retention.time.2[i] <- 
		treated_replicate2_dependent.peptides$Retention.time[j]
	treated_DP.constancy.results$Mass.difference.2[i] <- 
		treated_replicate2_dependent.peptides$DP.Mass.Difference[j]
	treated_DP.constancy.results$Mean.mass.difference[i] <- 
		mean(c(treated_replicate1_dependent.peptides$DP.Mass.Difference[i],
		treated_replicate2_dependent.peptides$DP.Mass.Difference[j]))}}

# Remove rows of NAs
treated_constant.DPs <- treated_DP.constancy.results[
	!is.na(treated_DP.constancy.results$Sequence), ]

# Perform second round of comparisons (treated <-> untreated), deposit outcomes
if(nrow(control_replicate1_dependent.peptides) > 0)
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 1:nrow(control_replicate1_dependent.peptides))
	if(!is.na(control_replicate1_dependent.peptides$DP.Base.Sequence[j]) &
		treated_constant.DPs$Sequence[i] ==
		control_replicate1_dependent.peptides$DP.Base.Sequence[j] &
		sign(treated_constant.DPs$Mean.mass.difference[i]) ==
		sign(control_replicate1_dependent.peptides$DP.Mass.Difference[j]) &
		treated_constant.DPs$Mean.mass.difference[i] + 
			tolerance >
		control_replicate1_dependent.peptides$DP.Mass.Difference[j] - 
			tolerance &
		treated_constant.DPs$Mean.mass.difference[i] - 
			tolerance <
		control_replicate1_dependent.peptides$DP.Mass.Difference[j] + 
			tolerance)
	treated_constant.DPs$In.control_replicate1[i] <- "Yes"}

# Perform third round of comparisons (treated <-> untreated), deposit outcomes
if(nrow(control_replicate2_dependent.peptides) > 0)
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 1:nrow(control_replicate2_dependent.peptides))
	if(!is.na(control_replicate2_dependent.peptides$DP.Base.Sequence[j]) &
		treated_constant.DPs$Sequence[i] ==
		control_replicate2_dependent.peptides$DP.Base.Sequence[j] &
		sign(treated_constant.DPs$Mean.mass.difference[i]) ==
		sign(control_replicate2_dependent.peptides$DP.Mass.Difference[j]) &
		treated_constant.DPs$Mean.mass.difference[i] + 
			tolerance >
		control_replicate2_dependent.peptides$DP.Mass.Difference[j] - 
			tolerance &
		treated_constant.DPs$Mean.mass.difference[i] - 
			tolerance <
		control_replicate2_dependent.peptides$DP.Mass.Difference[j] + 
			tolerance)
	treated_constant.DPs$In.control_replicate2[i] <- "Yes"}

# Add negative outcomes
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 7:8)
	if(is.na(treated_constant.DPs[i, j]))
	treated_constant.DPs[i, j] <- "No"}

# Isolate 'constantly conjoined' dependent peptides
treated_constantly.conjoined.DPs <- treated_constant.DPs[
	treated_constant.DPs$In.control_replicate1 == "No" &
	treated_constant.DPs$In.control_replicate2 == "No", ]

# Choose a graphics device and specify plotting-area dimensions
windows(7, 3)

# Define custom plot margins, layout and common axis limits
par(mar = c(5, 4, 1, 2) + 0.1)
layout(matrix(c(1, 1, 1, 1, 2, 2, 2), 
	nrow = 1, 
	ncol = 7), 
	respect = FALSE)
delta.mass.axis.limits <- c(-500, 500)

# Create an empty left-hand plot
plot(1,	0, 
	xlim = c(0, nchar(protein.sequence)),	
	ylim = delta.mass.axis.limits,
	type = "n",
	xlab = "Amino acid residue",
	ylab = expression(paste(Delta, "m / Da", sep = "")),
	las = 1)

# Add a dashed 'backbone'
lines(c(0.5, nchar(protein.sequence) + 0.5), 
	c(0, 0),
	lty = 3, 
	col = "red")

# Solidify observed regions of the backbone
for(i in 1:length(control_common.peptide.sequences))
	{for(j in 1:(nchar(protein.sequence) -
		(nchar(control_common.peptide.sequences[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(control_common.peptide.sequences[i]) - 1))
	if(segment == control_common.peptide.sequences[i])
	lines(c(j - 0.5,
		j + (nchar(control_common.peptide.sequences[i]) - 1) + 0.5),
		c(0, 0), 
		col = "red")}}

# Add dependent peptides
for(i in 1:nrow(treated_constantly.conjoined.DPs))
	{for(j in 1:(nchar(protein.sequence) - 
		(nchar(treated_constantly.conjoined.DPs$Sequence[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(treated_constantly.conjoined.DPs$Sequence[i]) - 1))
	if(segment == treated_constantly.conjoined.DPs$Sequence[i])
	{vertices <- matrix(data = c(j - 0.5,
		j + (nchar(treated_constantly.conjoined.DPs$Sequence[i]) - 1) + 0.5,
		j + (nchar(treated_constantly.conjoined.DPs$Sequence[i]) - 1) + 0.5,
		j - 0.5,
		0,
		0, 
		treated_constantly.conjoined.DPs$Mean.mass.difference[i], 
		treated_constantly.conjoined.DPs$Mean.mass.difference[i]),
		nrow = 4,
		ncol = 2)
	polygon(vertices, 
		border = gray(level = 0, alpha = 0.5))}}}

# Prepare histogram
histogram <- hist(treated_constantly.conjoined.DPs$Mean.mass.difference,
	breaks = seq(-500.5, 500.5, 1),
	right = FALSE,
	plot = FALSE)

# Draw histogram
plot(histogram, 
	main = NULL,
	xlim = delta.mass.axis.limits,
	ylim = c(0, max(histogram$counts) * 1.05),
	xlab = expression(paste(Delta, "m / Da", sep = "")),
	ylab = "Number of peptides",
	las = 1)
box()

# Specify how the histogram should be annotated
# > You could vary the annotation count
# > Type 'histogram.slices' to view the table
histogram.slices <- data.frame(Slice.number = 1:13, 
	First.cell = c(1, 101, 201, 301, 401, 451, 501, 
		502, 552, 602, 702, 802, 902),
	Last.cell = c(100, 200, 300, 400, 451, 500, 501, 
		551, 601, 701, 801, 901, 1001),
	Mass.difference.range = numeric(13),
	Annotation.count = c(0, 0, 0, 1, 3, 3, 1, 3, 3, 2, 1, 0, 0))

for(i in 1:nrow(histogram.slices))
histogram.slices$Mass.difference.range[i] <- paste(
	as.character(histogram$mids[histogram.slices$First.cell[i]]),
	"to",
	as.character(histogram$mids[histogram.slices$Last.cell[i]]))

# Annotate the histogram
for(i in 1:nrow(histogram.slices))
{slice_counts <- histogram$counts[histogram.slices$First.cell[i]:
	histogram.slices$Last.cell[i]]
if(sum(slice_counts) > 0 & histogram.slices$Annotation.count[i] > 0)
{cells.to.annotate <- histogram.slices$First.cell[i] + 
	order(-slice_counts)[1:histogram.slices$Annotation.count[i]] - 1
text(histogram$mids[cells.to.annotate], 
	histogram$counts[cells.to.annotate] + max(histogram$counts) / 19, 
	labels = as.character(histogram$mids[cells.to.annotate]))}}

# +-----------------------------------------------------------------------------------+