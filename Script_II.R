# +-----------------------------------------------------------------------------------+
# | Supporting information file S2 Script                                             |
# +-----------------------------------------------------------------------------------+
# | Script II                                                                         |
# |                                                                                   |
# | To accompany the report 'Visualisation tools for dependent peptide searches to    |
# | support the exploration of in vitro protein modifications' by G. W. Preston, L.   |
# | Yang, D. H. Phillips and C. S. Maier                                              |
# +-----------------------------------------------------------------------------------+

#   Usage notes:
#   1. This script filters the dependent peptides from a given analysis (e.g. an 
#      analysis of NEM-treated BSA), enriches them for modifications of interest (see
#      note 2) and generates a dependent-peptide localisation plot and a mass-shift 
#      frequency histogram.
#   2. In order to be visualised, a dependent peptide must: 
#      - have been detected in the given analysis;
#      - have gone undetected in a second analysis (e.g. an analysis of untreated 
#        BSA).
#   3. The script requires:
#      i.   R version 3.6.0 or later;
#      ii.  Windows (a Windows-specific graphics device will be used);
#      iii. package 'seqinR';
#      iv.  two allPeptides.txt files (e.g. one for NEM-treated BSA and one for
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

# Get search results (treated protein)
# > Specify the path to an allPeptides.txt file
treated_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", "George", "\\",
	"Documents", "\\", "Model data", "\\", "treated_replicate1_allPeptides.txt", 
	sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Get search results (untreated protein)
# > Specify the path to an allPeptides.txt file
control_all.peptides <- read.table(paste("C:", "\\", "Users", "\\", "George", "\\",
	"Documents", "\\", "Model data", "\\", "control_replicate1_allPeptides.txt", 
	sep = ""),
	sep = "\t",
	header = TRUE,
	fill = TRUE)

# Make a non-redundant list of observed peptide sequences
peptide.sequences <- unique(c(control_all.peptides$Sequence[
	control_all.peptides$Score > 80 &
		grepl(protein.ID, control_all.peptides$Proteins)],
	control_all.peptides$DP.Base.Sequence[control_all.peptides$DP.Score > 80 &	
		control_all.peptides$DP.PEP < 0.01 &
		control_all.peptides$DP.Decoy != "+" &
		grepl(protein.ID, control_all.peptides$DP.Proteins)]))

# Filter search results (treated protein)
treated_dependent.peptides <- treated_all.peptides[
	!is.na(treated_all.peptides$DP.Mass.Difference) &
	treated_all.peptides$DP.Mass.Difference < 500.5 &
	treated_all.peptides$DP.Mass.Difference >= -500.5 &
	treated_all.peptides$DP.Score > 80 &	
	treated_all.peptides$DP.PEP < 0.01 &
	treated_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, treated_all.peptides$DP.Proteins) &
	!(treated_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(treated_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	treated_all.peptides$DP.AA == "C") &
	!grepl("Cation:", treated_all.peptides$DP.Modification) &
	treated_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Filter search results (untreated protein)
control_dependent.peptides <- control_all.peptides[
	!is.na(control_all.peptides$DP.Mass.Difference) &
	control_all.peptides$DP.Mass.Difference < 500.5 &
	control_all.peptides$DP.Mass.Difference >= -500.5 &
	control_all.peptides$DP.Score > 80 &	
	control_all.peptides$DP.PEP < 0.01 &
	control_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, control_all.peptides$DP.Proteins) &
	!(control_all.peptides$DP.Modification %in% c("Carbamidomethyl",
		"Loss of ammonia",
		"Loss of water",
		"Unmodified")) &
	!(control_all.peptides$DP.Modification == "CarbamidomethylDTT" &
	control_all.peptides$DP.AA == "C") &
	!grepl("Cation:", control_all.peptides$DP.Modification) &
	control_all.peptides$DP.Peptide.Length.Difference >= 0, ]

# Make a vector that identifies which records to subtract
records.to.subtract <- as.numeric(rep(NA, nrow(treated_dependent.peptides)))
for(i in 1:nrow(treated_dependent.peptides))
{for(j in 1:nrow(control_dependent.peptides))
if(treated_dependent.peptides$DP.Base.Sequence[i] ==
		control_dependent.peptides$DP.Base.Sequence[j] &
	sign(treated_dependent.peptides$DP.Mass.Difference[i]) ==
		sign(control_dependent.peptides$DP.Mass.Difference[j]) &
	treated_dependent.peptides$DP.Mass.Difference[i] + tolerance >
		control_dependent.peptides$DP.Mass.Difference[j] - tolerance &
	treated_dependent.peptides$DP.Mass.Difference[i] - tolerance <
		control_dependent.peptides$DP.Mass.Difference[j] + tolerance)
records.to.subtract[i] <- i}

# Subtract the records
treated_unique.dependent.peptides <- treated_dependent.peptides[
	-records.to.subtract[!is.na(records.to.subtract)], ]

# Choose a graphics device and specify plot dimensions
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
for(i in 1:length(peptide.sequences))
	{for(j in 1:(nchar(protein.sequence) -
		(nchar(peptide.sequences[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(peptide.sequences[i]) - 1))
	if(segment == peptide.sequences[i])
	lines(c(j - 0.5,
		j + (nchar(peptide.sequences[i]) - 1) + 0.5),
		c(0, 0), 
		col = "red")}}

# Add dependent peptides
for(i in 1:nrow(treated_unique.dependent.peptides))
	{for(j in 1:(nchar(protein.sequence) - 
		(nchar(treated_unique.dependent.peptides$DP.Base.Sequence[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(
			treated_unique.dependent.peptides$DP.Base.Sequence[i]) - 1))
	if(segment == treated_unique.dependent.peptides$DP.Base.Sequence[i])
	{vertices <- matrix(data = c(j - 0.5,
		j + (nchar(
			treated_unique.dependent.peptides$DP.Base.Sequence[i]) -
				1) + 0.5,
		j + (nchar(treated_unique.dependent.peptides$DP.Base.Sequence[i]) - 
				1) + 0.5,
		j - 0.5, 0, 0, 
		treated_unique.dependent.peptides$DP.Mass.Difference[i], 
		treated_unique.dependent.peptides$DP.Mass.Difference[i]),
		nrow = 4,
		ncol = 2)
	polygon(vertices, 
		border = gray(level = 0, alpha = 0.5))}}}

# Prepare histogram
histogram <- hist(treated_unique.dependent.peptides$DP.Mass.Difference,
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
{cells.to.annotate <- histogram.slices$First.cell[i] + order(-slice_counts)[1:
	histogram.slices$Annotation.count[i]] - 1
text(histogram$mids[cells.to.annotate], 
	histogram$counts[cells.to.annotate] + max(histogram$counts) / 19, 
	labels = as.character(histogram$mids[cells.to.annotate]))}}

# +-----------------------------------------------------------------------------------+