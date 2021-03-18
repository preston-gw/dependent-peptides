# +-----------------------------------------------------------------------------------+
# | Supporting information file S4 Script                                             |
# +-----------------------------------------------------------------------------------+
# | Script IV                                                                         |
# |                                                                                   |
# | To accompany the report 'Visualisation tools for dependent peptide searches to    |
# | support the exploration of in vitro protein modifications' by G. W. Preston, L.   |
# | Yang, D. H. Phillips and C. S. Maier                                              |
# +-----------------------------------------------------------------------------------+

#   Usage notes:
#   1. This script isolates dependent peptides with a specified mass shift, enriches 
#      them on the basis of 'constant conjunction' (see note 2) and generates a 
#      probability localisation plot.
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
#           appears in the allPeptides.txt files;
#      vii. an expected mass shift and a tolerance. We use a tolerance of 0.01 Da.
#   4. Data from our BSA study, including the allPeptides.txt and *.fasta files 
#      specified in the script, are available via PRIDE/ProteomeXchange 
#      (accession PXD013040).
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

# > Input expected mass shift
expected.delta.mass <- 143.05825

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

# Make a non-redundant list of peptide sequences (untreated protein, analysis 1)
control_replicate1_peptide.sequences <- unique(c(
	control_replicate1_all.peptides$Sequence[grepl(protein.ID, 
		control_replicate1_all.peptides$Proteins) &
		control_replicate1_all.peptides$Score > 80],
	control_replicate1_all.peptides$DP.Base.Sequence[grepl(protein.ID, 
		control_replicate1_all.peptides$DP.Proteins) &
		control_replicate1_all.peptides$DP.Score > 80 &
		control_replicate1_all.peptides$DP.PEP < 0.01 &
		control_replicate1_all.peptides$DP.Decoy != "+"]))

# Make a non-redundant list of peptide sequences (untreated protein, analysis 2)
control_replicate2_peptide.sequences <- unique(c(
	control_replicate2_all.peptides$Sequence[grepl(protein.ID, 
		control_replicate2_all.peptides$Proteins) &
		control_replicate2_all.peptides$Score > 80],
	control_replicate2_all.peptides$DP.Base.Sequence[grepl(protein.ID, 
		control_replicate2_all.peptides$DP.Proteins) &
		control_replicate2_all.peptides$DP.Score > 80 &
		control_replicate2_all.peptides$DP.PEP < 0.01 &
		control_replicate2_all.peptides$DP.Decoy != "+"]))

# Compile sequences that are common to both lists
control_common.peptide.sequences <- control_replicate1_peptide.sequences[
	control_replicate1_peptide.sequences %in% 
		control_replicate2_peptide.sequences]

# Filter search results (treated protein, analysis 1)
treated_replicate1_dependent.peptides <- treated_replicate1_all.peptides[
	!is.na(treated_replicate1_all.peptides$DP.Mass.Difference) &
	treated_replicate1_all.peptides$DP.Mass.Difference < 
		expected.delta.mass + tolerance &
	treated_replicate1_all.peptides$DP.Mass.Difference > 
		expected.delta.mass - tolerance &
	treated_replicate1_all.peptides$DP.Score > 80 &	
	treated_replicate1_all.peptides$DP.PEP < 0.01 &
	treated_replicate1_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, treated_replicate1_all.peptides$DP.Proteins), ]

# Filter search results (treated protein, analysis 2)
treated_replicate2_dependent.peptides <- treated_replicate2_all.peptides[
	!is.na(treated_replicate2_all.peptides$DP.Mass.Difference) &
	treated_replicate2_all.peptides$DP.Mass.Difference < 
		expected.delta.mass + tolerance &
	treated_replicate2_all.peptides$DP.Mass.Difference > 
		expected.delta.mass - tolerance &
	treated_replicate2_all.peptides$DP.Score > 80 &	
	treated_replicate2_all.peptides$DP.PEP < 0.01 &
	treated_replicate2_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, treated_replicate2_all.peptides$DP.Proteins), ]

# Do not proceed any further if either of the last two data frames are empty
if(nrow(treated_replicate1_dependent.peptides) > 0 &
	nrow(treated_replicate2_dependent.peptides) > 0)

# Filter search results (untreated protein, analysis 1)
{control_replicate1_dependent.peptides <- control_replicate1_all.peptides[
	!is.na(control_replicate1_all.peptides$DP.Mass.Difference) &
	control_replicate1_all.peptides$DP.Mass.Difference < 
		expected.delta.mass + tolerance &
	control_replicate1_all.peptides$DP.Mass.Difference > 
		expected.delta.mass - tolerance &
	control_replicate1_all.peptides$DP.Score > 80 &	
	control_replicate1_all.peptides$DP.PEP < 0.01 &
	control_replicate1_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, control_replicate1_all.peptides$DP.Proteins), ]

# Filter search results (untreated protein, analysis 2)
control_replicate2_dependent.peptides <- control_replicate2_all.peptides[
	!is.na(control_replicate2_all.peptides$DP.Mass.Difference) &
	control_replicate2_all.peptides$DP.Mass.Difference < 
		expected.delta.mass + tolerance & 
	control_replicate2_all.peptides$DP.Mass.Difference > 
		expected.delta.mass - tolerance &
	control_replicate2_all.peptides$DP.Score > 80 &	
	control_replicate2_all.peptides$DP.PEP < 0.01 &
	control_replicate2_all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, control_replicate2_all.peptides$DP.Proteins), ]

# Make empty data frame for sequences, 'DP Probabilities' strings and test outcomes
treated_DP.constancy.results <- data.frame(
	DP.Base.Sequence = character(nrow(treated_replicate1_dependent.peptides)),
	DP.Probabilities.1 = character(nrow(treated_replicate1_dependent.peptides)),
	DP.Probabilities.2 = character(nrow(treated_replicate1_dependent.peptides)),
	In.control_replicate1 = character(nrow(treated_replicate1_dependent.peptides)),
	In.control_replicate2 = character(nrow(treated_replicate1_dependent.peptides)))
for(i in 1:nrow(treated_DP.constancy.results))
	{for(j in 1:ncol(treated_DP.constancy.results))
	treated_DP.constancy.results[i, j] <- NA}

# Perform first round of pairwise comparisons (treated <-> treated), deposit data
for(i in 1:nrow(treated_replicate1_dependent.peptides))
	{for(j in 1:nrow(treated_replicate2_dependent.peptides))
	if(treated_replicate1_dependent.peptides$DP.Base.Sequence[i] ==
		treated_replicate2_dependent.peptides$DP.Base.Sequence[j])
	{treated_DP.constancy.results$DP.Base.Sequence[i] <- 
		treated_replicate1_dependent.peptides$DP.Base.Sequence[i]
	treated_DP.constancy.results$DP.Probabilities.1[i] <- 
		treated_replicate1_dependent.peptides$DP.Probabilities[i]
	treated_DP.constancy.results$DP.Probabilities.2[i] <- 
		treated_replicate2_dependent.peptides$DP.Probabilities[j]}}

# Remove rows of NAs
treated_constant.DPs <- treated_DP.constancy.results[
	!is.na(treated_DP.constancy.results$DP.Base.Sequence), ]

# Perform second round of comparisons (treated <-> untreated), deposit outcomes
if(nrow(control_replicate1_dependent.peptides) > 0)
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 1:nrow(control_replicate1_dependent.peptides))
	if(!is.na(control_replicate1_dependent.peptides$DP.Base.Sequence[j]) &
		treated_constant.DPs$DP.Base.Sequence[i] ==
		control_replicate1_dependent.peptides$DP.Base.Sequence[j])
	treated_constant.DPs$In.control_replicate1[i] <- "Yes"}

# Perform third round of comparisons (treated <-> untreated), deposit outcomes
if(nrow(control_replicate2_dependent.peptides) > 0)
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 1:nrow(control_replicate2_dependent.peptides))
	if(!is.na(control_replicate2_dependent.peptides$DP.Base.Sequence[j]) &
		treated_constant.DPs$DP.Base.Sequence[i] ==
		control_replicate2_dependent.peptides$DP.Base.Sequence[j])
	treated_constant.DPs$In.control_replicate2[i] <- "Yes"}

# Add negative outcomes
for(i in 1:nrow(treated_constant.DPs))
	{for(j in 4:5)
	if(is.na(treated_constant.DPs[i, j]))
	treated_constant.DPs[i, j] <- "No"}

# Isolate 'constantly conjoined' dependent peptides
treated_constantly.conjoined.DPs <- treated_constant.DPs[
	treated_constant.DPs$In.control_replicate1 == "No" &
	treated_constant.DPs$In.control_replicate2 == "No", ]

# Prepare empty matrices
transparency <- matrix(nrow = nchar(protein.sequence), ncol = 11)
reservations <- matrix(nrow = nchar(protein.sequence), ncol = 11)
protein.probabilities <- matrix(0, nrow = nchar(protein.sequence), ncol = 11)

# Write data into matrices via sliding window
for(i in 1:nrow(treated_constantly.conjoined.DPs))
	{for(j in 1:(nchar(protein.sequence) - 
	(nchar(treated_constantly.conjoined.DPs$DP.Base.Sequence[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(
			treated_constantly.conjoined.DPs$DP.Base.Sequence[i]) - 1))

	if(segment == treated_constantly.conjoined.DPs$DP.Base.Sequence[i])
	
	# Define strips
	{footprint <- j:(j + (nchar(
		treated_constantly.conjoined.DPs$DP.Base.Sequence[i]) - 1))
	oversized.footprint <- ((j - 1):
		(j + length(footprint)))[((j - 1):(j + length(footprint))) %in%
			1:nchar(protein.sequence)]
	
	# Parse 'DP Probabilities'
	primary.fragments <- strsplit(
		treated_constantly.conjoined.DPs$DP.Probabilities.1[i],
		split = "(", 
		fixed = TRUE)[[1]]
	peptide.probabilities <- data.frame(
		Sequence = character(length(primary.fragments)), 
		Probability = numeric(length(primary.fragments)),
		Peptide.site = numeric(length(primary.fragments)),
		Protein.site = numeric(length(primary.fragments)))
	peptide.probabilities$Sequence[1] <- primary.fragments[1]
	peptide.probabilities$Peptide.site[1] <- nchar(primary.fragments[1])
		for(k in 2:length(primary.fragments))
		{secondary.fragments <- strsplit(primary.fragments[k], 
			split = ")", 
			fixed = TRUE)[[1]]
		peptide.probabilities$Probability[k - 1] <- as.numeric(
			secondary.fragments[1])
		peptide.probabilities$Sequence[k] <- secondary.fragments[2]
		peptide.probabilities$Peptide.site[k] <- 
			peptide.probabilities$Peptide.site[k - 1] + 
				nchar(peptide.probabilities$Sequence[k])
		peptide.probabilities$Protein.site <- 
			peptide.probabilities$Peptide.site + j - 1}
	
	# Allocate an unreserved strip, write data into matrices
	if(length(which(is.na(reservations[oversized.footprint, 5]))) == 
		length(oversized.footprint))
		{for(k in 1:nrow(peptide.probabilities))
			{for(m in 1:nrow(protein.probabilities))
			if(!is.na(peptide.probabilities$Protein.site[k]) & 
				peptide.probabilities$Protein.site[k] == m)
			protein.probabilities[m, 5] <- 
				peptide.probabilities$Probability[k]}
			reservations[oversized.footprint, 5] <- 1
			transparency[footprint, 5] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 4]))) == 
		length(oversized.footprint))
		{for(k in 1:nrow(peptide.probabilities))
			{for(m in 1:nrow(protein.probabilities))
			if(!is.na(peptide.probabilities$Protein.site[k]) & 
				peptide.probabilities$Protein.site[k] == m)
			protein.probabilities[m, 4] <- 
				peptide.probabilities$Probability[k]}
			reservations[oversized.footprint, 4] <- 1
			transparency[footprint, 4] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 3]))) == 
		length(oversized.footprint))
		{for(k in 1:nrow(peptide.probabilities))
			{for(m in 1:nrow(protein.probabilities))
			if(!is.na(peptide.probabilities$Protein.site[k]) & 
				peptide.probabilities$Protein.site[k] == m)
			protein.probabilities[m, 3] <- 
				peptide.probabilities$Probability[k]}
			reservations[oversized.footprint, 3] <- 1
			transparency[footprint, 3] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 2]))) == 
		length(oversized.footprint))
		{for(k in 1:nrow(peptide.probabilities))
			{for(m in 1:nrow(protein.probabilities))
			if(!is.na(peptide.probabilities$Protein.site[k]) & 
				peptide.probabilities$Protein.site[k] == m)
			protein.probabilities[m, 2] <- 
				peptide.probabilities$Probability[k]}
			reservations[oversized.footprint, 2] <- 1
			transparency[footprint, 2] <- 1} else
	if(length(which(is.na(reservations[oversized.footprint, 1]))) == 
		length(oversized.footprint))
		{for(k in 1:nrow(peptide.probabilities))
			{for(m in 1:nrow(protein.probabilities))
			if(!is.na(peptide.probabilities$Protein.site[k]) & 
				peptide.probabilities$Protein.site[k] == m)
			protein.probabilities[m, 1] <- 
				peptide.probabilities$Probability[k]}
			reservations[oversized.footprint, 1] <- 1
			transparency[footprint, 1] <- 1}}}}

# Calculate required number of panels 
panel.count <- ceiling(nchar(protein.sequence) / 100)

# Choose a graphics device and specify plotting-area dimensions
windows(14, 10)

# Define layout
if(panel.count < 6)
par(mfcol = c(6, 1)) else 
par(mfcol = c(panel.count, 1))

# Specify custom plot margins
par(mar = c(3, 5, 0, 5))

# Prepare the graphic, one panel at a time
for(i in 1:panel.count)

# Make composite image
{image(x = 1:nchar(protein.sequence), 
	z = protein.probabilities,
	breaks = seq(0, 1, 0.05), 
	col = gray(seq(1, 0.05, -0.05)),
	xlim = c(i * 100 - 99.5, 
	i * 100 + 0.5),
	axes = FALSE,
	ann = FALSE)
image(x = 1:nchar(protein.sequence), 
	z = transparency,
	breaks = seq(0, 1, 1), 
	col = rainbow(1, start = 2/3, alpha = 0.5), add = TRUE,
	xlim = c(i * 100 - 99.5, 
	i * 100 + 0.5),
	axes = FALSE,
	ann = FALSE)
axis(side = 1)

# Add a dashed 'backbone'
lines(c(0.5, nchar(protein.sequence) + 0.5), 
	c(0.5, 0.5),
	lty = 3)

	# Solidify observed regions of the backbone
	for(j in 1:length(control_common.peptide.sequences))
		{for(k in 1:(nchar(protein.sequence) -
		(nchar(control_common.peptide.sequences[j]) - 1)))
		{segment <- substr(protein.sequence, 
			start = k, 
			stop = k + (nchar(control_common.peptide.sequences[j]) - 1))
		if(segment == control_common.peptide.sequences[j])
		lines(c(k - 0.5, k + (nchar(
			control_common.peptide.sequences[j]) - 1) + 0.5),
				c(0.5, 0.5),
				lty = 1,
				col = "red")}}
	
	# Add annotations
	for(j in 1:nrow(protein.probabilities))
	if(sum(protein.probabilities[j, ]) > 0)
	text(j, 0.8, labels = paste(substr(protein.sequence, start = j, stop = j),
		as.character(j),
		sep = ""),
		srt = 90)

box()}} else
print("Expected delta mass not observed in at least one analysis of treated protein")

# +-----------------------------------------------------------------------------------+