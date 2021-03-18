# +-----------------------------------------------------------------------------------+
# | Supporting information file S5 Script                                             |
# +-----------------------------------------------------------------------------------+
# | Script V                                                                          |
# |                                                                                   |
# | To accompany the report 'Visualisation tools for dependent peptide searches to    |
# | support the exploration of in vitro protein modifications' by G. W. Preston, L.   |
# | Yang, D. H. Phillips and C. S. Maier                                              |
# +-----------------------------------------------------------------------------------+

#   Usage notes:
#   1. This script isolates dependent peptides with a specified mass shift and 
#      generates a probability localisation plot.
#   2. The script requires:
#      i.   R version 3.6.0 or later;
#      ii.  Windows (a Windows-specific graphics device will be used);
#      iii. package 'seqinR';
#      iv.  one allPeptides.txt file;
#      v.   a *.fasta file containing the sequence of a protein of interest;
#      vi.  the ID (e.g. UniProt ID) of the protein of interest, exactly as it 
#           appears in the allPeptides.txt file;
#      vii. an expected mass shift and a tolerance. We use a tolerance of 0.01 Da.  
#   3. Data from our BSA study, including the allPeptides.txt and *.fasta files 
#      specified in the script, are available via PRIDE/ProteomeXchange 
#      (accession PXD013040).
#   4. Further instructions are provided within the script. Comments beginning with a
#      greater-than symbol (>) are instructions for users.
#   5. Detailed explanatory notes can be found in Tables B and C of supporting 
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

# Get search results
# > Specify the path to an allPeptides.txt file
all.peptides <- read.table(paste("C:", "\\", "Users", "\\", "George", "\\",
	"Documents", "\\", "Model data", "\\", "treated_replicate1_allPeptides.txt", 
	sep = ""), 
	sep = "\t", 
	header = TRUE,
	fill = TRUE)

# Make a non-redundant list of observed peptide sequences
peptide.sequences <- unique(c(all.peptides$Sequence[grepl(protein.ID, 
	all.peptides$Proteins) &
	all.peptides$Score > 80],
	all.peptides$DP.Base.Sequence[
		all.peptides$DP.Score > 80 &
		all.peptides$DP.PEP < 0.01 &
		all.peptides$DP.Decoy != "+" &
		grepl(protein.ID, all.peptides$DP.Proteins)]))

# > Input expected mass shift
expected.delta.mass <- 143.05825

# > Specify tolerance
tolerance <- 0.01

# Filter search results
dependent.peptides <- all.peptides[
	!is.na(all.peptides$DP.Mass.Difference) &
	all.peptides$DP.Mass.Difference < 
		expected.delta.mass + tolerance &
	all.peptides$DP.Mass.Difference > 
		expected.delta.mass - tolerance &
	all.peptides$DP.Score > 80 &	
	all.peptides$DP.PEP < 0.01 &
	all.peptides$DP.Decoy != "+" &
	grepl(protein.ID, all.peptides$DP.Proteins), ]

# Do not proceed if there are no dependent.peptides
if(nrow(dependent.peptides) > 0)

# Prepare empty matrices
{transparency <- matrix(nrow = nchar(protein.sequence), ncol = 11)
reservations <- matrix(nrow = nchar(protein.sequence), ncol = 11)
protein.probabilities <- matrix(0, nrow = nchar(protein.sequence), ncol = 11)

# Write data into matrices via sliding window
for(i in 1:nrow(dependent.peptides))
	{for(j in 1:(nchar(protein.sequence) - 
	(nchar(dependent.peptides$DP.Base.Sequence[i]) - 1)))
	{segment <- substr(protein.sequence, 
		start = j, 
		stop = j + (nchar(
			dependent.peptides$DP.Base.Sequence[i]) - 1))

	if(segment == dependent.peptides$DP.Base.Sequence[i])

	# Define strips
	{footprint <- j:(j + (nchar(
		dependent.peptides$DP.Base.Sequence[i]) - 1))
	oversized.footprint <- ((j - 1):
		(j + length(footprint)))[((j - 1):(j + length(footprint))) %in%
			1:nchar(protein.sequence)]
	
	# Parse 'DP Probabilities'
	primary.fragments <- strsplit(dependent.peptides$DP.Probabilities[i],
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

# Define custom plot margins
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
	for(j in 1:length(peptide.sequences))
		{for(k in 1:(nchar(protein.sequence) -
		(nchar(peptide.sequences[j]) - 1)))
		{segment <- substr(protein.sequence, 
			start = k, 
			stop = k + (nchar(peptide.sequences[j]) - 1))
		if(segment == peptide.sequences[j])
		lines(c(k - 0.5, k + (nchar(
			peptide.sequences[j]) - 1) + 0.5),
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
print("Expected delta mass not observed")

# +-----------------------------------------------------------------------------------+