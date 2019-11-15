####################################################################################################################################
#		 	ASSIGNMENT ANALYSIS
# 						from
# "Biodiversity inventory of the grey mullets (Actinopterygii: Mugilidae) of the Indo-Australian Archipelago 
# through the iterative use of DNA-based species delimitation and specimen assignment methods"
# 
# Delrieu-Trottin E, Durand JD, Limmon G, Sukmono T, Kadarusman, Sugeha HY, Chen WJ, Busson F, Borsa P, Dahruddin H, 
# Sauri S, Fitriana Y, Zein MSA,  Hocdé R, Pouyaud L, Keith P, D Wowor, Steinke D, Hanner R, Hubert N
#
# 
#
######################################################################################################################################

# Install the necessary packages running:
install.packages("BarcodingR")
install.packages("adegenet")

# Load the two packages in R:
library("BarcodingR")
library("adegenet")

##############################################################################
#
# ASSIGNMENT ANALYSIS
#
##############################################################################


# One should prepare two FASTA format files, corresponding to your aligned 
# (1) reference sequences and (2) query sequences.

### import reference sequences using adegenet-- KNOWN , aligned; only 1 sequence per OTU
ref<-fasta2DNAbin("your_reference_sequences.fasta")

### import query sequences using adegenet -- UNKNOWN, aligned
que<-fasta2DNAbin("your_query_sequences.fasta")


################### 	Run FZKMER   #########################################
FZKMER <- barcoding.spe.identify2(ref, que, kmer = 5, optimization = "TRUE")

# to add a column specifying the method
FZKMER $output_identified$methods <-"FZKMER"

# to export results in a table that can be opened with Office:
write.csv(FZKMER$output_identified, "FZKMER.csv",row.names=FALSE)

################### 	Run FZ   #############################################
FZ <- barcoding.spe.identify(ref, que, method = "fuzzyId")

# to add a column specifying the method
FZ$output_identified$methods <-"FZ"

# to export results in a table that can be opened with Office:
write.csv(FZ$output_identified, "FZ.csv",row.names=FALSE)

################### 	Run BP   #############################################
### !!!! BEWARE - can be quite long to run !!!! 
BP <- barcoding.spe.identify(ref, que, method = "bpNewTraining")

# to add a column specifying the method
BP$output_identified$methods <-"BP"

# to export results in a table that can be opened with Office:
write.csv(BP$output_identified, "BP.csv",row.names=FALSE)

###############################################################################


# to export ALL results in a table that can be opened with Office:

data_all <- cbind(BP$output_identified, FZ$output_identified, FZKMER$output_identified)
write.csv(data_all, "Results_3_Assignment_methods.csv",row.names=FALSE)

#######################################################################################

