# crm_analysis
This is code and data of a manageable size for my shadow-enhancer-prediction project.

scripts in code file include:
	Python scripts for annotating known enhancers with genomic data pulled out of bedgraphs
	Python scripts for predicting enhancers based on clusters of transcription factor binding
	R scripts for training and testing how well enhancer activity can be predicted 
		currently titled  "crm_version number.R"-- I need to come up with something better

Machine learning input:
	-Genomic information (primarily taken from MacArthur et al. 2009 in the Eisen lab) for reporters identified by Kvon et al. 2014 in the Stark lab
	-Genomic information for putative enhancers, identified by overlap of Dorsal binding + one additional transcription factor
Machine learning output:
	-misc outputs from my R scripts

output_data
	-output from enhancer prediction scripts