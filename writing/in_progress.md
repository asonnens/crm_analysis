## Introduction

-Genetic function and molecular evolution are linked. 
-The connection between them is more complex in regulatory DNA, where the rules underpinning how sequence leads to function are less well understood. 

-In coding DNA, redundancy is frequently associated with very rapid evolution-- paralogs diverge extremely rapidly. 
-However, in regulatory DNA, redundancy appears to play a very different role. 
	-Binding site redundancy is a well known aspect of enhancer function, and redundancy of entire 	regulatory modules in some cases appears to be very well conserved between species, although 	in other cases appears to be associated with rapid phenotypic change. 

-Computational identification of enhancers-- history


-To get a better understanding of the role of redundant regulatory elements in Drosophila gene regulation, we used machine learning to predict active enhancers in early Drosophila development, and their expression patterns. 
	-Transcription factor binding and histone modification chip data was used as predictors of 	enhancer activity
	-A set of in-vivo reporters in Drosophila embryos was used as the basis for training a supervised 	random forest.  
	-Accuracy of these predictions was additionally assessed against well annotated enhancers not 	included in the training or testing sets.
	- Enhancers were then predicted in windows around genes known to be active in early 	development. 

-We found...

##Methods

### 1

### 2

### 3


## Results  

### 1 Analysis of features used as predictors
A number of features correlate with enhancer sequences in *Drosophila*, and with their resulting expression patterns. These include transcription factor occupancy [refs], specific histone marks[refs], DNase1 hypersensitivity [refs], evolutionary conservation [refs], and motif enrichment [refs]. Computational strategies that incorporate a combination of these features have been found to be the most successful[refs]. 

The transcription factors Dorsal, Snail and Twist are all known to co-regulate dorsal-ventral patterning in early embryonic development; clustered occupancy of these proteins has been successfully used to identify DV-enhancers [refs]. However, they have also been found to bind enhancers that regulate genes more typically associated with anterior-posterior development, such as gap genes *orthodenticle*, and *knirps*, and pair rule genes *runt* and *hairy* [refs]. For this reason, we also chose to use as predictors transcription factors specifically associated with anterior-posterior regulatory networks, including Bicoid, Caudal, Hairy, Knirps, Kruppel, and Giant [refs]. We also included the 'pioneer' transcription factor Zelda [refs], H3K27ac [refs] and H3K4me1 [refs] histone marks, and the histone acetyltransferase Nejire (also known as p300) [refs].

There are a large number of chromatin immunoprecipitation studies that have made their data publicly available. Notably, in 2009 the Berkeley *Drosophila* Transcription Project released data from ChIP-chip on twenty-one transcription factors involved in early embryonic development. Since then, several studies from the same lab [refs] and others [refs] have released overlapping datasets, frequently obtained from higher resolution[refs] techniques like ChIP-seq[refs]. The regions identified as bound by a given transcription factor using different techniques typically do not completely overlap (figure, refs). Even regions identified using the same technique by different labs, or using different antibodies even within a given experiment frequently produce variable results (figure, results). For this reason, we chose to include multiple datasets for some features. Our final list included twenty-three predictors, encompassing fourteen unique features (table1).

Interestingly, in addition to Dorsal, Snail and Twist binding to enhancers for AP development genes, transcription factors typically associated with AP development also appear to bind to enhancers for DV genes. It is possible these networks are both regulated by an overlapping set of transcription factors [refs]. It is also possible that this shared occupancy is also due to promiscuous binding at regions of open chromatin [refs].

***Figure 1*** ---------- table of features, and a couple of relevant venn diagrams


### 2 Characteristics of training set for supervised machine learning
-The size of the training set
-how Stark active and inactive reporters correlate with features
***Figure 2*** ---- heatmaps
-description of my annotations of expression patterns
-Example images (maybe should go in methods?)
***Figure 3---*** 

### 3 Predicting enhancer activity
-description of random forest
-which features are the most important in this analysis
-dropping any one feature does not improve results, also does not dramatically hurt results
***Figure 4*** -make this figure to reflect above two points
-including different subsets of the Stark reporters in training and testing dramatically influences results
-what results look like with no filtering
-what results look like after filtering based on presence of clustered TF binding
***Figure 5a*** - some kind of compilation figure of ROC and PR data?

### 4 Predicting enhancer expression patterns
-how methods for this do or don't differ from methods for predicting activity
-description of random forest
-which features are most important (surprisingly, mostly the same ones)
***Figure 5b***

De-novo predictions around genes known to be active
-potential enhancers are identified based on clustered TF occupancy
-enhancer defined as window of X-bp around center of TF cluster
***Figure 6*** --- flow chart! (this should maybe go in methods?)

How well predictions recover known regulatory elements
How well predictions recover known regulatory elements using different sets of assumptions
-Assumptions
	-the subset of the Stark dataset used for training the random forest
	-random forest settings
	-how clusters of TF binding are classified as putative enhancers
	-the size potential enhancers are set at
***Figure 7*** -some kind of compilation figure for above 


Using top three models, predictions around genes involving dorsal ventral activity
-I've only done this with one model...


## References

















