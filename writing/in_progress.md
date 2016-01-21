## Abstract


## Introduction


-Regulatory DNA architecture is a subject that is still not fully understood
- Most data points to it being far more complex than originally thought
	- superenhancers, shadow enhancers, etc.
- it's difficult to know non-coding function without extensive in-vivo testing
- enhancer prediction based on genomic cues is the future
- simplistic enhancer prediction can be misleading-- e.g. calling all DNase1 hypersensitive regions as enhancers
- Machine learning and computational enhancer prediction
- Stark dataset provides unique opportunity for testing limits of computational enhancer prediction
- etc. 


-We found...

##Methods




## Results  

### Identifying potential enhancers ###
We chose our genes of interest based on an annotated mesodermal or neurogenic ectodermal expression pattern, as based on in situs of the Berkeley Drosophila Genome Project. We identified potential enhancers de-novo based on clustering of transcription factors within a fifty kb window around genes of interest. We used several parameters as a threshold for 'clustering'-- the presence of at least two transcription factor peaks, the presence of at least three transcription factor peaks, or the presence of peaks for Dorsal, Snail and Twist. These regions were then classified as putative enhancers or background binding by a random forest algorithm trained regions from the Fly Enhancer Resource. We also identified and classified enhancers around several genes with very well understood regulation, with the intent of using these enhancers as a secondary validation step. 

### Training Set

The Stark Lab Fly Enhancers Resource includes X reporters, Y of which show activity. However, we filtered for regions that are highly active in stage 4-6 embryos. only Z of these show activity. The Fly Enhancers resource also includes annotations of expression patterns. However, at stage 4-6 embryos this included Q unique patterns, which led to small sample sizes per expression pattern. This is further complicated for the purposes of machine learning by the fact that many embryos were assigned multiple expression patterns. For our analysis we manually re-annotated the expression patterns for all highly active reporters, based on photos of in-situs from the Stark website, using broader categories (anterior, central and posterior expression, and ectoderm, mesoderm, and endoderm expression).

All reporters were cross-compared with various genomics datasets, to determine the characteristic features of inactive reporters and active reporters with various expression patterns. These included data for evolutionary conservation, transcription factor occupancy [refs], specific histone marks[refs], and motif enrichment [refs]-- Table 1. 

Phylogenetic footprinting has been an effective method of identifying enhancers in some species. Sequence conservation has also been observed in some Drosophila regulatory elements, but using it as a predictor is complicated by the generally rapid rate of divergence in regulatory sequences, and the compact and largely conserved nature of the Drosophila genome. Looking at pairwise conservation between various Drosophila species and Drosophila melanogaster in reporters, there is no visible relationship between conservation and activity. A comparison of reporter sequences that are active at all embryonic stages relative to sequences that are inactive across all embryonic stages does not show a visible signal. 
***Figure 1***

There are a large number of publicly available chromatin immunoprecipitation datasets, through ModEncode, the Berkeley Drosophila Genome Project, and others. We chose to use occupancy information from a range of transcription factors and several histone marks. Dorsal, Snail and Twist are all known to co-regulate dorsal-ventral patterning in early embryonic development; clustered occupancy of these proteins has been successfully used to identify DV-enhancers [refs]. However, they have also been found to bind enhancers that regulate genes more typically associated with anterior-posterior development, such as gap genes *orthodenticle*, and *knirps*, and pair rule genes *runt* and *hairy* [refs]. For this reason, we also chose to use as predictors transcription factors specifically associated with anterior-posterior regulatory networks, including Bicoid, Caudal, Hairy, Knirps, Kruppel, and Giant [refs]. We also included the 'pioneer' transcription factor Zelda [refs], H3K27ac [refs] and H3K4me1 [refs] histone marks, and the histone acetyltransferase Nejire (also known as p300) [refs]. There are a large number of chromatin immunoprecipitation studies that have made their data publicly available. Notably, in 2009 the Berkeley *Drosophila* Transcription Project released data from ChIP-chip on twenty-one transcription factors involved in early embryonic development. Since then, several studies from the same lab [refs] and others [refs] have released overlapping datasets, frequently obtained from higher resolution[refs] techniques like ChIP-seq[refs]. The regions identified as bound by a given transcription factor using different techniques typically do not completely overlap (figure 2-- part 2), and their resulting scores do not perfectly correlate (figure 2-- part 1). Even regions identified using the same technique by different labs, or using different antibodies even within a given experiment frequently produce variable results (refs, figure 2-- part 2). For this reason, we chose to include multiple datasets for some features. Our final list of ChIP datasets included twenty-three predictors, encompassing fourteen unique features. 

Interestingly, in addition to Dorsal, Snail and Twist binding to enhancers for AP development genes, transcription factors typically associated with AP development also appear to bind to enhancers for DV genes, and peak scores for these factors do not clearly correlate with AP enhancers than with DV enhancers. It is possible these networks are both regulated by an overlapping set of transcription factors [refs]. It is also possible that this shared occupancy is also due to promiscuous binding at regions of open chromatin [refs]. (figure 2 part 3)

Motifs have a relationship with binding of transcription factors and with regulatory function. We used PWMs associated with Dorsal, Snail and Twist binding from the JASPAR database, and assessed the number of high-affinity matches in the Stark regulatory regions using MAST. This led to a final set of predictors which included 39 unique features ( ***Table 1***) Interestingly, there was only minimal correlation between scores for different motifs of the same transcription factor, and even less between motifs for a transcription factor and actually occupancy in stage 4-6 embryos as measured by ChIP (Figure 2 part 4).

Although random forest algorithms are generally robust to co-linearity between features and non-normality of predictors, we scaled feature scores for all predictors around 0, and also created secondary datasets wherein features that are highly correlated with each other were dropped, and datasets in which the first 10 and 25 principal components were used instead of scores for individual features.
***Figure 3***, ***Table 2***
(table of covariation of features)


### Predicting enhancer activity
We tested the efficacy of different features using a random forest algorithm trained on a subset of the reporters from the Fly Enhancer Resource on both a test set based on reporters from the Fly Enhancer.
[Using the kitchen sink approach, in which two thirds of all highly active and inactive regions from stage 4-6 embryos were used to train a random forest algorithm, and one third was used as testing, results were moderately successful on the test dataset]
[Using various subsets, these variations performed the best on the test dataset]
[reduced correation or principal components 
[of the ones that perform the best, this is how leaving out one or many features influences performance]
[of the top few subsets, they performed X well on various version of predicted enhancers around well understood genes]


### Predicting enhancer expression patterns
-how methods for this do or don't differ from methods for predicting activity
-description of random forest
-which features are most important (surprisingly, mostly the same ones)

***[to complete-- how these algorithms perform at predicting expression patterns]***
***[likely compound error rates]***
***[how this compares to Furlong's methods accuracy]***


-Assumptions
	-the subset of the Stark dataset used for training the random forest
	-random forest settings
	-how clusters of TF binding are classified as putative enhancers
	-the size potential enhancers are set at
***Figure 7*** -some kind of compilation figure for above 

### Shadow enhancers etc. ###

## Discussion

## References

















