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
We identified potential enhancers de-novo based on clustering of transcription factors within a fifty kb window around genes of interest. Genes of interest were selected based on an annotated mesodermal or neurogenic ectodermal expression pattern, as based on in situs of the Berkeley Drosophila Genome Project. Several parameters of varying selectivity were used as a threshold for 'clustering', to generate sets of putative enhancers. Thresholds tested included the presence of at least two transcription factor peaks, the presence of at least three transcription factor peaks, or the presence of peaks for Dorsal, Snail and Twist. These regions were then classified as 'putative enhancers' or 'background binding' by a random forest algorithm trained on regions from the Fly Enhancer Resource. Enhancers were also identified and classified around several genes with very well understood regulation, with the intent of using these enhancers as a secondary validation step. 

### Training Set

The Stark Lab Fly Enhancers Resource includes X reporters, Y of which show activity. However, we filtered for regions that are highly active in stage 4-6 embryos, limiting our 'active' set to __ reporters. The Fly Enhancers resource also includes annotations of expression patterns; however, at stage 4-6 embryos this includes Q unique patterns, averaging __ reporters per expression pattern. To increase the size of the expression pattern training sets, we manually re-annotated the expression patterns for all highly active reporters based on photos of in-situs from the Stark website. We used using broader categories, dividing the embryo into thirds along the AP and DV axes. (***show image? cartoon? Make this***).

All reporters were cross-compared with various genomics datasets to determine the characteristic features of inactive reporters and active reporters associated with various expression patterns. These included data for evolutionary conservation, transcription factor occupancy [refs], specific histone marks[refs], and motif enrichment [refs]-- Table 1. 

Phylogenetic footprinting has been shown to be an effective means of identifying enhancers in some species. Sequence conservation has also been observed in some Drosophila regulatory elements, but using it as a predictor is complicated by the generally rapid rate of divergence in regulatory sequences, and the compact and largely conserved nature of the Drosophila genome. Looking at pairwise conservation between various Drosophila species and Drosophila melanogaster in reporters, there is no visible relationship between conservation and activity. A comparison of reporter sequences that are active at all embryonic stages relative to sequences that are inactive across all embryonic stages does not show a visible signal. 
***Figure 1*** - see figure-- shows conservation in simulans, yakuba, pseudoobsura, and grimshawi in all reporters, relative to melanogaster.
https://raw.githubusercontent.com/asonnens/crm_analysis/master/writing/figure1/Active_Inactive_conservation.png

There are a large number of publicly available genomic datasets, through ModEncode, UCSC, the Berkeley Drosophila Genome Project, and other open access resources. We chose to use occupancy information from a range of transcription factors and several histone marks. Dorsal, Snail and Twist are all known to co-regulate dorsal-ventral patterning in early embryonic development; clustered occupancy of these proteins has been successfully used to identify DV-enhancers [refs]. However, they have also been found to bind enhancers that regulate genes more typically associated with anterior-posterior development, such as gap genes *orthodenticle*, and *knirps*, and pair rule genes *runt* and *hairy* [refs]. For this reason, we also chose to use as predictors transcription factors specifically associated with anterior-posterior regulatory networks, including Bicoid, Caudal, Hairy, Knirps, Kruppel, and Giant [refs]. We also included the 'pioneer' transcription factor Zelda [refs], H3K27ac [refs] and H3K4me1 [refs] histone marks, and the histone acetyltransferase Nejire (also known as p300) [refs]. There are a large number of chromatin immunoprecipitation studies that have made their data publicly available. Notably, in 2009 the Berkeley *Drosophila* Transcription Project released data from ChIP-chip on twenty-one transcription factors involved in early embryonic development. Since then, several studies from the same lab [refs] and others [refs] have released overlapping datasets, frequently obtained from higher resolution[refs] techniques like ChIP-seq[refs]. The regions identified as bound by a given transcription factor using different techniques typically do not completely overlap (***figure 2*** part 1 representative venn diagrams), and their resulting scores do not perfectly correlate 
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_1/Venn_diagrams.png?raw=true
(***figure 2*** part 2- compilation of correlation heatmaps). 
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/AP_TFs_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/DV_TFs_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/Dl_motif_chip_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/Sna_motif_chip_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/Twi_motif_chip_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/chromatin_chip_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/evolution_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/motif_correlation_heatmap.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure2/Part_2/all_chip_correlation_heatmap.png
Even regions identified using the same technique by different labs, or using different antibodies even within a given experiment frequently produce variable results (refs, ***figure 2*** part 3 -make this figure-- venn diagrams for different antibodies same lab). For this reason, we chose to include multiple datasets for some features. Our final list of ChIP datasets included twenty-three predictors, encompassing fourteen unique features. 

Interestingly, in addition to Dorsal, Snail and Twist binding to enhancers for AP development genes, transcription factors typically associated with AP development also appear to bind to enhancers for DV genes, and peak scores for these factors do not clearly correlate with AP enhancers than with DV enhancers. It is possible these networks are both regulated by an overlapping set of transcription factors [refs]. It is also possible that this shared occupancy is also due to promiscuous binding at regions of open chromatin [refs]. (***figure 2*** part 3 make this figure-- screen shots of places where AP TFs bind at unexpected genes)

Motifs have a relationship with binding of transcription factors and with regulatory function. We used PWMs associated with Dorsal, Snail and Twist binding from the JASPAR database, and assessed the number of high-affinity matches in the Stark regulatory regions using MAST. This led to a final set of predictors which included 39 unique features ( ***Table 1*** table of features, and their provenance)   Interestingly, there was only minimal correlation between scores for different motifs of the same transcription factor, and even less between motifs for a transcription factor and actually occupancy in stage 4-6 embryos as measured by ChIP (***Figure 2*** part 2 with the rest of the correlation heat maps, currently).

Although random forest algorithms are generally robust to co-linearity between features and non-normality of predictors, we scaled feature scores for all predictors around 0, and also created secondary datasets wherein features that are highly correlated with each other were dropped, and datasets in which the first 10 and 25 principal components were used instead of scores for individual features.
***Table 2***, ***figure 3*** Plot of PCA1 x PCA2
(table of covariation of features)
https://github.com/asonnens/crm_analysis/blob/master/writing/figure3/PCA_reporters.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure3/scree_plot_supplement.png

### Predicting enhancer activity

We tested the efficacy of different features using a random forest algorithm trained and tested on a subsets of the reporters from the Fly Enhancer Resource. 
 
Using all highly active reporters, and all inactive reporters as a pool from which to sample training and test sets, a random forest algorithm was able to predict activity of the test set with X degree accuracy, and Y AUC. However, as there are ~7000 inactive reporters at this stage of development, relative to ~400 active reporters, this is a highly unbalanced dataset, which biases most machine learning algorithms to overestimate the occurrence of the majority event. A precision recall curve shows this to be the case-- area under precision recall curves is substantially lower than under the ROC curve (***figure 4*** part 1- ROC and PR curves for kitchen sink approach). 
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part1/ROC_all_data_scale.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part1/PR_all_data_scale.png
However, removing highly correlated elements or employing the first 25 principal components (which make up 99% of the predictors variation), have no significant impact on area under ROC or precision-recall curves, and only employing the first ten PCs substantially hurts performance on precision-recall (***figure 4*** part 2-- bargraphs of area under curves). 
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part2/AUC_various_variables.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part2/PRAUC_various_variables.png
As reducing the number of features did not improve predictions and made it more difficult to interpret biological significance, we used all 39 features for further analyses.
Interestingly, when using all features, or excluding highly correlated features, Zelda binding, Bicoid binding, and Snail binding are by far the most predictive features. Surprisingly (given figure 1), conservation with distantly related species of Drosophila also contributes to predictive accuracy, more than conservation in more closely related species or many other features. (***figure 4*** part 3-- importance plots).
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part3/Importance_all_data.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure4/Part3/Importance_reduced_data.png

Many of the inactive regions in stage 4-6 embryos have very little transcription factor occupancy. These regions are less useful for distinguishing between clusters resulting from background binding and truly active regions. As putative enhancers identified by clustered transcription factor occupancy by definition are bound by multiple transcription factors, we sought to improve predictive power by restricting the training set to regions that are similarly bound. We filtered the training and testing set by several parameters: to only include reporters that are co-bound by Dorsal, Snail and Twist, or to only include reporters that are bound by at least two transcription factors, at least one of which is Dorsal, Snail, Twist, or Zelda. In each case we also balanced uneven datasets by selecting equal numbers of the most highly bound active and inactive reporters (as measured by Euclidean Distance) from each filtered dataset. Based on AUC of ROC curves, balanced (blue) and unbalanced (gray) datasets were approximately equally effective at training the random forest, and more restrictive filtering by multiple transcription factors lowered overall accuracy. Based on AUC of precision-recall graphs, balanced datasets are significantly more effective than unbalanced datasets. ***figure 5*** bar charts of AUC for ROC and PR curves, various models
https://github.com/asonnens/crm_analysis/blob/master/writing/figure5/Figure5_part1/AUC_various_models.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure5/Figure5_part1/PRAUC_various_models.png

### Testing against validated enhancers
We tested each random forest on a list of enhancers predicted in a 50kb window around genes with well studied regulation in stage 4-6 embryos. This was done with the assumption that annotated regulatory regions around these genes are effective measures of true-positives, and that as the regions have been fairly well studied, additional enhancers identified within this window are likely to be false positives. 

Putative enhancers were identified by clustering of transcription factors. We used several different criteria-- clustered binding of Dorsal, Snail, and Twist, or binding of at least two transcription factors with one of them being Dorsal, Snail, Twist, or Zelda. Surprisingly, in every case random forest trained on a balanced datasets identified 100% of clusters as active. 

Random forests trained on unbalanced datasets did recover a significant number of false positives, although relatively fewer. Models trained on reporters filtered by the presence of Zelda binding were able to correctly identify the highest proportion of annotated enhancer sequence, relative to un-annotated regions. ***figure 5 part 2*** these bar graphs are hideous, I need to compile them into something more attractive.


### Analysis of feature importance

Importance analysis of this random forest shows that as in other models, Zelda, Bicoid and Snail binding are the three most important parameters for correctly identifying active regions, and conservation in distantly related Drosophila species is more important than conservation in closely related species. We evaluated the importance of categories of features by training random forests on the same pool of reporters, but leaving out categories of predictors (e.g. all conservation data, all motif data, etc.). Based on AUC of ROC curves, leaving out any single category of features had no effect on prediction success. However, based on AUC of precision-recall curves, leaving out Zelda binding or occupancy information for Dorsal-Ventral transcription factors (Dorsal, Snail, and Twist) caused a substantial drop in accuracy.
***figure 6*** bar charts of leave-factor-out analysis
https://github.com/asonnens/crm_analysis/blob/master/writing/figure6/Best_model_importance.png?raw=true
https://github.com/asonnens/crm_analysis/blob/master/writing/figure6/Dropping_features_AUC.png
https://github.com/asonnens/crm_analysis/blob/master/writing/figure6/Dropping_features_PRAUC.png


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

















