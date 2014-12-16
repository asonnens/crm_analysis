###############################################################################
#Anne Sonnenschein
#2014
#this script applies classes to find 'average' characteristics of a fixed group
#of well defined enhancers. It applies functions from the enhancer classes
#to identify putative enhancers around genes that are involved in a developmental
#pathway, and genes that are not involved, and act as negative control
###############################################################################

import glob
import enhancer_classes
import math

protein_list = ["Dorsal", "Snail", "Twist", "Zelda"]

def known_enhancers(infilename):
#########################################################
#this function takes a list of 'known enhancers' supplied
#by the user
#for each 'known enhancer it goes through the files of
#genomic information and extracts features associated
#with the genomic coordinates containing the enhancer
#########################################################
    known_enhancers = []
    infile = open(infilename)
    for line in infile:
        line_list = line.split(",")
        gene = line_list[0]
        if "chr" not in line_list[1]:
            chr = "chr" + line_list[1]
        else:
            chr = line_list[1]
        coord_1 = int(line_list[2])
        coord_2 = int(line_list[3])
        TSS = int(line_list[4])
        distance_to_TSS = 0
        if math.fabs(coord_1 - TSS) <= math.fabs(coord_2- TSS):
            distance_to_TSS =  math.fabs(coord_1 - TSS)
        else:
            distance_to_TSS =  math.fabs(coord_2 - TSS)
        expr = line_list[5].strip()
        each_enhancer = enhancer_classes.Enhancer(\
                                                  gene,str(chr + "\t" + str(coord_1) + "\t" + \
                                                  str(coord_2)),distance_to_TSS, expr,"confirmed")
        known_enhancers.append(each_enhancer)
    infile.close()
    for i in known_enhancers:
        for bedgraph in enhancer_classes.bedfiles:
            i.read_file(bedgraph)
        for wigfile in enhancer_classes.phast_files:
            i.read_file(wigfile)
        for wigfile in enhancer_classes.histone_files:
            i.read_file(wigfile)
    print "getting info for known enhancers"
    #this function gets the features for known enhancers
    known_enhancer_comparison = enhancer_classes.Enhancer_compare(known_enhancers, "known_enhancers")
    return known_enhancer_comparison

def genes_to_check(infilename):
#########################################################
#this function takes a list of genes supplied by the user
#the genes should all be associated with a specific 
#developmental pathway (preferably identified by their 
#in-situ expression
#it uses the coordinates for each gene region to search the 
#genomic information feature files to extract features in the
#that overlap, and identify putative enhancers
#########################################################
    genes_to_check = []
    putative_enhancers = []
    negative_controls = []
    infile = open(infilename)
    for line in infile:
        line_list = line.split(",")
        gene = line_list[0]
        if "chr" not in line_list[1]:
            chr = "chr" + line_list[1]
        else:
            chr = line_list[1]
        gene_coord_1 = int(line_list[2])
        gene_coord_2 = int(line_list[3])
        TSS = gene_coord_1
        if line_list[4] == "+":
            TSS = gene_coord_1
        else:
            TSS = gene_coord_2
        region_coord_1 = TSS - 25000
        region_coord_2 = TSS + 25000
        region_coords = str(chr + "\t" + str(region_coord_1) + "\t" + str(region_coord_2))
        status = line_list[5]
        expr = line_list[6].strip()
        print "finding putative enhancers"
        #this function finds features associated with the genes genomic region
        putative_enhancers = enhancer_classes.Enhancer_find(protein_list, region_coords, gene, TSS, expr, status)
        #this function creates putative enhancers from the list of features, finding features that cluster with
        #a seed feature. At the moment, the seed feature is set to Dorsal binding
        enhancer_classes.Enhancer_find.find_TF(putative_enhancers, "Dorsal")

      
known_enhancer_comparison = known_enhancers("/mnt/home/sonnens2/crm_analysis/input_data/input_files/known_enhancers.csv")
enhancer_classes.Enhancer_compare.print_enhancer_data(known_enhancer_comparison)

#finding potential enhancers for genes that are definitely involved in the developmental pathway
#genes_to_check("/mnt/home/sonnens2/crm_analysis/input_data/input_files/negative_enhancers.csv")
#finding potential enhancers for genes that ARE involved in the developmental pathway
genes_to_check("/mnt/home/sonnens2/crm_analysis/input_data/input_files/target_genes.csv")
    



