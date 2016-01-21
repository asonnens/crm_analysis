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
        for axt in enhancer_classes.axt_files:
            i.read_file(axt)
    print "getting info for known enhancers"
    #this function gets the features for known enhancers
    known_enhancer_comparison = enhancer_classes.Enhancer_compare(known_enhancers, "known_enhancers")
    return known_enhancer_comparison

def genes_to_check(infilename, outfilename, windowradius, seed):
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
    outputlist = []
    infile = open(infilename)
    for line in infile:
        line_list = line.split(",")
        print line
        gene = line_list[0]
        if "chr" not in line_list[1]:
            chr = "chr" + line_list[1]
        else:
            chr = line_list[1]
        gene_c_1 = int(line_list[2])
        gene_c_2 = int(line_list[3])
        gene_coord_1 = gene_c_1
        gene_coord_2 = gene_c_2
        if gene_c_1 > gene_c_2:
            gene_coord_1 = gene_c_2
            gene_coord_2 = gene_c_1
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
        print gene, region_coords
        print protein_list
        putative_enhancers = enhancer_classes.Enhancer_find(protein_list, region_coords, gene, TSS, expr, status, windowradius)
        #print putative_enhancers
        #this function creates putative enhancers from the list of features, finding features that cluster with
        #a seed feature. The seed feature default is set to Dorsal chip-chip binding
        return_list = enhancer_classes.Enhancer_find.find_TF(putative_enhancers, seed)
        for i in return_list:
            outputlist.append(i)
    print len(outputlist)
    output_output = enhancer_classes.Enhancer_compare(outputlist, outfilename)
    enhancer_classes.Enhancer_compare.print_enhancer_data(output_output)

      
#known_enhancer_comparison = known_enhancers("/mnt/home/sonnens2/crm_analysis/input_data/input_files/known_enhancers.csv")
#enhancer_classes.Enhancer_compare.print_enhancer_data(known_enhancer_comparison)
#crap = known_enhancers("/mnt/home/sonnens2/crm_analysis/update/input_data/input_files/my_crap_list.csv")
#enhancer_classes.Enhancer_compare.print_enhancer_data(crap)
#finding potential enhancers for genes that are definitely involved in the developmental pathway
#genes_to_check("/mnt/home/sonnens2/crm_analysis/input_data/input_files/negative_enhancers.csv")
#finding potential enhancers for genes that ARE involved in the developmental pathway
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Furlong_snail_target_500", 250, "snail")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Furlong_snail_prospective_500", 250, "snail")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Furlong_twist_target_500", 250, "twist")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Furlong_twist_prospective_500", 250, "twist")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Furlong_snail_target_1000", 500, "sna")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Furlong_snail_prospective_1000", 500, "sna")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Furlong_twist_target_1000", 500, "twist")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Furlong_twist_prospective_1000", 500, "twist")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Rushlow_Dorsal_target_500", 250, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Rushlow_Dorsal_prospective_500", 250, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Dorsal2009_target_500", 250, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Dorsal2009_prospective_500", 250, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TESTRushlow_Dorsal_target_1000", 500, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Rushlow_Dorsal_prospective_1000", 500, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TESTDorsal2009_target_1000", 500, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "TEST_Dorsal2009_prospective_1000", 500, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Zelda_target_1000", 500, "Zelda180")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Zelda_prospective_1000", 500, "Zelda180")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "Zelda_target_500", 250, "Zelda180")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "Zelda_prospective_500", 250, "Zelda180")

#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TESTRushlow_Dorsal_target_500", 250, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "TESTRushlow_Dorsal_prospective_500", 250, "DorsalR")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TESTDorsal2009_target_500", 250, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "TEST_Dorsal2009_prospective_500", 500, "dl")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TEST_Furlong_twist_target_500", 250, "twist")
#genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "TEST_Furlong_twist_prospective_500", 250, "twist")
genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/target_genes.csv", "TESTZelda_target_500", 250, "Zelda180")
genes_to_check("/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/prosp_genes.csv", "TESTZelda_prospective_500", 250, "Zelda180")



