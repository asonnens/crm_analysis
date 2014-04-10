#this script is intended to determine the 'average' characteristics of a fixed group of well defined enhancers

import glob
import enhancer_classes
import math

protein_list = ["Dorsal", "Snail", "Twist", "Zelda"]

def known_enhancers(infilename):
#########################################################
#
#
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
        for wigfile in enhancer_classes.wigfiles:
            i.read_file(wigfile)
    known_enhancer_comparison = enhancer_classes.Enhancer_compare(known_enhancers, "known_enhancers")
    return known_enhancer_comparison

def genes_to_check(infilename):
#########################################################
#
#
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
        region_coord_1 = TSS - 50000
        region_coord_2 = TSS + 50000
        region_coords = str(chr + "\t" + str(region_coord_1) + "\t" + str(region_coord_2))
        status = line_list[5]
        expr = line_list[6].strip()
        putative_enhancers = enhancer_classes.Enhancer_find(protein_list, region_coords, gene, TSS, expr, status)
        enhancer_classes.Enhancer_find.find_TF(putative_enhancers, "Dorsal")

      
known_enhancer_comparison = known_enhancers("/mnt/home/sonnens2/crm_analysis/input_files/known_enhancers.csv")
enhancer_classes.Enhancer_compare.print_enhancer_data(known_enhancer_comparison)


genes_to_check("/mnt/home/sonnens2/crm_analysis/input_files/negative_enhancers.csv")
genes_to_check("/mnt/home/sonnens2/crm_analysis/input_files/target_genes.csv")
    



