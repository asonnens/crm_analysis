#This script analyzes how different bedgraph files overlap with each other
#It can output quantitative comparisons of how the results of one dataset overlap with another
#Or turn these into venn diagrams
#It also can calculate how the peak scores for overlapping peaks correlate (Pearson's Correlation)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from pylab import *
import scipy
import numpy
import os
import glob



input_directory = "/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/"

#if performing operation on all bedgraphs, use this list
my_bedgraph_list = glob.glob(input_directory + "/*.bedgraph")


#if performing operation on only one bedgraph, select here 
#(if not using, just enter "none")
print "what file do you want to check?"
my_shorter_list = raw_input()
my_shorter_list = [input_directory  + my_shorter_list]
print "high or all or both sets of peaks? (print 'high' or 'all' or 'both'):"
level_peaks = raw_input()


########################################################
#This function gets information for a specific bedgraph
#It reutrns the file information/name, and the file size
########################################################
def get_file_info(myfile):
    statinfo1 = os.stat(myfile)
    thisfile = open(myfile)
    thisfile.readline()
    thisfile.readline()
    thisfile.readline()
    annotation_line = thisfile.readline()
    annotation_line = annotation_line.replace('track type=bedGraph name="', '')
    annotation_line = annotation_line.replace('" description="BedGraph format" visibility=full color=250,150,0 altColor=0,150,250 priority=21', '')
    annotation_line = annotation_line.replace('" description="BedGraph format" visibility=full color=250,150,0 altColor=0,150,250 priority=20', '')
    annotation_line = annotation_line.replace('" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20', '')
    annotation_line = annotation_line.replace('" description="BedGraph format" visibility=full color=300,400,0 altColor=0,400,300 priority=20', '')
    annotation_line = annotation_line.strip()
    thisfile.close()
    return(annotation_line, statinfo1)

########################################################
#This function gets information for a specific bedgraph
#It reutrns the file information/name, and the file size
########################################################
def get_90th_percentile(myfile):
    statinfo1 = os.stat(myfile)
    thisfile = open(myfile)
    thisfile.readline()
    thisfile.readline()
    thisfile.readline()
    annotation_line = thisfile.readline()
    number_list = []
    for line in thisfile:
        linelist = line.split("\t")
        myval = float(linelist[3])
        number_list.append(myval)
    numbers = np.asarray(number_list)
    return_value = np.mean(number_list) + np.std(number_list)
    thisfile.close()
    return(return_value)


########################################################
#This function gets the range of scores for a particular 
#bedgraph. It returns the minimum score for that file,
#and the value for 1 sd above the mean
########################################################
def get_cutoff(filename):
    myfile = open(filename)
    mylist = []
    first_line = myfile.readline()
    if "browser" in first_line:
        myfile.readline()
        myfile.readline()
    else:
        myfile.seek(0)
    for line in myfile:
        if "chr" in line:
            file1_info = line.split("\t")
            new_num = float(file1_info[3])
            mylist.append(new_num)
    mean_val = np.mean(np.asarray(mylist))
    print "mean: ", mean_val
    std_val = np.std(np.asarray(mylist)) + mean_val
    print "Mean + 1 Sd: ", std_val
    min_val =  min(mylist)
    myfile.close()
    return(min_val, std_val)


########################################################
#This function determines how much the features in
#two bedgraphs overlap. It also determines how much
#scores correlate in corresponding regions between files
#It only looks at overlapping regions with scores above
#the predetermined threshold (minimum score for that file,
#i.e. no threshold, or 1sd above mean)
########################################################
def peak_overlap(smaller_file_name, bigger_file_name, cutoff_val, cutoff_val2):    
    smaller_file = open(smaller_file_name)
    bigger_file = open(bigger_file_name)
    print smaller_file_name, bigger_file_name
    smaller_list = []
    bigger_list = []
    overlap_count = 0
    filesize1 = 0
    filesize2 = 0
    overlap_list = []
    first_line = smaller_file.readline()
    if "browser" in first_line:
        smaller_file.readline()
        smaller_file.readline()
    else:
        smaller_file.seek(0)
    for line in smaller_file:
        filesize2 = 0
        each_line_overlap = 0
        second_first_line = bigger_file.readline()
        if "browser" in second_first_line:
            bigger_file.readline()
            bigger_file.readline()
        else:
            bigger_file.seek(0)
        if "chr" in line and float(line.split("\t")[3]) >= cutoff_val:
            filesize1 = filesize1 + 1
            file1_info = line.split("\t")
            file1_chr = file1_info[0]
            file1_coord1 = int(file1_info[1])
            file1_coord2 = int(file1_info[2])
            for line2 in bigger_file:
                if "chr" in line2 and float(line2.split("\t")[3]) >= cutoff_val2:
                    file2_info = line2.split("\t")
                    filesize2 = filesize2 + 1
                    if (file2_info[0] != file1_chr):
                        pass
                    else:
                        if (((int(file2_info[1]) > file1_coord2)) or (file1_coord1 > int(file2_info[2]))):
                            pass
                        else:
                            each_line_overlap = each_line_overlap + 1
                            smaller_list.append(float(file1_info[3].strip()))
                            bigger_list.append(float(file2_info[3].strip()))
        if each_line_overlap >  0:
            overlap_count = overlap_count + 1
            overlap_list.append([line, line2])
        bigger_file.seek(0)
    print "overlap ", overlap_count
    print "smaller list ",len(smaller_list)
    print "bigger list ", len(bigger_list)
    smaller_file.seek(0)
    bigger_file.seek(0)
    if overlap_count > 0:
        smaller_list = np.asarray(smaller_list)
        bigger_list = np.asarray(bigger_list)
        correlation = np.corrcoef(smaller_list, bigger_list)[0,1] * 100
        overlap = (float(overlap_count)/filesize1) * 100
    else:
        correlation = 0
        overlap = 0
    return(overlap, correlation)

for file in my_shorter_list:
    filename1, size1 = get_file_info(file)
    outfilename1 = filename1 + "_overlap.txt"
    outfilename2 = filename1 + "_corr.txt"
    min_val, cutoff_val = get_cutoff(file)
    if (level_peaks == "high") or (level_peaks == "both"):    
        outfilename1_high = filename1 + "_high_overlap.txt"
        outfilename2_high = filename1 + "_high_corr.txt"
        outfile_temp1_high = open(outfilename1_high, "w")
        outfile_temp1_high.write("dataset" + "\t" + filename1 +  "\r\n")
        outfile_temp2_high = open(outfilename2_high, "w")
        outfile_temp2_high.write("dataset" + "\t" + filename1 +  "\r\n")
    if (level_peaks == "all") or (level_peaks == "both"):
        outfile_temp1 = open(outfilename1, "w")
        outfile_temp1.write("dataset" + "\t" + filename1 +  "\r\n")
        outfile_temp2 = open(outfilename2, "w")
        outfile_temp2.write("dataset" + "\t" + filename1 +  "\r\n")
    for second_file in my_bedgraph_list:
        min_val2, cutoff_val2 = get_cutoff(second_file)
        second_name, size2 = get_file_info(second_file)
        print second_name
        print min_val, min_val2
        if size2 <= size1:
             if (level_peaks == "all") or (level_peaks == "both"):
                 overlap_info = peak_overlap(second_file, file, min_val2, min_val)
                 print filename1, second_name, overlap_info[0], overlap_info[1]
                 outfile_temp1.write(second_name + "\t" +  str(overlap_info[0]) + "\r\n")
                 outfile_temp2.write(second_name + "\t" +  str(overlap_info[1]) + "\r\n")
             if (level_peaks == "high") or (level_peaks == "both"):    
                 overlap_info = peak_overlap(second_file, file, cutoff_val2, cutoff_val)         
                 print filename1, second_name, overlap_info[0], overlap_info[1]
                 outfile_temp1_high.write(second_name + "\t" +  str(overlap_info[0]) + "\r\n")
                 outfile_temp2_high.write(second_name + "\t" +  str(overlap_info[1]) + "\r\n")
        else:
             if (level_peaks == "all") or (level_peaks == "both"):
                 overlap_info = peak_overlap(file, second_file, min_val, min_val2)
                 print filename1, second_name, overlap_info[0], overlap_info[1]
                 outfile_temp1.write(second_name + "\t" +  str(overlap_info[0]) + "\r\n")
                 outfile_temp2.write(second_name + "\t" +  str(overlap_info[1]) + "\r\n")
             if (level_peaks == "high") or (level_peaks == "both"):    
                 overlap_info = peak_overlap(file, second_file, cutoff_val, cutoff_val2)
                 print filename1, second_name, overlap_info[0], overlap_info[1]
                 outfile_temp1_high.write(second_name + "\t" +  str(overlap_info[0]) + "\r\n")
                 outfile_temp2_high.write(second_name + "\t" +  str(overlap_info[1]) + "\r\n")
    if (level_peaks == "all") or (level_peaks == "both"):
        outfile_temp1.close()
        outfile_temp2.close()
    if (level_peaks == "high") or (level_peaks == "both"):
        outfile_temp1_high.close()
        outfile_temp2_high.close()


