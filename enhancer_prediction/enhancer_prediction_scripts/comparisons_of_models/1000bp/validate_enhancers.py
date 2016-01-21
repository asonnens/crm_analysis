
#function 1- convert output to bedgraphs
#function 2- enter annotated enhancers as a dict/list
#function 3- compare each dataset with annotated enhancers and output a dict or list of success
#function 4- scatterplot of results

import matplotlib
matplotlib.use('Agg')
import glob
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import scipy
import numpy



def function1(filelist):
    for i in filelist:
        infile = open(i)
        outname = i.replace(".tsv", ".bedgraph")
        outfile = open(outname, "w")
        outfile.write("browser hide all" + "\r\n")
        outfile.write("browser pack refGene encodeRegions" + "\r\n")
        outfile.write("browser full altGraph" + "\r\n")
        outfile.write('track type=bedGraph name="' + i.replace(".tsv", "") + '" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20' + "\r\n")
        for line in infile:
            line = line.replace('"', '')
            linelist = line.split("\t")
            each_line = linelist[1].split("_")
            if len(each_line) > 1:
                outstring = each_line[1] + "\t" + each_line[2] + "\t" + each_line[3]
                if linelist[3].strip() == "on":
                    outstring = outstring + "\t" + "5\r\n"
                elif linelist[3].strip() == "off":
                    outstring = outstring + "\t" + "-5\r\n"
                else:
                    outstring = "pass"
                if outstring != "pass":
                    outfile.write(outstring) 
        infile.close()
        outfile.close()

def function2(answerfile):
    answer_list = []
    answers = open(answerfile)
    for line in answers:
        if line[0:3] == "chr":
            linelist = line.split("\t")
            thisline = [linelist[0], int(linelist[1]), int(linelist[2])]
            answer_list.append(thisline) 
    answers.close()
    return(answer_list)

def function3(filelist, answerfile, annotated_list):
    count = 0
    for i in filelist:
        if i != answerfile:
            print "working"
            count = count + 1
            color_me = "black"
            if "unbalanced" in i:
                color_me = "black"
            else:
                color_me = "blue"
            ifile = open(i)
            total_ann = 0
            total_match = 0
            total_mismatch = 0
            for line in ifile:
                if (line[0:3] == "chr" in line) and ("-5" not in line):
                    match = 0
                    mismatch = 0
                    thisline = line.split("\t")
                    chr = thisline[0]
                    coord1 = int(thisline[1])
                    coord2 = int(thisline[2])
                    total_ann = 0
                    for ann in annotated_list:
                        total_ann = total_ann + ann[2] - ann[1]
                        if (chr != ann[0]) or (ann[1] > coord2) or (coord1 > ann[2]):
                                pass
                        else:
                            match = ann[2] - ann[1]
                            if (coord2 >= (ann[2] + 50)):
                                mismatch = mismatch + (coord2 - (ann[2] + 50))
                            if (ann[2] > coord2):
                                match = match - (ann[2] - coord2)
                            if (ann[1] >= (coord1 + 50)):
                                mismatch = mismatch + (ann[1] - (coord1 + 50))
                            if (coord1 > ann[1]):
                                match = match - (coord1 - ann[1])                            
                    if match == 0:
                        total_mismatch = total_mismatch + (coord2- coord1)
                    else:
                        total_mismatch = total_mismatch + mismatch
                        total_match = total_match + match
            thisfile_accuracy = [float(total_match)/float(total_ann), float(total_mismatch)/float(total_ann), i]
            plt.scatter(thisfile_accuracy[1],thisfile_accuracy[0], color = color_me)
            x = thisfile_accuracy[1]
            y = thisfile_accuracy[0]
            text_x = thisfile_accuracy[1] + -.1
            text_y = thisfile_accuracy[0] - 0.5
            #plt.annotate(i.replace(".bedgraph", ""), 
            #               xy = (x, y), xytext = (count, count),
            #               textcoords = 'offset points', ha = 'right', va = 'bottom',
            #               bbox = dict(boxstyle = 'round,pad=0.1', fc = color_me, alpha = 0.2),
            #               arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'), size = 8)  
        plt.ylim([-1,1])
        plt.xlim([-1,5]) 
        plt.xlabel("False Positive")
        plt.ylabel("True Positive") 
        savefig("temp")
                    



filelist1 = glob.glob("*.tsv")
filelist2 = glob.glob("*.bedgraph")
answer_file = "Annotated_enhancers.bedgraph"

function1(filelist1)
annotated = function2(answer_file)
function3(filelist2, answer_file, annotated)
