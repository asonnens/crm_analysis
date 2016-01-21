
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
    FSna_list = []
    FTwi_list = []
    Twi09_list = []
    Dl15_list = []
    Dl09_list = []
    Zld_list = []
    all_data_list = []
    FSna_list2 = []
    FTwi_list2 = []
    Twi09_list2 = []
    Dl15_list2 = []
    Dl09_list2 = []
    Zld_list2 = []
    all_data_list2 = []
    for i in filelist:
        if i != answerfile:
            print "working"
            count = count + 1
            color_me = "purple"
            ifile = open(i)
            print i
            total_ann = 0
            total_match = 0
            total_mismatch = 0
            for line in ifile:
                print line
                if (line[0:3] == "chr") and ("-5" not in line):
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
                            if (coord2 >= (ann[2] + 100)):
                                mismatch = mismatch + (coord2 - (ann[2] + 100))
                            if (ann[2] > coord2):
                                match = match - (ann[2] - coord2)
                            if (ann[1] >= (coord1 + 100)):
                                mismatch = mismatch + (ann[1] - (coord1 + 100))
                            if (coord1 > ann[1]):
                                match = match - (coord1 - ann[1])                            
                    if match == 0:
                        total_mismatch = total_mismatch + (coord2- coord1)
                    else:
                        total_mismatch = total_mismatch + mismatch
                        total_match = total_match + match
            if total_match > 0:
                thisfile_accuracy = float(total_match)/float(total_ann)
                thisfile_inaccuracy = float(total_mismatch)/float(total_ann)
            else:
                thisfile_accuracy = float(0)
                thisfile_inaccuracy = float(0)
            if "FSna" in i:
                color_me = "darkgreen"
                FSna_list.append(thisfile_accuracy)
                FSna_list2.append(thisfile_inaccuracy)
                print "Fsna"
            elif "FTwi" in i:
                color_me = "darkblue"
                FTwi_list.append(thisfile_accuracy)
                FTwi_list2.append(thisfile_inaccuracy)
                print "FTwi"
            elif "Twi09" in i:
                color_me = "cadetblue"
                Twi09_list.append(thisfile_accuracy)
                Twi09_list2.append(thisfile_inaccuracy)
                print "Twi09"
            elif "Dl15" in i:
                color_me = "darkred"
                Dl15_list.append(thisfile_accuracy)
                Dl15_list2.append(thisfile_inaccuracy)
                print "Dl15"
            elif "Dl09" in i:
                color_me = "orange"
                Dl09_list.append(thisfile_accuracy)
                Dl09_list2.append(thisfile_inaccuracy)
                print "Dl09"
            elif "Zld" in i:
                color_me = "yellow"
                Zld_list.append(thisfile_accuracy)
                Zld_list2.append(thisfile_inaccuracy)
                print "zelda"
            else:
                print "all data"
                all_data_list.append(thisfile_accuracy)
                all_data_list2.append(thisfile_inaccuracy)
                color_me = "black"
            
            #plt.scatter(thisfile_accuracy[1],thisfile_accuracy[0], color = color_me, alpha = 1)
            #x = thisfile_accuracy[1]
            #y = thisfile_accuracy[0]
            #text_x = thisfile_accuracy[1] + -.1
            #text_y = thisfile_accuracy[0] - 0.5
            #plt.annotate(i.replace(".bedgraph", ""), 
            #               xy = (x, y), xytext = (count, count),
            #               textcoords = 'offset points', ha = 'right', va = 'bottom',
            #               bbox = dict(boxstyle = 'round,pad=0.1', fc = color_me, alpha = 0.2),
            #               arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'), size = 8)  
        
    list_of_models = [all_data_list, Dl15_list, Dl09_list, FSna_list, FTwi_list, Twi09_list, Zld_list]
    for i in list_of_models:
        print i
    list_of_means = np.asarray([np.mean(np.asarray(i)) for i in list_of_models])
    print list_of_means
    list_of_sd = np.asarray([np.std(np.asarray(i), ddof = 2) for i in list_of_models])
    list_of_models2 = [all_data_list2, Dl15_list2, Dl09_list2, FSna_list2, FTwi_list2, Twi09_list2, Zld_list2]
    list_of_means2 = np.asarray([np.mean(np.asarray(i)) for i in list_of_models2])
    list_of_sd2 = np.asarray([np.std(np.asarray(i), ddof =2) for i in list_of_models2])
    N = len(list_of_means)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, list_of_means, width, color='r', yerr=list_of_sd)
    rects2 = ax.bar(ind + width, list_of_means2, width, color='y', yerr=list_of_sd2)
    ax.set_ylabel('Rate')
    ax.set_title('Percentage Match and Mismatch by model- all')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('all_data', 'Dl15', 'Dl09', 'Sna14', 'Twi14', 'Twi09', 'Zld11'))
    ax.legend((rects1[0], rects2[0]), ('Match', 'Mismatch'))

    savefig("comparing_unbalanced_models_bargraph")
                    



filelist1 = glob.glob("*.tsv")

answer_file = "Annotated_enhancers.bedgraph"

function1(filelist1)
print "function 1 complete"
annotated = function2(answer_file)
print "function 2 complete"
filelist2 = glob.glob("*.bedgraph")
function3(filelist2, answer_file, annotated)
