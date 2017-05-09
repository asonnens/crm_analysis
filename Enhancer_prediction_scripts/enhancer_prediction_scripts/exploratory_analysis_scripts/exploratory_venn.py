import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from pylab import *
import scipy
import numpy
import os

input_directory = "/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/"

def make_hist(fileinput, filename):
    filetitle = filename + ".png"
    peak_score_list = []
    filein = open(fileinput)
    filein.readline()
    filein.readline()
    filein.readline()
    filein.readline()
    for line in filein:
        mylist = line.split("\t")
        if len(mylist) == 4:
            peak_score_list.append(float(mylist[3].strip()))
    plt.hist(peak_score_list)
    number_of_peaks = str(len(peak_score_list))
    print filename
    mean_val = str(numpy.mean(peak_score_list))
    median_val = str(numpy.median(peak_score_list))
    stdev_peaks = str(numpy.std(peak_score_list))
    print "number of peaks: ", number_of_peaks
    print "mean: ", mean_val
    print "sd: ", stdev_peaks
    peaks_info = number_of_peaks + " peaks \nmean " + mean_val + "\nmedian " + median_val + "\nstdev " + stdev_peaks
    axis_len = max(peak_score_list) * 0.5 
    axis_lev = len(peak_score_list) * float(stdev_peaks) * 0.7
    print axis_len, axis_lev
    plt.title(filename)
    plt.xlabel("frequency")
    plt.xlabel("peakscore")
    plt.text(axis_len, axis_lev, peaks_info, fontsize=14)
    savefig(filetitle)
    plt.close()
    return(peak_score_list)

def make_boxplot(inputscores,inputscorenames, filename):
    filetitle = filename + ".png"
    data_to_plot = inputscores
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    ax.set_xticklabels(inputscorenames)
    bp = ax.boxplot(data_to_plot)
    fig.savefig(filetitle, bbox_inches='tight')
    plt.close()

def peak_overlap(smaller_file, bigger_file):
    overlap_count = 0
    filesize1 = 0
    filesize2 = 0
    first_line = smaller_file.readline()
    if "browser" in first_line:
        smaller_file.readline()
        smaller_file.readline()
    else:
        smaller_file.seek(0)
    for line in smaller_file:
        filesize1 = filesize1 + 1
        filesize2 = 0
        each_line_overlap = 0
        second_first_line = bigger_file.readline()
        if "browser" in second_first_line:
            bigger_file.readline()
            bigger_file.readline()
        else:
            bigger_file.seek(0)
        if "chr" in line:
            file1_info = line.split("\t")
            file1_chr = file1_info[0]
            file1_coord1 = int(file1_info[1])
            file1_coord2 = int(file1_info[2])
            for line2 in bigger_file:
                filesize2 = filesize2 + 1
                file2_info = line2.split("\t")
                if file2_info[0] != file1_chr:
                    pass
                else:
                    if (((int(file2_info[1]) > file1_coord2)) or (file1_coord1 > int(file2_info[2]))):
                        pass
                    else:
                        each_line_overlap = each_line_overlap + 1
        if each_line_overlap >  0:
            overlap_count = overlap_count + 1
        bigger_file.seek(0)
    print "overlap_count, ", overlap_count
    print "filesize1 ", filesize1, ", filesize2 ", filesize2
    return(overlap_count, filesize1, filesize2)

def make_venn(direc,file1,file2,filename,text1,text2,color1,color2,color3):
    smallername = ""
    biggername = ""
    fileinput1 = direc + file1
    fileinput2 = direc + file2
    statinfo1 = os.stat(fileinput1)
    statinfo2 = os.stat(fileinput2)
    if statinfo1.st_size > statinfo2.st_size:
        smallerfile = open(fileinput2)
        biggerfile = open(fileinput1)
        smallername = text2
        biggername = text1
    else:
        smallerfile = open(fileinput1)
        biggerfile = open(fileinput2)
        smallername = text1
        biggername = text2
    overlap_info = peak_overlap(smallerfile, biggerfile)
    print overlap_info
    file1_nums = overlap_info[1] - overlap_info[0]
    file2_nums = overlap_info[2] - overlap_info[0]
    print file1_nums, file2_nums
    if biggername == text1:
        v= venn2(subsets = (file2_nums, file1_nums, overlap_info[0]), set_labels = ("", ""))
      #  v= venn2(subsets = (file2_nums, file1_nums, overlap_info[0]), set_labels = (biggername, smallername))
    else:
        v= venn2(subsets = (file1_nums, file2_nums, overlap_info[0]), set_labels = ("", ""))
       # v= venn2(subsets = (file1_nums, file2_nums, overlap_info[0]), set_labels = (smallername, biggername))
    print file2_nums, file1_nums, overlap_info[0]
    v.get_patch_by_id('10').set_color(color1)
    v.get_patch_by_id('01').set_color(color2)
    v.get_patch_by_id('11').set_color(color3)
    v.get_patch_by_id('10').set_alpha(1.0)
    v.get_patch_by_id('01').set_alpha(1.0)
    v.get_patch_by_id('11').set_alpha(0.7)
#    plt.title(filename)
#    plt.annotate(text1)
    savefig(filename)
    plt.close()


make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "bcd_2013_mel.bedgraph" , "RushlowDorsal_vs_Bicoid2013", "Rushlow Dorsal", "Bicoid 2013", '#993333', 'lime', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "bcd_macarthur.bedgraph" , "RushlowDorsal_vs_Bicoid", "Rushlow Dorsal", "Bicoid MacArthur", '#993333', 'green', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "CAD_2010_mel.bedgraph" , "RushlowDorsal_vs_Cad2010", "Rushlow Dorsal", "Caudal 2010", '#993333', 'orchid', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "cad_macarthur.bedgraph" , "RushlowDorsal_vs_Cad", "Rushlow Dorsal", "Caudal MacArthur", '#993333', 'purple', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "dl_all_macarthur.bedgraph" , "RushlowDorsal_vs_Dorsal", "Rushlow Dorsal", "Dorsal MacArthur", '#993333', 'red', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "giant_chip_chip.bedgraph" , "RushlowDorsal_vs_Gt", "Rushlow Dorsal", "Giant MacArthur", '#993333', 'orange', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "gt_2013_mel.bedgraph" , "RushlowDorsal_vs_Gt2013", "Rushlow Dorsal", "Giant 2013", '#993333', 'goldenrod', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "H3K27ac_White.bedgraph" , "RushlowDorsal_vs_WhiteH3K27", "Rushlow Dorsal", "White H3K27ac", '#993333', 'khaki', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "hairy_chip_chip.bedgraph" , "RushlowDorsal_vs_Hry", "Rushlow Dorsal", "Hairy MacArthur", '#993333', 'gray', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "hb_2013_mel.bedgraph" , "RushlowDorsal_vs_Hb2013", "Rushlow Dorsal", "Hunchback 2013", '#993333', 'brown', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "hunchback_macarthur.bedgraph" , "RushlowDorsal_vs_Hb", "Rushlow Dorsal", "Hunchback MacArthur", '#993333', 'darkslateblue', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "KNI_2010_mel.bedgraph" , "RushlowDorsal_vs_Kni2010", "Rushlow Dorsal", "Knirps 2010", '#993333', 'aquamarine', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "kr_2013_mel.bedgraph" , "RushlowDorsal_vs_Kr2013", "Rushlow Dorsal", "Kruppel 2013", '#993333', 'mediumpurple', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "RushlowDorsal_vs_KKH3K27", "Rushlow Dorsal", "KK H3K27ac", '#993333', 'gold', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "RushlowDorsal_vs_KKH3K4", "Rushlow Dorsal", "KK H3K4me1", '#993333', 'peachpuff', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "sna_all_macarthur.bedgraph" , "RushlowDorsal_vs_Snail", "Rushlow Dorsal", "Snail MacArthur", '#993333', 'dimgray', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "snail_furlong.bedgraph" , "RushlowDorsal_vs_SnailFurlong", "Rushlow Dorsal", "Snail Furlong", '#993333', 'silver', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "tw_all_macarthur.bedgraph" , "RushlowDorsal_vs_TwistMacArthur", "Rushlow Dorsal", "Twist MacArthur", '#993333', 'darkolivegreen', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "Twist_Furlong.bedgraph" , "RushlowDorsal_vs_TwistFurlong", "Rushlow Dorsal", "Twist Furlong", '#993333', 'olive', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "RushlowDorsal_vs_TwistZeit", "Rushlow Dorsal", "Twist Zeitlinger", '#993333', 'darkseagreen', 'pink')
make_venn(input_directory, "dl_rushlow_chipseq.bedgraph", "zelda_eisen_30_180min.bedgraph" , "RushlowDorsal_vs_Zelda", "Rushlow Dorsal", "Zelda", '#993333', 'peru', 'pink')

make_venn(input_directory, "dl_all_macarthur.bedgraph", "bcd_macarthur.bedgraph" , "Dorsal_MacArthur_vs_Bcd", "Dorsal MacArthur", "Bicoid MacArthur", 'tomato', 'green', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "bcd_2013_mel.bedgraph" , "Dorsal_MacArthur_vs_Bcd2013", "Dorsal MacArthur", "Bicoid 2013", 'tomato', 'lime', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "CAD_2010_mel.bedgraph" , "Dorsal_MacArthur_vs_Cad2010", "Dorsal MacArthur", "Caudal 2010", 'tomato', 'orchid', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "cad_macarthur.bedgraph" , "Dorsal_MacArthur_vs_Cad", "Dorsal MacArthur", "Cadaul MacArthur", 'tomato', 'purple', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "giant_chip_chip.bedgraph" , "Dorsal_MacArthur_vs_Gt", "Dorsal MacArthur", "Giant MacArthur", 'tomato', 'orange', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "gt_2013_mel.bedgraph" , "Dorsal_MacArthur_vs_Gt2013", "Dorsal MacArthur", "Giant 2013", 'tomato', 'goldenrod', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "Dorsal_MacArthur_vs_WhiteH3K27", "Dorsal MacArthur", "White H3K27ac", 'tomato', 'khaki', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "hairy_chip_chip.bedgraph" , "Dorsal_MacArthur_vs_Hry", "Dorsal MacArthur", "Hairy MacArthur", 'tomato', 'gray', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "Dorsal_MacArthur_vs_Hb2013", "Dorsal MacArthur", "Hunchback 2013", 'tomato', 'cadetblue', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "hunchback_macarthur.bedgraph" , "Dorsal_MacArthur_vs_Hb", "Dorsal MacArthur", "Hunchback MacArthur", 'tomato', 'darkslateblue', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "Dorsal_MacArthur_vs_Kni2010", "Dorsal MacArthur", "Knirps 2010", 'tomato', 'aquamarine', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "Dorsal_MacArthur_vs_Kr2013", "Dorsal MacArthur", "Kruppel 2013", 'tomato', 'mediumpurple', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Dorsal_MacArthur_vs_KKH3K27", "Dorsal MacArthur", "KK H3K27ac", 'tomato', 'gold', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Dorsal_MacArthur_vs_KKH3K4", "Dorsal MacArthur", "KK H3K4me1", 'tomato', 'peachpuff', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "sna_all_macarthur.bedgraph" , "Dorsal_MacArthur_vs_Snail", "Dorsal MacArthur", "Snail MacArthur", 'tomato', 'dimgray', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "snail_furlong.bedgraph" , "Dorsal_MacArthur_vs_SnailFurlong", "Dorsal MacArthur", "Snail Furlong", 'tomato', 'silver', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "tw_all_macarthur.bedgraph" , "Dorsal_MacArthur_vs_TwistMacArthur", "Dorsal MacArthur", "Twist MacArthur", 'tomato', 'darkolivegreen', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "Twist_Furlong.bedgraph" , "Dorsal_MacArthur_vs_TwistFurlong", "Dorsal MacArthur", "Twist Furlong", 'tomato', 'olive', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "Dorsal_MacArthur_vs_TwistZeit", "Dorsal MacArthur", "Twist Zeitlinger", 'tomato', 'darkseagreen', 'salmon')
make_venn(input_directory, "dl_all_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Dorsal_MacArthur_vs_Zelda", "Dorsal MacArthur", "Zelda", 'tomato', 'peru', 'salmon')

make_venn(input_directory, "bcd_macarthur.bedgraph", "bcd_2013_mel.bedgraph" , "MacArthurBicoid_vs_Bicoid", "MacArthur Bicoid", "Bicoid MacArthur", 'darkgreen', 'lime', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "CAD_2010_mel.bedgraph" , "MacArthurBicoid_vs_Cad2010", "MacArthur Bicoid", "Caudal 2010", 'darkgreen', 'orchid', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "cad_macarthur.bedgraph" , "MacArthurBicoid_vs_Cad", "MacArthur Bicoid", "Caudal MacArthur", 'darkgreen', 'purple', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "dl_all_macarthur.bedgraph" , "MacArthurBicoid_vs_Dorsal", "MacArthur Bicoid", "Dorsal MacArthur", 'darkgreen', 'red', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "giant_chip_chip.bedgraph" , "MacArthurBicoid_vs_Gt", "MacArthur Bicoid", "Giant MacArthur", 'darkgreen', 'orange', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "gt_2013_mel.bedgraph" , "MacArthurBicoid_vs_Gt2013", "MacArthur Bicoid", "Giant 2013", 'darkgreen', 'goldenrod', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "MacArthurBicoid_vs_WhiteH3K27", "MacArthur Bicoid", "White H3K27ac", 'darkgreen', 'khaki', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "hairy_chip_chip.bedgraph" , "MacArthurBicoid_vs_Hry", "MacArthur Bicoid", "Hairy MacArthur", 'darkgreen', 'gray', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "MacArthurBicoid_vs_Hb2013", "MacArthur Bicoid", "Hunchback 2013", 'darkgreen', 'cadetblue', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "hunchback_macarthur.bedgraph" , "MacArthurBicoid_vs_Hb", "MacArthur Bicoid", "Hunchback MacArthur", 'darkgreen', 'darkslategray', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "MacArthurBicoid_vs_Kni2010", "MacArthur Bicoid", "Knirps 2010", 'darkgreen', 'aquamarine', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "MacArthurBicoid_vs_Kr2013", "MacArthur Bicoid", "Kruppel 2013", 'darkgreen', 'mediumpurple', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "MacArthurBicoid_vs_KKH3K27", "MacArthur Bicoid", "KK H3K27ac", 'darkgreen', 'gold', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "MacArthurBicoid_vs_KKH3K4", "MacArthur Bicoid", "KK H3K4me1", 'darkgreen', 'peachpuff', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "sna_all_macarthur.bedgraph" , "MacArthurBicoid_vs_Snail", "MacArthur Bicoid", "Snail MacArthur", 'darkgreen', 'dimgray', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "snail_furlong.bedgraph" , "MacArthurBicoid_vs_SnailFurlong", "MacArthur Bicoid", "Snail Furlong", 'darkgreen', 'silver', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "tw_all_macarthur.bedgraph" , "MacArthurBicoid_vs_TwistMacArthur", "MacArthur Bicoid", "Twist MacArthur", 'darkgreen', 'darkolivegreen', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "Twist_Furlong.bedgraph" , "MacArthurBicoid_vs_TwistFurlong", "MacArthur Bicoid", "Twist Furlong", 'darkgreen', 'olive', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "MacArthurBicoid_vs_TwistZeit", "MacArthur Bicoid", "Twist Zeitlinger", 'darkgreen', 'darkseagreen', 'yellowgreen')
make_venn(input_directory, "bcd_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "MacArthurBicoid_vs_Zelda", "MacArthur Bicoid", "Zelda", 'darkgreen', 'peru', 'yellowgreen')


make_venn(input_directory, "bcd_2013_mel.bedgraph", "cad_macarthur.bedgraph" , "bicoid_2013_vs_Cad", "Bicoid 2013", "Caudal MacArthur", 'lime', 'purple', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "dl_all_macarthur.bedgraph" , "bicoid_2013_vs_Dorsal", "Bicoid 2013", "Dorsal MacArthur", 'lime', 'red', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "giant_chip_chip.bedgraph" , "bicoid_2013_vs_Gt", "Bicoid 2013", "Giant MacArthur", 'lime', 'orange', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "gt_2013_mel.bedgraph" , "bicoid_2013_vs_Gt2013", "Bicoid 2013", "Giant 2013", 'lime', 'goldenrod', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "H3K27ac_White.bedgraph" , "bicoid_2013_vs_WhiteH3K27", "Bicoid 2013", "White H3K27ac", 'lime', 'khaki', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "hairy_chip_chip.bedgraph" , "bicoid_2013_vs_Hry", "Bicoid 2013", "Hairy MacArthur", 'lime', 'gray', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "hb_2013_mel.bedgraph" , "bicoid_2013_vs_Hb2013", "Bicoid 2013", "Hunchback 2013", 'lime', 'cadetblue', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "hunchback_macarthur.bedgraph" , "bicoid_2013_vs_Hb", "Bicoid 2013", "Hunchback MacArthur", 'lime', 'darkslateblue', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "KNI_2010_mel.bedgraph" , "bicoid_2013_vs_Kni2010", "Bicoid 2013", "Knirps 2010", 'lime', 'aquamarine', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "kr_2013_mel.bedgraph" , "bicoid_2013_vs_Kr2013", "Bicoid 2013", "Kruppel 2013", 'lime', 'mediumpurple', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "bicoid_2013_vs_KKH3K27", "Bicoid 2013", "KK H3K27ac", 'lime', 'gold', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "bicoid_2013_vs_KKH3K4", "Bicoid 2013", "KK H3K4me1", 'lime', 'peachpuff', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "sna_all_macarthur.bedgraph" , "bicoid_2013_vs_Snail", "Bicoid 2013", "Snail MacArthur", 'lime', 'dimgray', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "snail_furlong.bedgraph" , "bicoid_2013_vs_SnailFurlong", "Bicoid 2013", "Snail Furlong", 'lime', 'silver', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "tw_all_macarthur.bedgraph" , "bicoid_2013_vs_TwistMacArthur", "Bicoid 2013", "Twist MacArthur", 'lime', 'darkolivegreen', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "Twist_Furlong.bedgraph" , "bicoid_2013_vs_TwistFurlong", "Bicoid 2013", "Twist Furlong", 'lime', 'olive', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "bicoid_2013_vs_TwistZeit", "Bicoid 2013", "Twist Zeitlinger", 'lime', 'darkseagreen', 'green')
make_venn(input_directory, "bcd_2013_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "bicoid_2013_vs_Zelda", "Bicoid 2013", "Zelda", 'lime', 'peru', 'green')

make_venn(input_directory, "CAD_2010_mel.bedgraph", "giant_chip_chip.bedgraph" , "CAD_2010_mel_vs_Gt", "Caudal 2010", "Giant MacArthur", 'orchid', 'orange', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "gt_2013_mel.bedgraph" , "CAD_2010_mel_vs_Gt2013", "Caudal 2010", "Giant 2013", 'orchid', 'goldenrod', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "H3K27ac_White.bedgraph" , "CAD_2010_mel_vs_WhiteH3K27", "Caudal 2010", "White H3K27ac", 'orchid', 'khaki', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "hairy_chip_chip.bedgraph" , "CAD_2010_mel_vs_Hry", "Caudal 2010", "Hairy MacArthur", 'orchid', 'gray', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "hb_2013_mel.bedgraph" , "CAD_2010_mel_vs_Hb2013", "Caudal 2010", "Hunchback 2013", 'orchid', 'cadetblue', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "hunchback_macarthur.bedgraph" , "CAD_2010_mel_vs_Hb", "Caudal 2010", "Hunchback MacArthur", 'orchid', 'darkslateblue', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "KNI_2010_mel.bedgraph" , "CAD_2010_mel_vs_Kni2010", "Caudal 2010", "Knirps 2010", 'orchid', 'aquamarine', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "kr_2013_mel.bedgraph" , "CAD_2010_mel_vs_Kr2013", "Caudal 2010", "Kruppel 2013", 'orchid', 'mediumpurple', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "CAD_2010_mel_vs_KKH3K27", "Caudal 2010", "KK H3K27ac", 'orchid', 'gold', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "CAD_2010_mel_vs_KKH3K4", "Caudal 2010", "KK H3K4me1", 'orchid', 'peachpuff', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "sna_all_macarthur.bedgraph" , "CAD_2010_mel_vs_Snail", "Caudal 2010", "Snail MacArthur", 'orchid', 'dimgray', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "snail_furlong.bedgraph" , "CAD_2010_mel_vs_SnailFurlong", "Caudal 2010", "Snail Furlong", 'orchid', 'silver', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "tw_all_macarthur.bedgraph" , "CAD_2010_mel_vs_TwistMacArthur", "Caudal 2010", "Twist MacArthur", 'orchid', 'darkolivegreen', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "Twist_Furlong.bedgraph" , "CAD_2010_mel_vs_TwistFurlong", "Caudal 2010", "Twist Furlong", 'orchid', 'olive', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "CAD_2010_mel_vs_TwistZeit", "Caudal 2010", "Twist Zeitlinger", 'orchid', 'darkseagreen', 'pink')
make_venn(input_directory, "CAD_2010_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "CAD_2010_mel_vs_Zelda", "Caudal 2010", "Zelda", 'orchid', 'peru', 'pink')

make_venn(input_directory, "cad_macarthur.bedgraph", "CAD_2010_mel.bedgraph" , "cad_macarthur_vs_Cad", "Caudal MacArthur", "Caudal MacArthur", 'purple', 'orchid', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "giant_chip_chip.bedgraph" , "cad_macarthur_vs_Gt", "Caudal MacArthur", "Giant MacArthur", 'purple', 'orange', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "gt_2013_mel.bedgraph" , "cad_macarthur_vs_Gt2013", "Caudal MacArthur", "Giant 2013", 'purple', 'goldenrod', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "cad_macarthur_vs_WhiteH3K27", "Caudal MacArthur", "White H3K27ac", 'purple', 'khaki', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "hairy_chip_chip.bedgraph" , "cad_macarthur_vs_Hry", "Caudal MacArthur", "Hairy MacArthur", 'purple', 'gray', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "cad_macarthur_vs_Hb2013", "Caudal MacArthur", "Hunchback 2013", 'purple', 'cadetblue', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "hunchback_macarthur.bedgraph" , "cad_macarthur_vs_Hb", "Caudal MacArthur", "Hunchback MacArthur", 'purple', 'darkslateblue', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "cad_macarthur_vs_Kni2010", "Caudal MacArthur", "Knirps 2010", 'purple', 'aquamarine', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "cad_macarthur_vs_Kr2013", "Caudal MacArthur", "Kruppel 2013", 'purple', 'mediumpurple', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "cad_macarthur_vs_KKH3K27", "Caudal MacArthur", "KK H3K27ac", 'purple', 'gold', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "cad_macarthur_vs_KKH3K4", "Caudal MacArthur", "KK H3K4me1", 'purple', 'peachpuff', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "sna_all_macarthur.bedgraph" , "cad_macarthur_vs_Snail", "Caudal MacArthur", "Snail MacArthur", 'purple', 'dimgray', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "snail_furlong.bedgraph" , "cad_macarthur_vs_SnailFurlong", "Caudal MacArthur", "Snail Furlong", 'purple', 'silver', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "tw_all_macarthur.bedgraph" , "cad_macarthur_vs_TwistMacArthur", "Caudal MacArthur", "Twist MacArthur", 'purple', 'darkolivegreen', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "Twist_Furlong.bedgraph" , "cad_macarthur_vs_TwistFurlong", "Caudal MacArthur", "Twist Furlong", 'purple', 'olive', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "cad_macarthur_vs_TwistZeit", "Caudal MacArthur", "Twist Zeitlinger", 'purple', 'darkseagreen', 'lightcoral')
make_venn(input_directory, "cad_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "cad_macarthur_vs_Zelda", "Caudal MacArthur", "Zelda", 'purple', 'peru', 'lightcoral')

make_venn(input_directory, "sna_all_macarthur.bedgraph", "giant_chip_chip.bedgraph" , "sna_vs_Gt", "Snail MacArthur", "Giant MacArthur", 'dimgray', 'orange', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "gt_2013_mel.bedgraph" , "sna_vs_Gt2013", "Snail MacArthur", "Giant 2013", 'dimgray', 'goldenrod', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "sna_vs_WhiteH3K27", "Snail MacArthur", "White H3K27ac", 'dimgray', 'khaki', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "hairy_chip_chip.bedgraph" , "sna_vs_Hry", "Snail MacArthur", "Hairy MacArthur", 'dimgray', 'gray', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "sna_vs_Hb2013", "Snail MacArthur", "Hunchback 2013", 'dimgray', 'cadetblue', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "hunchback_macarthur.bedgraph" , "sna_vs_Hb", "Snail MacArthur", "Hunchback MacArthur", 'dimgray', 'darkslateblue', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "sna_vs_Kni2010", "Snail MacArthur", "Knirps 2010", 'dimgray', 'aquamarine', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "sna_vs_Kr2013", "Snail MacArthur", "Kruppel 2013", 'dimgray', 'mediumpurple', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "sna_vs_KKH3K27", "Snail MacArthur", "KK H3K27ac", 'dimgray', 'gold', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "sna_vs_KKH3K4", "Snail MacArthur", "KK H3K4me1", 'dimgray', 'peachpuff', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "snail_furlong.bedgraph" , "sna_vs_SnailFurlong", "Snail MacArthur", "Snail Furlong", 'dimgray', 'silver', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "tw_all_macarthur.bedgraph" , "sna_vs_TwistMacArthur", "Snail MacArthur", "Twist MacArthur", 'dimgray', 'darkolivegreen', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "Twist_Furlong.bedgraph" , "sna_vs_TwistFurlong", "Snail MacArthur", "Twist Furlong", 'dimgray', 'olive', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "sna_vs_TwistZeit", "Snail MacArthur", "Twist Zeitlinger", 'dimgray', 'darkseagreen', 'black')
make_venn(input_directory, "sna_all_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "sna_vs_Zelda", "Snail MacArthur", "Zelda", 'dimgray', 'peru', 'black')

make_venn(input_directory, "snail_furlong.bedgraph", "giant_chip_chip.bedgraph" , "snaFurlong_vs_Gt", "Snail Furlong", "Giant MacArthur", 'silver', 'orange', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "gt_2013_mel.bedgraph" , "snaFurlong_vs_Gt2013", "Snail Furlong", "Giant 2013", 'silver', 'goldenrod', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "H3K27ac_White.bedgraph" , "snaFurlong_vs_WhiteH3K27", "Snail Furlong", "White H3K27ac", 'silver', 'khaki', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "hairy_chip_chip.bedgraph" , "snaFurlong_vs_Hry", "Snail Furlong", "Hairy MacArthur", 'silver', 'gray', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "hb_2013_mel.bedgraph" , "snaFurlong_vs_Hb2013", "Snail Furlong", "Hunchback 2013", 'silver', 'cadetblue', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "hunchback_macarthur.bedgraph" , "snaFurlong_vs_Hb", "Snail Furlong", "Hunchback MacArthur", 'silver', 'darkslateblue', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "KNI_2010_mel.bedgraph" , "snaFurlong_vs_Kni2010", "Snail Furlong", "Knirps 2010", 'silver', 'aquamarine', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "kr_2013_mel.bedgraph" , "snaFurlong_vs_Kr2013", "Snail Furlong", "Kruppel 2013", 'silver', 'mediumpurple', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "snaFurlong_vs_KKH3K27", "Snail Furlong", "KK H3K27ac", 'silver', 'gold', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "snaFurlong_vs_KKH3K4", "Snail Furlong", "KK H3K4me1", 'silver', 'peachpuff', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "tw_all_macarthur.bedgraph" , "snaFurlong_vs_TwistMacArthur", "Snail Furlong", "Twist MacArthur", 'silver', 'darkolivegreen', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "Twist_Furlong.bedgraph" , "snaFurlong_vs_TwistFurlong", "Snail Furlong", "Twist Furlong", 'silver', 'olive', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "snaFurlong_vs_TwistZeit", "Snail Furlong", "Twist Zeitlinger", 'silver', 'darkseagreen', 'darkslategray')
make_venn(input_directory, "snail_furlong.bedgraph", "zelda_eisen_30_180min.bedgraph" , "snaFurlong_vs_Zelda", "Snail Furlong", "Zelda", 'silver', 'peru', 'darkslategray')

make_venn(input_directory, "tw_all_macarthur.bedgraph", "giant_chip_chip.bedgraph" , "Twi_vs_Gt", "Twist MacArthur", "Giant MacArthur", 'darkolivegreen', 'orange', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "gt_2013_mel.bedgraph" , "Twi_vs_Gt2013", "Twist MacArthur", "Giant 2013", 'darkolivegreen', 'goldenrod', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "Twi_vs_WhiteH3K27", "Twist MacArthur", "White H3K27ac", 'darkolivegreen', 'khaki', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "hairy_chip_chip.bedgraph" , "Twi_vs_Hry", "Twist MacArthur", "Hairy MacArthur", 'darkolivegreen', 'gray', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "Twi_vs_Hb2013", "Twist MacArthur", "Hunchback 2013", 'darkolivegreen', 'cadetblue', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "hunchback_macarthur.bedgraph" , "Twi_vs_Hb", "Twist MacArthur", "Hunchback MacArthur", 'darkolivegreen', 'darkslateblue', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "Twi_vs_Kni2010", "Twist MacArthur", "Knirps 2010", 'darkolivegreen', 'aquamarine', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "Twi_vs_Kr2013", "Twist MacArthur", "Kruppel 2013", 'darkolivegreen', 'mediumpurple', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Twi_vs_KKH3K27", "Twist MacArthur", "KK H3K27ac", 'darkolivegreen', 'gold', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Twi_vs_KKH3K4", "Twist MacArthur", "KK H3K4me1", 'darkolivegreen', 'peachpuff', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "Twist_Furlong.bedgraph" , "Twi_vs_TwistFurlong", "Twist MacArthur", "Twist Furlong", 'darkolivegreen', 'olive', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "Twi_vs_TwistZeit", "Twist MacArthur", "Twist Zeitlinger", 'darkolivegreen', 'darkseagreen', 'palegreen')
make_venn(input_directory, "tw_all_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Twi_vs_Zelda", "Twist MacArthur", "Zelda", 'darkolivegreen', 'peru', 'palegreen')


make_venn(input_directory, "Twist_Furlong.bedgraph", "giant_chip_chip.bedgraph" , "TwiFurlong_vs_Gt", "Twist Furlong", "Giant MacArthur", 'olive', 'orange', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "gt_2013_mel.bedgraph" , "TwiFurlong_vs_Gt2013", "Twist Furlong", "Giant 2013", 'olive', 'goldenrod', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "H3K27ac_White.bedgraph" , "TwiFurlong_vs_WhiteH3K27", "Twist Furlong", "White H3K27ac", 'olive', 'khaki', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "hairy_chip_chip.bedgraph" , "TwiFurlong_vs_Hry", "Twist Furlong", "Hairy MacArthur", 'olive', 'gray', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "hb_2013_mel.bedgraph" , "TwiFurlong_vs_Hb2013", "Twist Furlong", "Hunchback 2013", 'olive', 'cadetblue', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "hunchback_macarthur.bedgraph" , "TwiFurlong_vs_Hb", "Twist Furlong", "Hunchback MacArthur", 'olive', 'darkslateblue', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "KNI_2010_mel.bedgraph" , "TwiFurlong_vs_Kni2010", "Twist Furlong", "Knirps 2010", 'olive', 'aquamarine', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "kr_2013_mel.bedgraph" , "TwiFurlong_vs_Kr2013", "Twist Furlong", "Kruppel 2013", 'olive', 'mediumpurple', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "TwiFurlong_vs_KKH3K27", "Twist Furlong", "KK H3K27ac", 'olive', 'gold', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "TwiFurlong_vs_KKH3K4", "Twist Furlong", "KK H3K4me1", 'olive', 'peachpuff', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "Zeitlinger_tw_chipseq.bedgraph" , "TwiFurlong_vs_TwistZeit", "Twist Furlong", "Twist Zeitlinger", 'olive', 'darkseagreen', 'green')
make_venn(input_directory, "Twist_Furlong.bedgraph", "zelda_eisen_30_180min.bedgraph" , "TwiFurlong_vs_Zelda", "Twist Furlong", "Zelda", 'olive', 'peru', 'green')

make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "giant_chip_chip.bedgraph" , "TwiZeitlinger_vs_Gt", "Twist Zeitlinger", "Giant MacArthur", 'darkseagreen', 'orange', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "gt_2013_mel.bedgraph" , "TwiZeitlinger_vs_Gt2013", "Twist Zeitlinger", "Giant 2013", 'darkseagreen', 'goldenrod', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "H3K27ac_White.bedgraph" , "TwiZeitlinger_vs_WhiteH3K27", "Twist Zeitlinger", "White H3K27ac", 'darkseagreen', 'khaki', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "hairy_chip_chip.bedgraph" , "TwiZeitlinger_vs_Hry", "Twist Zeitlinger", "Hairy MacArthur", 'darkseagreen', 'gray', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "hb_2013_mel.bedgraph" , "TwiZeitlinger_vs_Hb2013", "Twist Zeitlinger", "Hunchback 2013", 'darkseagreen', 'cadetblue', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "hunchback_macarthur.bedgraph" , "TwiZeitlinger_vs_Hb", "Twist Zeitlinger", "Hunchback MacArthur", 'darkseagreen', 'darkslateblue', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "KNI_2010_mel.bedgraph" , "TwiZeitlinger_vs_Kni2010", "Twist Zeitlinger", "Knirps 2010", 'darkseagreen', 'aquamarine', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "kr_2013_mel.bedgraph" , "TwiZeitlinger_vs_Kr2013", "Twist Zeitlinger", "Kruppel 2013", 'darkseagreen', 'mediumpurple', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "TwiZeitlinger_vs_KKH3K27", "Twist Zeitlinger", "KK H3K27ac", 'darkseagreen', 'gold', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "TwiZeitlinger_vs_KKH3K4", "Twist Zeitlinger", "KK H3K4me1", 'darkseagreen', 'peachpuff', 'skyblue')
make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "zelda_eisen_30_180min.bedgraph" , "TwiZeitlinger_vs_Zelda", "Twist Zeitlinger", "Zelda", 'darkseagreen', 'peru', 'skyblue')


make_venn(input_directory, "giant_chip_chip.bedgraph", "gt_2013_mel.bedgraph" , "Giant_MacArthur_vs_Gt2013", "Giant MacArthur", "Giant 2013", 'orange', 'goldenrod', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "H3K27ac_White.bedgraph" , "Giant_MacArthur_vs_WhiteH3K27", "Giant MacArthur", "White H3K27ac", 'orange', 'khaki', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "hairy_chip_chip.bedgraph" , "Giant_MacArthur_vs_Hry", "Giant MacArthur", "Hairy MacArthur", 'orange', 'gray', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "hb_2013_mel.bedgraph" , "Giant_MacArthur_vs_Hb2013", "Giant MacArthur", "Hunchback 2013", 'orange', 'cadetblue', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "hunchback_macarthur.bedgraph" , "Giant_MacArthur_vs_Hb", "Giant MacArthur", "Hunchback MacArthur", 'orange', 'darkslateblue', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "KNI_2010_mel.bedgraph" , "Giant_MacArthur_vs_Kni2010", "Giant MacArthur", "Knirps 2010", 'orange', 'aquamarine', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "kr_2013_mel.bedgraph" , "Giant_MacArthur_vs_Kr2013", "Giant MacArthur", "Kruppel 2013", 'orange', 'mediumpurple', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Giant_MacArthur_vs_KKH3K27", "Giant MacArthur", "KK H3K27ac", 'orange', 'gold', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Giant_MacArthur_vs_KKH3K4", "Giant MacArthur", "KK H3K4me1", 'orange', 'peachpuff', 'salmon')
make_venn(input_directory, "giant_chip_chip.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Giant_MacArthur_vs_Zelda", "Giant MacArthur", "Zelda", 'orange', 'peru', 'salmon')


make_venn(input_directory, "gt_2013_mel.bedgraph", "H3K27ac_White.bedgraph" , "Gt2013_vs_WhiteH3K27", "Giant 2013", "White H3K27ac", 'goldenrod', 'khaki', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "hairy_chip_chip.bedgraph" , "Gt2013_vs_Hry", "Giant 2013", "Hairy MacArthur", 'goldenrod', 'gray', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "hb_2013_mel.bedgraph" , "Gt2013_vs_Hb2013", "Giant 2013", "Hunchback 2013", 'goldenrod', 'cadetblue', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "hunchback_macarthur.bedgraph" , "Gt2013_vs_Hb", "Giant 2013", "Hunchback MacArthur", 'goldenrod', 'darkslateblue', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "KNI_2010_mel.bedgraph" , "Gt2013_vs_Kni2010", "Giant 2013", "Knirps 2010", 'goldenrod', 'aquamarine', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "kr_2013_mel.bedgraph" , "Gt2013_vs_Kr2013", "Giant 2013", "Kruppel 2013", 'goldenrod', 'mediumpurple', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Gt2013_vs_KKH3K27", "Giant 2013", "KK H3K27ac", 'goldenrod', 'gold', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Gt2013_vs_KKH3K4", "Giant 2013", "KK H3K4me1", 'goldenrod', 'peachpuff', 'yellow')
make_venn(input_directory, "gt_2013_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Gt2013_vs_Zelda", "Giant 2013", "Zelda", 'goldenrod', 'peru', 'yellow')

make_venn(input_directory, "hairy_chip_chip.bedgraph", "H3K27ac_White.bedgraph" , "Hairy_MacArthur_vs_WhiteH3K27", "Hairy MacArthur", "White H3K27ac", 'sienna', 'khaki', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "hb_2013_mel.bedgraph" , "Hairy_MacArthur_vs_Hb2013", "Hairy_MacArthur", "Hunchback 2013", 'sienna', 'cadetblue', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "hunchback_macarthur.bedgraph" , "Hairy_MacArthur_vs_Hb", "Hairy_MacArthur", "Hunchback MacArthur", 'sienna', 'darkslateblue', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "KNI_2010_mel.bedgraph" , "Hairy_MacArthur_vs_Kni2010", "Hairy_MacArthur", "Knirps 2010", 'sienna', 'aquamarine', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "kr_2013_mel.bedgraph" , "Hairy_MacArthur_vs_Kr2013", "Hairy_MacArthur", "Kruppel 2013", 'sienna', 'mediumpurple', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Hairy_MacArthur_vs_KKH3K27", "Hairy_MacArthur", "KK H3K27ac", 'sienna', 'gold', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Hairy_MacArthur_vs_KKH3K4", "Hairy_MacArthur", "KK H3K4me1", 'sienna', 'peachpuff', 'tan')
make_venn(input_directory, "hairy_chip_chip.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Hairy_MacArthur_vs_Zelda", "Hairy MacArthur", "Zelda", 'sienna', 'peru', 'tan')

make_venn(input_directory, "hunchback_macarthur.bedgraph", "H3K27ac_White.bedgraph" , "Hunchback_MacArthur_vs_WhiteH3K27", "Hunchback MacArthur", "White H3K27ac", 'darkslateblue', 'khaki', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "hb_2013_mel.bedgraph" , "Hunchback_MacArthur_vs_Hb2013", "Hunchback_MacArthur", "Hunchback 2013", 'darkslateblue', 'cadetblue', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "KNI_2010_mel.bedgraph" , "Hunchback_MacArthur_vs_Kni2010", "Hunchback_MacArthur", "Knirps 2010", 'darkslateblue', 'aquamarine', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "kr_2013_mel.bedgraph" , "Hunchback_MacArthur_vs_Kr2013", "Hunchback_MacArthur", "Kruppel 2013", 'darkslateblue', 'mediumpurple', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Hunchback_MacArthur_vs_KKH3K27", "Hunchback_MacArthur", "KK H3K27ac", 'darkslateblue', 'gold', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Hunchback_MacArthur_vs_KKH3K4", "Hunchback_MacArthur", "KK H3K4me1", 'darkslateblue', 'peachpuff', 'blue')
make_venn(input_directory, "hunchback_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Hunchback_MacArthur_vs_Zelda", "Hunchback MacArthur", "Zelda", 'darkslateblue', 'peru', 'blue')

make_venn(input_directory, "hb_2013_mel.bedgraph", "H3K27ac_White.bedgraph" , "Hb_2013_vs_WhiteH3K27", "Hunchback 2013", "White H3K27ac", 'cadetblue', 'khaki', 'springgreen')
make_venn(input_directory, "hb_2013_mel.bedgraph", "KNI_2010_mel.bedgraph" , "Hb_2013_vs_Kni2010", "Hb_2013", "Knirps 2010", 'cadetblue', 'aquamarine', 'springgreen')
make_venn(input_directory, "hb_2013_mel.bedgraph", "kr_2013_mel.bedgraph" , "Hb_2013_vs_Kr2013", "Hb_2013", "Kruppel 2013", 'cadetblue', 'mediumpurple', 'springgreen')
make_venn(input_directory, "hb_2013_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Hb_2013_vs_KKH3K27", "Hb_2013", "KK H3K27ac", 'cadetblue', 'gold', 'springgreen')
make_venn(input_directory, "hb_2013_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Hb_2013_vs_KKH3K4", "Hb_2013", "KK H3K4me1", 'cadetblue', 'peachpuff', 'springgreen')
make_venn(input_directory, "hb_2013_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Hb_2013_vs_Zelda", "Hb 2013", "Zelda", 'cadetblue', 'peru', 'springgreen')

make_venn(input_directory, "KNI_2010_mel.bedgraph", "H3K27ac_White.bedgraph" , "Kni_2010_vs_WhiteH3K27", "Kni_2010", "White H3K27ac", 'aquamarine', 'khaki', 'seagreen')
make_venn(input_directory, "KNI_2010_mel.bedgraph", "kr_2013_mel.bedgraph" , "Kni_2010_vs_Kr2013", "Kni_2010", "Kruppel 2013", 'aquamarine', 'mediumpurple', 'seagreen')
make_venn(input_directory, "KNI_2010_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Kni_2010_vs_KKH3K27", "Kni_2010", "KK H3K27ac", 'aquamarine', 'gold', 'seagreen')
make_venn(input_directory, "KNI_2010_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Kni_2010_vs_KKH3K4", "Kni_2010", "KK H3K4me1", 'aquamarine', 'peachpuff', 'seagreen')
make_venn(input_directory, "KNI_2010_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Kni_2010_vs_Zelda", "Hb 2013", "Zelda", 'aquamarine', 'peru', 'seagreen')

make_venn(input_directory, "kr_2013_mel.bedgraph", "H3K27ac_White.bedgraph" , "Kr_2013_vs_WhiteH3K27", "Kr_2013", "White H3K27ac", 'mediumpurple', 'khaki', 'crimson')
make_venn(input_directory, "kr_2013_mel.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Kr_2013_vs_KKH3K27", "Kr_2013", "KK H3K27ac", 'mediumpurple', 'gold', 'crimson')
make_venn(input_directory, "kr_2013_mel.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Kr_2013_vs_KKH3K4", "Kr_2013", "KK H3K4me1", 'mediumpurple', 'peachpuff', 'crimson')
make_venn(input_directory, "kr_2013_mel.bedgraph", "zelda_eisen_30_180min.bedgraph" , "Kr_2013_vs_Zelda", "Hb 2013", "Zelda", 'mediumpurple', 'peru', 'crimson')

make_venn(input_directory, "zelda_eisen_30_180min.bedgraph", "H3K27ac_White.bedgraph" , "Zelda_vs_WhiteH3K27", "Zelda", "White H3K27ac", 'peru', 'khaki', 'saddlebrown')
make_venn(input_directory, "zelda_eisen_30_180min.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "Zelda_vs_KKH3K27", "Zelda", "KK H3K27ac", 'peru', 'gold', 'saddlebrown')
make_venn(input_directory, "zelda_eisen_30_180min.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "Zelda_vs_KKH3K4", "Zelda", "KK H3K4me1", 'peru', 'peachpuff', 'saddlebrown')

make_venn(input_directory, "H3K27ac_White.bedgraph", "Kurtulus_H3K27ac.bedgraph" , "WhiteH3K27ac_vs_KKH3K27", "WhiteH3K27ac", "KK H3K27ac", 'khaki', 'gold', 'salmon')
make_venn(input_directory, "H3K27ac_White.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "WhiteH3K27ac_vs_KKH3K4", "WhiteH3K27ac", "KK H3K4me1", 'khaki', 'peachpuff', 'salmon')

make_venn(input_directory, "Kurtulus_H3K27ac.bedgraph", "Kurtulus_H3K4me1.bedgraph" , "KKH3K27ac_vs_KKH3K4", "WhiteH3K27ac", "KK H3K4me1", 'gold', 'peachpuff', 'red')
