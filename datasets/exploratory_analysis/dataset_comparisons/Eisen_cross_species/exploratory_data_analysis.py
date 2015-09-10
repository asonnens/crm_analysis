import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from pylab import *
import scipy
import numpy
import os


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
    smaller_file.readline()
    smaller_file.readline()
    smaller_file.readline()
    for line in smaller_file:
        filesize1 = filesize1 + 1
        filesize2 = 0
        each_line_overlap = 0
        bigger_file.readline()
        bigger_file.readline()
        bigger_file.readline()
        bigger_file.readline()
        file1_info = line.split("\t")
        file1_chr = file1_info[0]
        if len(file1_info) ==4:
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

def make_venn(fileinput1, fileinput2, filename):
    statinfo1 = os.stat(fileinput1)
    statinfo2 = os.stat(fileinput2)
    if statinfo1.st_size > statinfo2.st_size:
        smallerfile = open(fileinput2)
        biggerfile = open(fileinput1)
    else:
        smallerfile = open(fileinput1)
        biggerfile = open(fileinput2)
    overlap_info = peak_overlap(smallerfile, biggerfile)
    print overlap_info
    file1_nums = overlap_info[1] - overlap_info[0]
    file2_nums = overlap_info[2] - overlap_info[0]
    print file1_nums, file2_nums
    if statinfo1.st_size > statinfo2.st_size:
        venn2(subsets = (file2_nums, file1_nums, overlap_info[0]))
        print file2_nums, file1_nums, overlap_info[0]
        plt.title(filename)
        savefig(filename)
        plt.close()
    else:
        venn2(subsets = (file1_nums, file2_nums, overlap_info[0]))
        plt.title(filename)
        savefig(filename)
        plt.close()

    
bcd_mel_2010 = make_hist("BCD_2010_mel.bedgraph", "bicoid_2010mel_eisen")
bcd_yak_2010 = make_hist("BCD_2010_yak.bedgraph", "bicoid_2010yak_eisen")
bcd_mel_2013 = make_hist("bcd_2013_mel.bedgraph", "bicoid_2013mel_eisen")
bcd_yak_2013 = make_hist("bcd_2013_yak.bedgraph", "bicoid_2013yak_eisen")
bcd_pse_2013 = make_hist("bcd_2013_pse.bedgraph", "bicoid_2013pse_eisen")
bcd_vir_2013 = make_hist("bcd_2013_vir.bedgraph", "bicoid_2013vir_eisen")
cad_mel_2010 = make_hist("CAD_2010_mel.bedgraph", "cad_2010mel_eisen")
gt_mel_2010 = make_hist("GT_2010_mel.bedgraph", "giant_2010mel_eisen")
kni_mel_2010 = make_hist("KNI_2010_mel.bedgraph", "knirps_2010mel_eisen")
hb1_mel_2010 = make_hist("HB1_2010_mel.bedgraph", "hb1_2010mel_eisen")
hb2_mel_2010 = make_hist("HB2_2010_mel.bedgraph", "hb2_2010mel_eisen")
kr1_mel_2010 = make_hist("KR1_2010_mel.bedgraph", "kr1_2010mel_eisen")
kr2_mel_2010 = make_hist("KR2_2010_mel.bedgraph", "kr2_2010mel_eisen")
gt_mel_2013 = make_hist("gt_2013_mel.bedgraph", "gt_2013mel_eisen")
hb_mel_2013 = make_hist("hb_2013_mel.bedgraph", "hb_2013mel_eisen")
kr_mel_2013 = make_hist("kr_2013_mel.bedgraph", "kr_2013mel_eisen")

#bcd25 = make_hist("bicoid_macarthur25.bedgraph", "bicoid_25percent_macarthur")
#bcd1 = make_hist("bcd_macarthur.bedgraph", "bicoid_1percent_macarthur")
#my_hist_data = make_hist("dl_all_macarthur.bedgraph", "dl_1percent_macarthur")
#dl_Rushlow = make_hist("dl_rushlow_chipseq.bedgraph", "Rushlow_dl")
#Twist_Furlong = make_hist("Twist_Furlong.bedgraph", "Furlong_tw")
#twi_1 = make_hist("tw_all_macarthur.bedgraph", "twist_1percent_macarthur graph4")
#twi_25 = make_hist("twist_macarthur25.bedgraph", "twist_25percent_macarthur graph4")
#zeitlinger_twist = make_hist("Zeitlinger_tw_chipseq.bedgraph", "twist_1percent_zeitlinger")
#zelda_120 = make_hist("zelda_eisen_30_120min.bedgraph", "zelda30_120percent_Eisen")
#zelda_180 = make_hist("zelda_eisen_30_180min.bedgraph", "zelda30_180_Eisen")
#White_ac = make_hist("H3K27ac_White.bedgraph", "White_H3K27ac")
#White_me = make_hist("H3K4me1_White.bedgraph", "White_H3K4me1")
#Kurt_ac = make_hist("Kurtulus_H3K27ac.bedgraph", "Kurt_H3K27ac")
#Kurt_me = make_hist("Kurtulus_H3K4me1.bedgraph", "Kurt_H3K4me1")

eisen_mel_bcd = [bcd_mel_2010, bcd_mel_2013]
eisen_mel_names = ["bcdmel2010", "bcdmel2013"]
eisen_yak_bcd = [bcd_yak_2010, bcd_yak_2013]
eisen_yak_names = ["bcdyak2010", "bcdyak2013"]
eisen_bcd_2010 = [bcd_mel_2010, bcd_yak_2010]
eisen_bcd_2010_names = ["bcdmel2010", "bcdyak2010"]
eisen_bcd_2013 = [bcd_mel_2013, bcd_yak_2013, bcd_pse_2013, bcd_vir_2013]
eisen_bcd_2013_names = ["bcdmel2013", "bcdyak13", "bcdpse13", "bcdvir13"]
eisen_mel_2010 = [cad_mel_2010, kni_mel_2010]
eisen_mel_2010_names = ["cadmel2010", "knimel2010"]
eisen_mel_2013 = [bcd_mel_2013, gt_mel_2013,hb_mel_2013,kr_mel_2013]
eisen_mel_2013_names = ["bcdmel2013","gtmel2013", "hbmel2013", "krmel2013"] 

make_boxplot(eisen_mel_2010, eisen_mel_2010_names, "eisen 2010 datasets")
make_boxplot(eisen_mel_2013, eisen_mel_2013_names, "eisen 2013 datasets")


#make_boxplot(eisen_scores, eisen_names, "eisen_boxplot1")
#make_boxplot(White_scores, White_names, "White_boxplot1")
#make_boxplot(Kurt_scores, Kurt_names, "Kurtulus boxplot1")
#make_boxplot(eisen_mel_bcd, eisen_mel_names, "eisen mel bicoid")
#make_boxplot(eisen_yak_bcd, eisen_yak_names, "eisen yak bicoid")
#make_boxplot(eisen_bcd_2010, eisen_bcd_2010_names, "eisen bicoid 2010")
#make_boxplot(eisen_bcd_2013, eisen_bcd_2013_names, "eisen bicoid 2013")
#make_venn("BCD_2010_mel.bedgraph", "BCD_2010_yak.bedgraph", "BCD 2010 mel vs yak")
#make_venn("bcd_2013_mel.bedgraph", "bcd_2013_yak.bedgraph", "BCD 2013 mel vs yak")
#make_venn("bcd_2013_mel.bedgraph", "bcd_2013_pse.bedgraph", "BCD 2013 mel vs pse")
#make_venn("bcd_2013_mel.bedgraph", "bcd_2013_vir.bedgraph", "BCD 2013 mel vs vir")
#make_venn("dl_all_macarthur.bedgraph", "dl_rushlow_chipseq.bedgraph", "Dorsal venn diagram")
#make_venn("H3K27ac_White.bedgraph", "Kurtulus_H3K27ac.bedgraph", "H3K27ac venn diagram")
#make_venn("H3K4me1_White.bedgraph", "Kurtulus_H3K4me1.bedgraph", "H3K4me1 venn diagram")
#make_venn("dl_all_macarthur.bedgraph", "H3K27ac_White.bedgraph", "Dorsal MacArthur H3K27ac")
#make_venn("dl_all_macarthur.bedgraph", "zelda_eisen_30_180min.bedgraph", "Dorsal MacArthur Zelda")
#make_venn("dl_all_macarthur.bedgraph", "sna_all_macarthur.bedgraph", "Dorsal MacArthur Snail")
#make_venn("dl_all_macarthur.bedgraph", "tw_all_macarthur.bedgraph", "Dorsal Macarthur Twist")
#make_venn("dl_rushlow_chipseq.bedgraph", "snail_furlong.bedgraph", "Dorsal Rushlow-Furlong Snail")
#make_venn("dl_rushlow_chipseq.bedgraph", "Zeitlinger_tw_chipseq.bedgraph", "Dorsal Rushlow-Zeitlinger Twist")
#make_venn("snail_furlong.bedgraph", "Zeitlinger_tw_chipseq.bedgraph", "Snail Furlong-Zeitlinger Twist")
#make_venn("sna_all_macarthur.bedgraph", "tw_all_macarthur.bedgraph", "Snail Macarthur Twist")
#make_venn("zelda_eisen_30_180min.bedgraph", "H3K27ac_White.bedgraph", "Zelda White H3K27ac")
#make_venn("zelda_eisen_30_180min.bedgraph", "H3K4me1_White.bedgraph", "Zelda White H3K4me1")
#make_venn("zelda_eisen_30_180min.bedgraph", "Kurtulus_H3K27ac.bedgraph", "Zelda Kurtulus H3K27ac")
#make_venn("zelda_eisen_30_180min.bedgraph", "Kurtulus_H3K4me1.bedgraph", "Zelda Kurtulus H3K4me1")
#make_venn("dl_rushlow_chipseq.bedgraph", "H3K27ac_White.bedgraph", "Dorsal Rushlow-White H3K27ac")
#make_venn("dl_rushlow_chipseq.bedgraph", "H3K4me1_White.bedgraph", "Dorsal Rushlow-White H3K4me1")

