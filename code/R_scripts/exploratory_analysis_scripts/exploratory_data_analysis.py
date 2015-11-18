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
    overlap_list = []
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
            overlap_list.append([line, line2])
        bigger_file.seek(0)
    print "overlap_count, ", overlap_count
    print "filesize1 ", filesize1, ", filesize2 ", filesize2
    smaller_file.seek(0)
    bigger_file.seek(0)
    return(overlap_count, filesize1, filesize2, overlap_list)


def peak_overlap3(smallest_file,second_file,bigger_file):
    overlap1_info = peak_overlap(smallest_file, second_file)
    overlap2_info = peak_overlap(smallest_file, bigger_file)
    overlap3_info = peak_overlap(second_file, bigger_file)
    overlap_list1 = overlap1_info[3]
    overlap_list2 = overlap2_info[3]
    overlap_list3 = overlap3_info[3]
    overlap_1_1 = [i[0] for i in overlap_list1]
    overlap_1_2 = [i[1] for i in overlap_list1]
    overlap_2_1 = [i[0] for i in overlap_list2]
    overlap_2_3 = [i[1] for i in overlap_list2]
    overlap_3_2 = [i[0] for i in overlap_list3]
    overlap_3_3 = [i[1] for i in overlap_list3]
    one_2 = overlap1_info[0]
    one_3 = overlap2_info[0]
    two_3 = overlap3_info[0]
    print one_2, one_3, two_3
    one_2_3 = 0
    print "lengths ", len(overlap_1_1), len(overlap_2_1)
    for i in overlap_1_1:
        if i in overlap_2_1:
            one_2_3 = one_2_3 + 1
    one_2 = one_2 -one_2_3
    one_3 = one_3 - one_2_3
    two_3 = two_3 - one_2_3
    one_file = overlap1_info[1] - (one_2 + one_2_3 + one_3)
    two_file = overlap1_info[2] - (one_2 + one_2_3 + two_3)
    three_file = overlap2_info[2] - (one_3 + one_2_3 + two_3)
    return(one_file, two_file, one_2, three_file, one_3, two_3, one_2_3)

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
    file1_nums = overlap_info[1] - overlap_info[0]
    file2_nums = overlap_info[2] - overlap_info[0]
    print file1_nums, file2_nums
    if biggername == text1:
        v= venn2(subsets = (file2_nums, file1_nums, overlap_info[0]), set_labels = (biggername, smallername))
    else:
        v= venn2(subsets = (file1_nums, file2_nums, overlap_info[0]), set_labels = (smallername, biggername))
    print file2_nums, file1_nums, overlap_info[0]
    v.get_patch_by_id('10').set_color(color1)
    v.get_patch_by_id('01').set_color(color2)
    v.get_patch_by_id('11').set_color(color3)
    v.get_patch_by_id('10').set_alpha(1.0)
    v.get_patch_by_id('01').set_alpha(1.0)
    v.get_patch_by_id('11').set_alpha(0.7)
    plt.title(filename)
    savefig(filename)
    plt.close()

def make_venn3(direc, file1, file2, file3, filename, text1, text2, text3, color1):
    #s = (
    #    50,    # Abc
    #    50,    # aBc
    #    0,    # ABc
    #    50,    # abC
    #    10,    # AbC
    #    10,  # aBC
    #    0,    # ABC
    #)

    fileinput1 = direc + file1
    fileinput2 = direc + file2
    fileinput3 = direc + file3
    statinfo1 = os.stat(fileinput1)
    statinfo2 = os.stat(fileinput2)
    statinfo3 = os.stat(fileinput3)
    if statinfo1.st_size > statinfo2.st_size:
        if statinfo2.st_size > statinfo3.st_size:
            biggest_file = open(fileinput1)
            big_text = text1
            smallest_file = open(fileinput3)
            small_text = text3
            med_file = open(fileinput2)
            med_text = text2
        else:
            if statinfo1.st_size > statinfo3.st_size:
                biggest_file = open(fileinput1)
                big_text = text1
                smallest_file = open(fileinput2)
                small_text = text2
                med_file = open(fileinput3)
                med_text = text3
            else:
                biggest_file = open(fileinput3)
                big_text = text3
                smallest_file = open(fileinput2)
                small_text = text2
                med_file = open(fileinput1)
                med_text = text1           
    else:
        if statinfo1.st_size > statinfo3.st_size:
            biggest_file = open(fileinput2)
            big_text = text2
            smallest_file = open(fileinput3)
            small_text = text3
            med_file = open(fileinput1)
            med_text = text1
        else:
            if statinfo2.st_size > statinfo3.st_size:
                biggest_file = open(fileinput2)
                big_text = text2
                smallest_file = open(fileinput1)
                small_text = text1
                med_file = open(fileinput3)
                med_text = text3
            else:
                biggest_file = open(fileinput3)
                big_text = text3
                smallest_file = open(fileinput1)
                small_text = text1
                med_file = open(fileinput2)
                med_text = text2
    overlap_info = peak_overlap3(smallest_file, med_file, biggest_file)
    print overlap_info
    #matplotlib.rcParams.update({'font.size': 22})
    v = venn3(subsets = (overlap_info[0], overlap_info[1], overlap_info[2], overlap_info[3], overlap_info[4], overlap_info[5],overlap_info[6]), set_labels = (small_text, med_text, big_text))
    v.get_patch_by_id('100').set_color("green")
    v.get_patch_by_id('100').set_alpha(0.8)
    v.get_patch_by_id('010').set_color("cyan")
    v.get_patch_by_id('010').set_alpha(0.8)
    v.get_patch_by_id('001').set_color(color1)
    v.get_patch_by_id('001').set_alpha(1.0)
    #v = venn3(subsets=(1,1,0,1,0,0,0))
    plt.title(filename)
    savefig(filename)
    plt.close()
   

input_directory = "/mnt/home/sonnens2/crm_analysis/datasets/exploratory_analysis/dataset_comparisons/venn_diagrams/"


#make_venn(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "zelda_eisen_30_180min.bedgraph", "test1", "test1", "test1", "red", "red", "red")
#make_venn3(input_directory, "dl_rushlow_chipseq.bedgraph", "zelda_eisen_30_180min.bedgraph", "Zeitlinger_tw_chipseq.bedgraph", "Dl_Twi_Zelda", "Dl_seq", "Zelda_seq", "Twi_seq", "#993333")
#make_venn3(input_directory, "dl_rushlow_chipseq.bedgraph", "snail_furlong.bedgraph", "Twist_Furlong.bedgraph", "Dl_Sna_Twi", "Dl_seq", "Sna_2014_chip", "Twi_2014_chip", "#993333")
#make_venn3(input_directory, "DV_TF.bedgraph", "zelda_eisen_30_180min.bedgraph", "Bcd_Cad.bedgraph", "DV_Zeld_BcdCad", "Dl_Sna_Twi", "Zelda_seq", "Bcd_Cad", "#993333")
#make_venn3(input_directory, "DV_TF.bedgraph", "zelda_eisen_30_180min.bedgraph", "KK_histones.bedgraph", "DV_Zld_histones", "Dl_Sna_Twi", "Zelda_seq", "Histone marks", "#993333")
#make_venn3(input_directory, "KK_histones.bedgraph", "zelda_eisen_30_180min.bedgraph", "Bcd_Cad.bedgraph", "Histones_Zld_BcdCad", "Histone marks", "Zelda_seq", "Bcd_Cad", "#993333")

#make_venn3(input_directory, "bcd_macarthur.bedgraph", "Stark_expression", "Stark_off", "Bcd 2009", "Express", "Off", "Bcd 2009", "darkgreen")
#make_venn3(input_directory, "bcd_2013_mel.bedgraph", "Stark_expression", "Stark_off", "Bcd 2013", "Express", "Off", "Bcd 2013", "lime")
#make_venn3(input_directory, "cad_macarthur.bedgraph", "Stark_expression", "Stark_off", "Cad 2009", "Express", "Off", "Cad 2009", "purple")
#make_venn3(input_directory, "CAD_2010_mel.bedgraph", "Stark_expression", "Stark_off", "Cad 2010", "Express", "Off", "Cad 2010", "orchid")
#make_venn3(input_directory, "dl_all_macarthur.bedgraph", "Stark_expression", "Stark_off", "Dl 2009", "Express", "Off", "Dl 2009", "tomato")
#make_venn3(input_directory, "giant_chip_chip.bedgraph", "Stark_expression", "Stark_off", "Gt 2009", "Express", "Off", "Gt 2009", "orange")
#make_venn3(input_directory, "gt_2013_mel.bedgraph", "Stark_expression", "Stark_off", "Gt 2013", "Express", "Off", "Gt 2013", "goldenrod")
#make_venn3(input_directory, "H3K27ac_White.bedgraph", "Stark_expression", "Stark_off", "White H3K27ac", "Express", "Off", "W_H3K27ac", "khaki")
#make_venn3(input_directory, "hairy_chip_chip.bedgraph", "Stark_expression", "Stark_off", "Hry 2009", "Express", "Off", "Hry 2009", "sienna")
#make_venn3(input_directory, "hb_2013_mel.bedgraph", "Stark_expression", "Stark_off", "Hb 2013", "Express", "Off", "Hb 2013", "cadetblue")
#make_venn3(input_directory, "hunchback_macarthur.bedgraph", "Stark_expression", "Stark_off", "Hb 2009", "Express", "Off", "Hb 2009", "darkslategray")
#make_venn3(input_directory, "KNI_2010_mel.bedgraph", "Stark_expression", "Stark_off", "Kni 2010", "Express", "Off", "Kni 2010", "aquamarine")
#make_venn3(input_directory, "KNI_2010_mel.bedgraph", "Stark_expression", "Stark_off", "Kni 2010", "Express", "Off", "Kni 2010", "aquamarine")
#make_venn3(input_directory, "kr_2013_mel.bedgraph", "Stark_expression", "Stark_off", "Kr 2013", "Express", "Off", "Kr 2013", "mediumpurple")
#make_venn3(input_directory, "Kurtulus_H3K27ac.bedgraph", "Stark_expression", "Stark_off", "KK H3K27ac", "Express", "Off", "KK H3K27ac", "gold")
#make_venn3(input_directory, "Kurtulus_H3K4me1.bedgraph", "Stark_expression", "Stark_off", "KK H3K4me1", "Express", "Off", "KK H3K4me1", "peachpuff")
#make_venn3(input_directory, "sna_all_macarthur.bedgraph", "Stark_expression", "Stark_off", "Sna 2009", "Express", "Off", "Sna 2009", "dimgray")
#make_venn3(input_directory, "snail_furlong.bedgraph", "Stark_expression", "Stark_off", "Snail Furlong", "Express", "Off", "Snail Furlong", "silver")
#make_venn3(input_directory, "tw_all_macarthur.bedgraph", "Stark_expression", "Stark_off", "Twi 2009", "Express", "Off", "Twi 2009", "darkolivegreen")
#make_venn3(input_directory, "Twist_Furlong.bedgraph", "Stark_expression", "Stark_off", "Furlong Twi", "Express", "Off", "Furlong Twi", "olive")
#make_venn3(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "Stark_expression", "Stark_off", "Zeit Twi", "Express", "Off", "Zeit Twi", "darkseagreen")
#make_venn3(input_directory, "zelda_eisen_30_180min.bedgraph", "Stark_expression", "Stark_off", "Zelda", "Express", "Off", "Zelda", "peru")



#make_venn3(input_directory, "Zeitlinger_tw_chipseq.bedgraph", "Twist_Furlong.bedgraph", "tw_all_macarthur.bedgraph", "TwiZeit", "Twi2009", "TwistFurlong", "Twist_Something","peru")

#make_venn3(input_directory, "p300_White.bedgraph", "Stark_expression", "Stark_off", "p300", "Express", "Off", "p300", "blue")






  
#bcd_mel_2010 = make_hist("BCD_2010_mel.bedgraph", "bicoid_2010mel_eisen")
#bcd_yak_2010 = make_hist("BCD_2010_yak.bedgraph", "bicoid_2010yak_eisen")
#bcd_mel_2013 = make_hist("bcd_2013_mel.bedgraph", "bicoid_2013mel_eisen")
#bcd_yak_2013 = make_hist("bcd_2013_yak.bedgraph", "bicoid_2013yak_eisen")
#bcd_pse_2013 = make_hist("bcd_2013_pse.bedgraph", "bicoid_2013pse_eisen")
#bcd_vir_2013 = make_hist("bcd_2013_vir.bedgraph", "bicoid_2013vir_eisen")
#cad_mel_2010 = make_hist("CAD_2010_mel.bedgraph", "cad_2010mel_eisen")
#gt_mel_2010 = make_hist("GT_2010_mel.bedgraph", "giant_2010mel_eisen")
#kni_mel_2010 = make_hist("KNI_2010_mel.bedgraph", "knirps_2010mel_eisen")
#hb1_mel_2010 = make_hist("HB1_2010_mel.bedgraph", "hb1_2010mel_eisen")
#hb2_mel_2010 = make_hist("HB2_2010_mel.bedgraph", "hb2_2010mel_eisen")
#kr1_mel_2010 = make_hist("KR1_2010_mel.bedgraph", "kr1_2010mel_eisen")
#kr2_mel_2010 = make_hist("KR2_2010_mel.bedgraph", "kr2_2010mel_eisen")
#gt_mel_2013 = make_hist("gt_2013_mel.bedgraph", "gt_2013mel_eisen")
#hb_mel_2013 = make_hist("hb_2013_mel.bedgraph", "hb_2013mel_eisen")
#kr_mel_2013 = make_hist("kr_2013_mel.bedgraph", "kr_2013mel_eisen")

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
#furlong_snail = make_hist("snail_furlong.bedgraph", "Furlong_Snail")
#furlong_twist = make_hist("Twist_Furlong.bedgraph", "Furlong_Twist")

#eisen_mel_bcd = [bcd_mel_2010, bcd_mel_2013]
#eisen_mel_names = ["bcdmel2010", "bcdmel2013"]
#eisen_yak_bcd = [bcd_yak_2010, bcd_yak_2013]
#eisen_yak_names = ["bcdyak2010", "bcdyak2013"]
#eisen_bcd_2010 = [bcd_mel_2010, bcd_yak_2010]
#eisen_bcd_2010_names = ["bcdmel2010", "bcdyak2010"]
#eisen_bcd_2013 = [bcd_mel_2013, bcd_yak_2013, bcd_pse_2013, bcd_vir_2013]
#eisen_bcd_2013_names = ["bcdmel2013", "bcdyak13", "bcdpse13", "bcdvir13"]
#eisen_mel_2010 = [cad_mel_2010, kni_mel_2010]
#eisen_mel_2010_names = ["cadmel2010", "knimel2010"]
#eisen_mel_2013 = [bcd_mel_2013, gt_mel_2013,hb_mel_2013,kr_mel_2013]
#eisen_mel_2013_names = ["bcdmel2013","gtmel2013", "hbmel2013", "krmel2013"] 


#Furlong_data = [furlong_snail, furlong_twist]
#Furlong_names = ["FurlongSnail", "FurlongTwist"]

#make_boxplot(eisen_mel_2010, eisen_mel_2010_names, "eisen 2010 datasets")
#make_boxplot(eisen_mel_2013, eisen_mel_2013_names, "eisen 2013 datasets")
#make_boxplot(Furlong_data, Furlong_names, "Furlong chip-chip data")

#make_boxplot(eisen_scores, eisen_names, "eisen_boxplot1")
#make_boxplot(White_scores, White_names, "White_boxplot1")
#make_boxplot(Kurt_scores, Kurt_names, "Kurtulus boxplot1")
#make_boxplot(eisen_mel_bcd, eisen_mel_names, "eisen mel bicoid")
#make_boxplot(eisen_yak_bcd, eisen_yak_names, "eisen yak bicoid")
#make_boxplot(eisen_bcd_2010, eisen_bcd_2010_names, "eisen bicoid 2010")
#make_boxplot(eisen_bcd_2013, eisen_bcd_2013_names, "eisen bicoid 2013")


DVvals = (68, 40, 11, 36, 95, 5, 46)
APvals = (109, 10, 16, 33, 12, 33, 89)

DV = venn3(subsets = DVvals, set_labels = ("Ect", "End", "Mes"))

DV.get_patch_by_id('100').set_color("red")
DV.get_patch_by_id('010').set_color("yellow")
DV.get_patch_by_id('001').set_color("blue")
DV.get_patch_by_id('100').set_alpha(0.8)
DV.get_patch_by_id('010').set_alpha(0.8)
DV.get_patch_by_id('001').set_alpha(0.8)
plt.title("DV venn diagram")
savefig("DV_venn_diagram")
plt.close()


AP = venn3(subsets = APvals, set_labels = ("A", "C", "P"))

AP.get_patch_by_id('100').set_color("green")
AP.get_patch_by_id('010').set_color("purple")
AP.get_patch_by_id('001').set_color("brown")
AP.get_patch_by_id('100').set_alpha(0.8)
AP.get_patch_by_id('010').set_alpha(0.8)
AP.get_patch_by_id('001').set_alpha(0.8)
plt.title("AP venn diagram")
savefig("AP_venn_diagram")
plt.close()


