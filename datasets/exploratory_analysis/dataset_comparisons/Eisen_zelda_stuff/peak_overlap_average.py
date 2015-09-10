#misc functions for analyizing other peoples chip-chip and chip-seq output, 
#and determining the percent of peak overlap as well as
#summary statistics


import numpy


#function for finding and counting the number of overlapping peaks
#this is fed into the file_stats function, which actually outputs the information

def peak_overlap(smaller_file, bigger_file):
    overlap_count = 0
    smaller_file.readline()
    smaller_file.readline()
    smaller_file.readline()
    for line in smaller_file:
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
    return(overlap_count)


#function for calculating summary statistics
#including difference in number of peaks between two antibodies
#percent overlap between two antibodies

def file_stats(filename1, filename2, fileinfo = ""):
    difference_filesize = 0
    percent_overlap = 0
    file1 = open(filename1)
    file2 = open(filename2)
    file1_length = sum(1 for line in file1)
    file2_length = sum(1 for line in file2)
    print "file 1 is ", file1_length, " file 2 is ", file2_length
    file1.seek(0)
    file2.seek(0)
    if file2_length > file1_length:
        difference_filesize = 100 * (file2_length - file1_length)/file1_length
        percent_overlap = 100 * peak_overlap(file1, file2)/file1_length
    else:
        difference_filesize = 100 *  (file1_length - file2_length)/file2_length
        percent_overlap = 100 * peak_overlap(file2, file1)/file2_length
    print fileinfo, " difference in filesize is ", difference_filesize, "peak overlap percentage is ", percent_overlap  
    file1.close()
    file2.close()
    return(percent_overlap, difference_filesize)

#within an individual file, calculates the average score for a chip-peak, and the standard deviation around that mean
def each_file_stats(filename, fileinfo = ""):
    file1 = open(filename)
    file1.readline()
    bound_region_list = []
    peak_score_list = []
    for line in file1:
        line = line.strip()
        info = line.split("\t")
        if len(info) == 4:
            bound_region_list.append(float(info[3]))
        #peak_score_list.append(float(info[6]))
    print fileinfo, numpy.average(bound_region_list), numpy.std(bound_region_list)
    file1.close()





file_stats("zelda_eisen_cycle13.bedgraph", "zelda_eisen_cycle14.bedgraph", "Eisen 13 14 comparison")
#file_stats("sna_1_081506-sym-1_table.txt", "sna_2_081506-sym-1_table.txt", "snail files")
#file_stats("hry_1_040108-sym-1_table.txt", "hry_2_040108-sym-1_table.txt", "hairy files")
#file_stats("bcd_1_012505-sym-1_table.txt", "bcd_2_092005-sym-1_table.txt", "bicoid files")
#file_stats("hb_1_012505-sym-1_table.txt", "hb_2_092305-sym-1_table.txt", "hunchback files")


each_file_stats("zelda_eisen_cycle13.bedgraph", "Eisen cycle 13 Zelda")
each_file_stats("zelda_eisen_cycle14.bedgraph", "Eisen cycle 14 Zelda")
#each_file_stats("sna_1_081506-sym-1_table.txt", "snail file 1")
#each_file_stats("sna_2_081506-sym-1_table.txt", "snail file 2")
#each_file_stats("hry_1_040108-sym-1_table.txt", "hairy file 1")
#each_file_stats("hry_2_040108-sym-1_table.txt", "hairy file 2")
#each_file_stats("bcd_1_012505-sym-1_table.txt", "bicoid file 1")
#each_file_stats("bcd_2_092005-sym-1_table.txt", "bicoid file 2")
#each_file_stats("hb_1_012505-sym-1_table.txt", "hunchback file 1")
#each_file_stats("hb_2_092305-sym-1_table.txt", "hunchback file 2")
#each_file_stats("gt_2_020107-sym-1_table.txt", "giant file")
#each_file_stats("cad_1_020107-sym-1_table.txt", "caudal file")
#each_file_stats("dl_3_120106-sym-1_table.txt", "dorsal file")


#file_stats("twi_1_081506-sym-25_table.txt", "twi_2_081506-sym-25_table.txt", "twist files 25")
#file_stats("sna_1_081506-sym-25_table.txt", "sna_2_081506-sym-25_table.txt", "snail files 25")
#file_stats("hry_1_040108-sym-25_table.txt", "hry_2_040108-sym-25_table.txt", "hairy files 25")
#file_stats("bcd_1_012505-sym-25_table.txt", "bcd_2_092005-sym-25_table.txt", "bicoid files 25")
#file_stats("hb_1_012505-sym-25_table.txt", "hb_2_092305-sym-25_table.txt", "hunchback files 25")

#each_file_stats("Zeitlinger_tw_chipseq.bedgraph", "Zeitlinger file 1")
#each_file_stats("Zeitlinger_tw_chipseq2.bedgraph", "Zeitlinger file 2")
#each_file_stats("sna_1_081506-sym-25_table.txt", "snail file 25 1")
#each_file_stats("sna_2_081506-sym-25_table.txt", "snail file 25 2")
#each_file_stats("hry_1_040108-sym-25_table.txt", "hairy file 25 1")
#each_file_stats("hry_2_040108-sym-25_table.txt", "hairy file 25 2")
#each_file_stats("bcd_1_012505-sym-25_table.txt", "bicoid file 25 1")
#each_file_stats("bcd_2_092005-sym-25_table.txt", "bicoid file 25 2")
#each_file_stats("hb_1_012505-sym-25_table.txt", "hunchback file 25 1")
#each_file_stats("hb_2_092305-sym-25_table.txt", "hunchback file 25 2")
#each_file_stats("gt_2_020107-sym-25_table.txt", "giant file 25 ")
#each_file_stats("cad_1_020107-sym-25_table.txt", "caudal file 25")
#each_file_stats("dl_3_120106-sym-25_table.txt", "dorsal file 25")



