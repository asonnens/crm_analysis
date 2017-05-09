# This is a script for taking formatted output from read_stark.py and 
# determining the average characterisitcs of each putative enhancer, filtering
# for pre-specified criteria, measuring euclidean distance between enhancers,
# and undersampling the majority category
# it also replaces "0" 'no-binding' values with a value that's 1 standard deviation lower
# than the lowest value for that category

import matplotlib
from random import shuffle
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
Nope_file = open("Inactive_unedited.tsv")
Yep_file = open("Active_unedited.tsv")



#This function converts each file into a dictionary 
#dictionary keys are each enhancers name
def make_dict(file1, file2):
    expr_dict = {}
    line1 = file1.readline()
    file2.readline()
    keylist = line1.split("\t")
    for line in file1:
        linelist = line.split("\t")
        for i in linelist:
            i.strip()
        expr_dict[linelist[0]] = linelist
    for line in file2:
        linelist = line.split("\t")
        for i in linelist:
            i.strip()
        expr_dict[linelist[0]] = linelist
    return(expr_dict, keylist)


#this function takes the dict from make_dict and 
#creates a new dictionary that's filtered for enhancers
#that have binding in the specific key
def filter_dict(dict, keyinfo):
    print "\nselect a key from the following options, or type 'none':\n " 
    print ', '.join(str(k) for k in keyinfo)
    filter_dict1 = {}
    min_val = 0
    val_list = []
    keyname = raw_input()
    print "you have selected ",keyname
    print "what samplesize do you want to use?"
    samplesize = raw_input()
    keyinfo = map(str.strip, keyinfo)
    if keyname != "none":
        i = keyinfo.index(keyname)
        for keys, values in dict.items():
            val_list.append(float(values[i]))
    min_val = min(val_list)
    if keyname != "none":
        i = keyinfo.index(keyname)
        for keys, values in dict.items():
            if float(values[i]) > min_val:  
                filter_dict1[keys] = values
    else:
        filter_dict1 = dict
    return(filter_dict1, keyname, samplesize)


#right now this function prints the average of all enhancer 
#characteristics for each filtered enhancer, and 
#makes histograms of these scores
def average_values(keyinfo, newdict):
    mean_list_on = []
    mean_list_off = []
    for keys, values in newdict.items():
        if values[1] == "on":
            a = np.array(map(float, values[3:]))
            #print np.mean(a)
            mean_list_on.append(np.mean(a))
        else:
            a = np.array(map(float, values[3:]))
            #print np.mean(a)
            mean_list_off.append(np.mean(a))           
    b = np.array(map(float, mean_list_on))
    c = np.array(map(float, mean_list_off))
    plt.hist(c, bins = 40, histtype = 'stepfilled', color = 'b', label = "Off subset")
    plt.hist(b, bins = 40, histtype = 'stepfilled', color = 'r', alpha = 0.5, label = "On subset")
    plt.title("Zeitlinger_Twist filtered histogram")
    plt.xlabel("enhancer average characteristic score")
    plt.ylabel("frequency")
    plt.legend()
    savefig("ZeitlingerTwist_filter")
    plt.close()
        

#This function plots the histogram of value scores for on vs off subsets
#by a user-selected feature
def given_value(keyinfo, newdict):
    value_list_on = []
    value_list_off = []
    print "which value to you want to plot? type 'none' if none\n"
    keyname = raw_input()
    if keyname != "none":
        keyinfo = map(str.strip, keyinfo)
        for keys, values in newdict.items():
            i = keyinfo.index(keyname)
            a = np.array(float(values[i]))
            if values[1] == "on":  
                value_list_on.append(a)
            if values[1] == "off":
                value_list_off.append(a)
        b = np.array(map(float, value_list_on))
        c = np.array(map(float, value_list_off))
        plt.hist(c, bins = 40, histtype = 'stepfilled', color = 'b', label = "Off subset")
        plt.hist(b, bins = 40, histtype = 'stepfilled', color = 'r', alpha = 0.5, label = "On subset")
        plt.title(keyname + " histogram")
        plt.xlabel(keyname + " score")
        plt.ylabel("frequency")
        plt.legend()
        savefig(keyname + " scores")
        plt.close()
    

def euclid_dist(keyinfo, newdict, samplenum):
    samplesize = int(samplenum)
    list_on = [] 
    list_off = []
    sum_list_on = [0] * 24
    sum_list_on_2 = [0] * 24
    euc_on = []
    euc_off = []
    for keys, values in newdict.items():
        if values[1] == "on":
            list_on.append(values)
        else:
            list_off.append(values)
    print len(list_on), len(list_off)
    #conor addition    
    temp_list_on=[[y.replace('\n','') for y in x] for x in list_on]
    temp_list_on=np.array(temp_list_on)
    temp_list_on=temp_list_on[:,2:].astype('float')
    #zscore cutting
    temp_list_on_scaled=(temp_list_on-np.mean(temp_list_on,axis=0))/np.std(temp_list_on,axis=0)
    median_list_on=np.mean(temp_list_on_scaled,axis=0)
    '''
    for vallist in list_on:
        for i in range(3,len(vallist)):
            sum_list_on[i] += float(vallist[i])
    mean_list_on = [ x/len(list_on) for x in sum_list_on ]
    '''
    for vN,vallist in enumerate(list_on):
        simplesum = 0
        for i in range(3,len(vallist)):
            simplesum = simplesum + ((median_list_on[i-3] - (float(temp_list_on_scaled[vN,i-3])))**2)
        euc_on.append(simplesum**0.5)
        vallist.append(simplesum**0.5)  
    list_on.sort(key=lambda x: x[-1]) 
    new_list = list_on[0:samplesize]
    outlier_on = int(0.1 * len(list_on))
    all_on = list_on[0:(len(list_on))]
    euc_on_300 = [i[-1] for i in new_list]
    temp_list_on_scaled=(temp_list_on-np.mean(temp_list_on,axis=0))/np.std(temp_list_on,axis=0)
    median_list_on_short=np.mean(temp_list_on_scaled,axis=0)
    #now do it to list off
    temp_list_off=[[y.replace('\n','') for y in x] for x in list_off]
    temp_list_off=np.array(temp_list_off)
    temp_list_off=temp_list_off[:,2:].astype('float')
    #zscore cutting IMPORTANT SCALE TO SAME VALUES AS LIST ON (it might make sense to scale them all together whatever)
    temp_list_off=(temp_list_off-np.mean(temp_list_on,axis=0))/np.std(temp_list_on,axis=0)
    '''
    for vallist in new_list:
        for i in range(3,len(vallist)-1):
            sum_list_on_2[i] += float(vallist[i])
    mean_list_on_short = [ x/len(new_list) for x in sum_list_on_2 ]
    '''
    for vN,vallist in enumerate(list_off):
        simplesum = 0
        for i in range(3,len(vallist)):
            simplesum = simplesum + ((median_list_on_short[i-3] - (float(temp_list_off[vN,i-3])))**2)
        euc_off.append(simplesum**0.5) 
        vallist.append(simplesum**0.5)
    list_off.sort(key = lambda x: x[-1])
    r_list = list_off[50:-50]
    shuffle(r_list)
    if samplesize == 300:
        conservative_list = list_off[50:350]
        lax_list = list_off[-350:-50]
        random_list = r_list[0:300]
        all_list = list_off
    elif samplesize == 250:
        conservative_list = list_off[50:300]
        lax_list = list_off[-300:-50]
        random_list = r_list[0:250]
        all_list = list_off
    elif samplesize == 100:
        conservative_list = list_off[50:150]
        lax_list = list_off[-200:-100]
        random_list = r_list[0:100]
        all_list = list_off
    elif samplesize == 70:
        conservative_list = list_off[0:70]
        lax_list = list_off[-200:-130]
        random_list = r_list[0:70]
        all_list = list_off
    else:
        print "Error: sample size is: ", samplesize
    print len(conservative_list), len(lax_list), len(random_list), len(all_list)
    return(new_list, all_on, conservative_list, lax_list, random_list, all_list)


def make_euc_plot(conservative_list, lax_list, random_list, all_list, dataname):
    cons_nums = np.array([x[-1] for x in conservative_list])
    lax_nums = np.array([x[-1] for x in lax_list])
    rand_nums = np.array([x[-1] for x in random_list])
    all_nums = np.array([x[-1] for x in all_list])
    plt.hist(all_nums, bins = 8, histtype = 'step', color = 'k', label = "All")
    plt.hist(cons_nums, bins = 3, histtype = 'step', color = 'r', label = "Conservative")
    plt.hist(lax_nums, bins = 3, histtype = 'step', color = 'b', label = "Lax")
    plt.hist(rand_nums, bins = 3, histtype = 'step', color = 'g', label = "Random")
    plt.ylim(ymax = 1000)
    plt.title("Euclidean Subsets_" + dataname)
    plt.xlabel("Euclidean distance")
    plt.ylabel("frequency")
    plt.legend()
    savefig("Euclidean_comparison_" + dataname)
    plt.close()

def make_output(on_list, off_list, keywords, filename):
    outfile = open(filename, "w")
    keywords_string = "\t".join(keywords)
    outfile.write(keywords_string)
    for i in on_list:
        out_i = i[0:-1]
        out_i = [x if x != "0" else "NA" for x in out_i]
        out_string = "\t".join(out_i)
        out_string = out_string.strip()
        outfile.write(out_string)
        outfile.write("\r\n")
    for i in off_list:
        out_i = i[0:-1]
        out_i = [x if x != "0" else "NA" for x in out_i]
        out_string = "\t".join(out_i)
        outfile.write(out_string)
    outfile.close()

def replace_unbound(active_file, inactive_file):
    active_file.seek(0)
    first_line = active_file.readline()
    first_line_list = first_line.split("\t")
    list_of_mins = ["enhancer", "status"]
    for i in range(2,len(first_line_list)):
        temp_list = []
        active_file.seek(0)
        inactive_file.seek(0)
        active_file.readline()
        inactive_file.readline()
        for line in active_file:
            linelist = line.split("\t")
            temp_list.append(float(linelist[i]))
        for line in inactive_file:
            linelist = line.split("\t")
            temp_list.append(float(linelist[i]))
        print np.min(temp_list), np.std(temp_list), str(np.min(temp_list) - np.std(temp_list))
        list_of_mins.append([np.min(temp_list), str(np.min(temp_list) - np.std(temp_list))])
    return(list_of_mins)

def drop_zeros(file_list, list_of_mins):
    for i_list in file_list:
        i = i_list[0]
        i.seek(0)
        first_line = i.readline()
        outfile = open(i_list[1], "w")
        outfile.write(first_line)
        for line in i: 
            count = 0
            out_i = []       
            line_list = line.split("\t")
            for index_count in range(0, len(line_list)):
                if count < 2:
                    out_i.append(str(line_list[count]))
                else:
                    if float(line_list[count]) != 0:
                        out_i.append(str(line_list[count]))
                    else:
                        out_i.append(str(list_of_mins[count][1]))
                count = count + 1
            out_string = "\t".join(out_i)
            outfile.write(out_string)
        outfile.close()
               

min_list = replace_unbound(Yep_file, Nope_file)
target_file = open("target_unedited.tsv")
#file_list = [[Yep_file, "Active.tsv"], [Nope_file, "Inactive.tsv"], [target_file, "target_file.tsv"]]
file_list = [[target_file, "target_file.tsv"]]
drop_zeros(file_list, min_list)
Yep_file.close()
Nope_file.close()
Yep_file = open("Active.tsv")
Nope_file = open("Inactive.tsv")
reporter_dict, keylist  = make_dict(Yep_file, Nope_file)
new_dict, keyname, samplenum = filter_dict(reporter_dict, keylist)        
average_values(keylist, new_dict)
given_value(keylist, new_dict)
list_on, all_on, cons_list, lax_list, rand_list, all_list = euclid_dist(keylist, new_dict, samplenum)
make_euc_plot(cons_list, lax_list, rand_list, all_list, keyname)
outputname1 = "cons_train_" + keyname + ".tsv"
outputname2 = "lax_train_" + keyname + ".tsv"
outputname3 = "rand_train_" + keyname + ".tsv"
outputname4 = "all_train_" + keyname + ".tsv"
make_output(list_on, cons_list, keylist, outputname1)
make_output(list_on, lax_list, keylist, outputname2)
make_output(list_on, rand_list, keylist, outputname3)
make_output(all_on, all_list, keylist, outputname4)

Nope_file.close()
Yep_file.close()
