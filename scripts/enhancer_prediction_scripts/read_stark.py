#functions for reading and parsing the data from Kvon 2014
#because it uses enhancer_classes requires scipy and numpy

import enhancer_classes

input = "/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/Stark_data"

expression_input = "/mnt/home/sonnens2/crm_analysis/datasets/enhancer_prediction_input_data/enhancer_prediction_input_files/Stark_expression"

        
#this function reads the file, and converts each line into a list
#with expression information as a sublist

def read_Stark(filename, dataset_name):
    reporter_list = []
    myfile = open(filename)
    for line in myfile:
        line = line.strip()
        myreporter = line.split("\t")
        chr = myreporter[0]
        start_coord = myreporter[1]
        end_coord = myreporter[2]
        coords = [chr, start_coord, end_coord]
        expression_info = myreporter[4]
        expression_list = expression_info.split("|")
        expression = []
        for i in expression_list:
            if i == "NA":
                expression = expression + [["NA",0]]
            elif i == "on":
                expression = expression + [["on","NA"]]
            else:
                expr = i.split(";")
                expression = expression + [[expr[0],expr[1]]]
        coords_string = "\t".join(coords)
        each_reporter = enhancer_classes.Enhancer("Kvon2014", coords_string, "NA", expression, "reporter")
        reporter_list.append(each_reporter)            
    myfile.close()
    return(reporter_list)

my_reporters = read_Stark(input, "Kvon2014")
#my_reporters = read_Stark(expression_input, "Kvon2014")


#this function takes the list output from the function "read_Stark"
#any reporter with expression "NA" is added to NA list
#any reporter with only one expression pattern of score > 3 is added to one_list
#any reporter with at least one expression pattern of score > 3 is added to all_expression_list
#any reporter with more than one expression pattern, with at least one of score > 3, is added to multi_list
#the enhancer classes are used to add features that overlap their coordinates to each line
#output from this is passed to the enhancer classes to be printed

def print_Stark(my_reporters):
    all_expression_list = []
    NA_list = []
    one_list = []
    multi_list = []
    expression_list = []
    count = 0
    for i in my_reporters:
        count = count + 1
        if count%10 == 0:
            print count,"/",len(my_reporters)
        for j in i.expr:
            if len(i.expr) == 1:
                if j[1] == "NA":
                    for bedgraph in enhancer_classes.bedfiles:
                        i.read_file(bedgraph)
                    expression_list.append(i)
                elif int(j[1]) == 0:
                    for bedgraph in enhancer_classes.bedfiles:
                        i.read_file(bedgraph)
                    NA_list.append(i)
                elif int(j[1]) >= 3:
                    for bedgraph in enhancer_classes.bedfiles:
                        i.read_file(bedgraph)
                    one_list.append(i)
                    all_expression_list.append(i)
                else:
                    pass
            else:
                if int(j[1]) < 3:
                    pass
                else:
                    if i not in all_expression_list:
                        for bedgraph in enhancer_classes.bedfiles:
                            i.read_file(bedgraph)
                        all_expression_list.append(i)
                        multi_list.append(i)
                    else:
                        pass
    return(one_list, multi_list, all_expression_list, NA_list, expression_list)

all_lists = print_Stark(my_reporters)
                    
#(self, gene, coords, distance, expr, status)
one_list = all_lists[0]
multi_list = all_lists[1]
all_expression_list = all_lists[2]
NA_list = all_lists[3]
expression_list = all_lists[4]


#stark_all_comparison = enhancer_classes.Enhancer_compare(expression_list, "Stark_expression")
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_all_comparison)
#stark_one_reporter_comparison = enhancer_classes.Enhancer_compare(one_list, "Stark_reporters_e")
#stark_NA_reporter_comparison = enhancer_classes.Enhancer_compare(NA_list, "Stark_reporters_NA_Aug28")
stark_all_reporter_comparison = enhancer_classes.Enhancer_compare(all_expression_list, "Stark_reporters_all_test")
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_one_reporter_comparison)
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_NA_reporter_comparison)
enhancer_classes.Enhancer_compare.print_enhancer_data(stark_all_reporter_comparison)
