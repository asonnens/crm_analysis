import enhancer_classes
input = "../input_data/input_files/Stark_data"

#this has been really modified to get everything into one text file...it throws out expression info


class Reporter(object):
    """This class includes attributes of reporters with no assigned genes"""
    ###########################################################################
    #Reporter attributes:  coordinates (list), expression(list of lists), dataset(string) 
    #experiement in which identified
    ###########################################################################
    def __init__(self, coords, expression, dataset):
        self.coords = coords
        self.dataset = dataset
        self.expression = expression

        

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
            else:
                expr = expression_list[0].split(";")
                expression = [[expr[0],expr[1]]]
            #else:
            #    expr = i.split(";")
            #    expression = expression + [[expr[0],expr[1]]]
        my_reporter = Reporter(coords, expression, "Stark")
        coords_string = "\t".join(coords)
        each_reporter = enhancer_classes.Enhancer("Kvon2014", coords_string, "NA", expression, "reporter")
        reporter_list.append(each_reporter)            
    myfile.close()
    return(reporter_list)

my_reporters = read_Stark(input, "Kvon2014")
one_list = []
NA_list = []
multi_list = []
count = 0
for i in my_reporters:
    count = count + 1
    if count%10 == 0:
        print count,"/",len(my_reporters)
    for j in i.expr:
        if len(i.expr) == 1 and int(j[1]) >=1:
            for bedgraph in enhancer_classes.bedfiles:
                i.read_file(bedgraph)
            one_list.append(i)
        elif len(i.expr) == 1 and j[1] == 0:
            for bedgraph in enhancer_classes.bedfiles:
                i.read_file(bedgraph)
            NA_list.append(i)
  #      elif len(i.expr) >= 2:
  #          big_list = []
  #          for j in i.expr:
  #              if int(j[1]) > 3:
  #                  big_list.append(j)
  #          if len(big_list) == 1:
  #              i.expr = big_list
  #              for bedgraph in enhancer_classes.bedfiles:
  #                  i.read_file(bedgraph)
  #              multi_list.append(i)
        else:
            pass


#stark_all_comparison = enhancer_classes.Enhancer_compare(all_list, "Stark_reporters_all")
#stark_one_reporter_comparison = enhancer_classes.Enhancer_compare(one_list, "Stark_reporters_one")
stark_NA_reporter_comparison = enhancer_classes.Enhancer_compare(NA_list, "Stark_reporters_NA")
#stark_multi_reporter_comparison = enhancer_classes.Enhancer_compare(multi_list, "Stark_reporters_multi")
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_one_reporter_comparison)
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_all_comparison)
enhancer_classes.Enhancer_compare.print_enhancer_data(stark_NA_reporter_comparison)
#enhancer_classes.Enhancer_compare.print_enhancer_data(stark_multi_reporter_comparison)
