###############################################################################
#Anne Sonnenschein
#2014
#classes for finding and annotating transcription factors and enhancers
###############################################################################


###############################################################################
import glob
import scipy, numpy
import math
from collections import defaultdict
###############################################################################


bedfiles = glob.glob("/mnt/home/sonnens2/crm_analysis/input_data/input_files/*.bedgraph")
phast_files = glob.glob("/mnt/home/sonnens2/crm_analysis/other_files/phastcons/*.pp")
histone_files = glob.glob("/mnt/home/sonnens2/crm_analysis/other_files/Furlong_bits/*.wig")


class TF(object):
    """This class includes attributes of putative transcription factors"""
    ###########################################################################
    #TF attributes:  TFprotein, coordinates, method of identification, 
    #experiement in which identified
    ###########################################################################
    def __init__(self, protein, coords, score, method, dataset):
        self.protein = protein
        self.coords = coords
        self.score = score
        self.method = method
        self.dataset = dataset
        self.name = str(dataset) + "_" + str(protein) + "_" + str(coords)
        chromosome_region = self.coords.split("\t")
        self.chr = chromosome_region[0]
        self.coord1 = int(chromosome_region[1])
        self.coord2 = int(chromosome_region[2])
        self.feature_type = self.dataset + "," + self.protein + "," + self.method


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



class Enhancer(object):
    """This class includes attributes and functions of putative enhancers"""
    ###########################################################################
    #Attributes: associated gene, coordinates, distance from TSS of associated 
    #gene, expression pattern of the gene it drives, and associated features/TF
    ###########################################################################
    def __init__(self, gene, coords, distance, expr, status):
        """initiates enhancer class-- gene, coordinates, expr, etc."""
        self.gene = gene
        self.coords = coords
        self.distance = distance
        self.expr = expr
        self.status = status
        chromosome_region = self.coords.split("\t")
        self.chr = chromosome_region[0]
        self.coord1 = int(chromosome_region[1])
        self.coord2 = int(chromosome_region[2])
        self.TF_dict = {}
        self.phast_cons = 0
        self.H3 = 0
        self.chromatin = {}
        self.feature_list = []
        self.score = 5
        self.name = gene + "_" + self.chr + "_" + str(self.coord1) + "_" + str(self.coord2)

    ###########################################################################
    #Functions for reading files into enhancer, and printing enhancer info
    ########################################################################### 

    def add_TF(self, TFactor):
        #########################################################
        #function adds TF objects to enhancer's TF_dict
        #########################################################
        """This function reads bedfiles adds overlapping TFs"""  
        self.TF_dict[TFactor.name] = TFactor

    def check_file_type(self, infile):
        #########################################################
        #function checks filetype based on header lines
        #########################################################
        filetype = ""
        first_line = infile.readline()
        if "browser hide all" in first_line:
            filetype = "bedgraph"
        elif "fixedStep" in first_line:
            filetype = "wig"
        elif "type=wiggle" in first_line:
            filetype = "wig"
            first_line = infile.readline()
        elif "##gff" in first_line:
            filetype = "gff"
        else:
            filetype = "unknown file type"
        return filetype, first_line
     
     ###########################################################################
     #Functions for bedfiles
     ########################################################################### 

    def filter_bedfile(self, chr, TF_to_check, dataset, method, protein):
        #########################################################
        #adds TFs to the enhancer based on overlap
        #first function for reading bedfiles
        #checks determines if bedfile line falls within enhancer
        #checks chromosome, if chromosome matches, checks for overlap
        #########################################################
        bed_line = TF_to_check.split("\t") 
        if "chr" not in chr:
            chrline = "chr" + str(chr) #making sure nomenclature is consistent
        else:
            chrline = chr
        if bed_line[0] == chrline:
            TF_coord1 = int(bed_line[1])
            TF_coord2 = int(bed_line[2])
            if (self.coord1 > TF_coord2) or (TF_coord1 > self.coord2): 
               #conditions of NOT overlapping features
                pass
            else:
                 TF_coords = chrline + "\t" + str(TF_coord1) + "\t" + str(TF_coord2)
                 score = bed_line[3].strip()  
                 TFactor = TF(protein, TF_coords, score, method, dataset)
                 self.add_TF(TFactor)    

    def read_bedgraph(self, infile):
        #########################################################
        #second function for reading bedgraphs
        #processes header info, and checks filetype, before sending
        #to filter_bedfile function
        #########################################################
        for line in infile:
            if "browser" in line:
                pass
            elif "track type" in line:
                line = line.strip()
                header_info = line.split("=")
                filename_info = header_info[2].split(" ")
                dataset = filename_info[0].strip('\"')
                method = filename_info[1].strip('\"')
                protein = filename_info[2].strip('\"')    
                feature_type = dataset + "," + protein + "," + method
                if feature_type not in self.feature_list:
                    self.feature_list.append(feature_type)    
            else:
                self.filter_bedfile(self.chr, line, dataset, method, protein)

     ###########################################################################
     #Functions for wigfiles
     ########################################################################### 

    def get_wigfile_info(self, wig_info_line):
        #########################################################
        #first function for reading wigfiles
        #for wigfile, if information/header line
        #returns the chromosome, start, and step info
        #########################################################
        wig_info = wig_info_line.split(" ")
        chrom_info = wig_info[1]
        wig_chr = chrom_info.split("=")[1]
        start_info = wig_info[2]
        wig_start = int(start_info.split("=")[1])
        step_info = wig_info[3]
        wig_step = int(step_info.split("=")[1])
        return wig_chr, wig_start, wig_step

    def read_wigfile(self, infile, first_line, filename):
        #########################################################
        #second function for reading wigfiles
        #if header line, gets header info with 
        #get_wigfile_info function
        #########################################################
        file_chr, start, step, = self.get_wigfile_info(first_line)
        filename_extension = filename.split(".")[-1]
        filename_ID = filename.split(".")[0] 
        current_coord = start
        wig_info = []
        wig_sum = 0
        for line in infile:
            if "fixedStep" in line:
                file_chr, start, step, = self.get_wigfile_info(line)
                current_coord = start
            else:
                current_coord = current_coord + step
                if (self.chr == file_chr) and \
                   (current_coord >= self.coord1) and (current_coord <= self.coord2): 
                    current_val = line.strip()
                    wig_sum = wig_sum + float(current_val)
                    line_info =  str(self.chr) + "\t" + str(current_coord-step) + "\t" + \
                                 str(current_coord) + "\t" + str(current_val) + "\n"   
                    #wig_info.append(line_info)        
                else:
                    pass
        wig_average = wig_sum/(self.coord2 - self.coord1)
        if filename_extension == "pp":
            self.phast_cons = wig_average
        elif filename_extension == "wig":
            if "H3-subtracted" in filename_ID:
                filename_ID = filename_ID.strip("_H3-subtracted")
                self.chromatin[filename_ID] = wig_average 
            else:
                if "H3" in filename_ID:   
                    self.H3 = wig_average      
        else:
            print "unknown file type? ", filename

    def read_file(self, filename):
        #########################################################
        #reads files into enhancer
        #checks filetype using check_file_type function
        #depending on filetype, feeds file into appropriate function
        #########################################################
        infile = open(filename)
        filetype, first_line = self.check_file_type(infile)
        filename_ID = filename.split("/")[-1]
        #print "reading ", filename_ID + " for " + self.status + " enhancer of " + self.gene
        dataset = ""
        method = ""
        protein = ""
        if filetype == "bedgraph":
            self.read_bedgraph(infile)
        if filetype == "wig" or "wig2":
            if self.chr in infile.read():
                infile.seek(0)
                self.read_wigfile(infile, first_line, filename_ID)
        infile.close()

    def print_TFs(self):
        #########################################################
        #for each TF in TF_dict, prints dataset, TFprotein, coordinates, and 'score'
        #########################################################
        print "Enhancer of ", self.gene, " , ", self.distance, " from TSS , coordinates: ", self.coords
        for keys in self.TF_dict:
            print self.TF_dict[keys].dataset, self.TF_dict[keys].protein, self.TF_dict[keys].chr, \
                  self.TF_dict[keys].coord1, self.TF_dict[keys].coord2, self.TF_dict[keys].score
     #   for filename, score in self.other_features_dict:
     #       print filename, score

    def print_enhancer_bedfile(self, outfile):
        #########################################################
        #print enhancer coordinates, and score in bedfile format
        #########################################################
        my_line = self.chr + "\t" + str(self.coord1) + "\t" + \
                  str(self.coord2) + "\t" +  str(self.score) + "\r\n"
        outfile.write(my_line)



class Enhancer_compare(object):
    """this class takes in a set of enhancers, and provides \
       average values and variance for all numerical features"""
    ###########################################################################
    #This creates the object of an enhancer comparison-- basically a list of 
    #enhancers, and some summary statistics of their characteristcs
    #attributes of Enhancer_compare: the list of enhancers, any 
    #unincorporated feature files associated with these enhancers
    ###########################################################################

    def __init__(self, enhancer_list, list_name, input_info = bedfiles):
        #########################################################
        #attributes of Enhancer_compare object include:
        #enhancer list, and unincorporated feature files
        #########################################################
        print "comparing enhancers ", list_name
        self.list_name = list_name
        self.enhancer_list = enhancer_list
        self.input_info = input_info

    ###########################################################################
    #Enhancer_compare functions:
    #-calculate average length, distance from TSS within list
    #-calculate average scores, counts for a bedfile feature
    #-print summary info for enhancer to a csv file
    ###########################################################################

    def get_basic_stats(self, my_stat, my_list):
        #########################################################
        #calculates min, max, mean, and sd for a feature
        #########################################################
        my_array = numpy.array(my_list)
        average_val = numpy.mean(my_array)
        min_val = numpy.amin(my_array)
        max_val = numpy.amax(my_array)
        sd_val = numpy.std(my_array)
        print "average ", my_stat, ": ", average_val
        print "min ", my_stat, ": ", min_val, ", max ", my_stat, ": ", max_val
        print "sd ", my_stat, ": ", sd_val

    def average_enhancer_length(self):
        #########################################################
        #determines enhancer length, and 
        #gets average length for a set of enhancers
        #########################################################
        my_list = []
        for i in self.enhancer_list:
            enhancer_length = (i.coord2 - i.coord1)
            my_list.append(enhancer_length)
        get_basic_stats("length", my_list)

    def average_enhancer_TSSdist(self):
        #########################################################
        #determines distance from TSS, and average for a set of
        #enhancers
        #(distance is a required attribute for each enhancer)
        #########################################################
        my_list = []
        for i in self.enhancer_list:
            my_list.append(i.distance)
        get_basic_stats("distance from TSS", my_list)  

    def list_of_features(self):
        #########################################################
        #lists all the features names in the feature files
        #########################################################
        feature_list = []
        for i in self.enhancer_list:
            for j in i.feature_list:
                if j not in feature_list:
                    feature_list.append(j)   
        feature_list.sort()
        print len(feature_list)
        return feature_list

    def average_enhancer_feature(self):
        #########################################################
        #gets mean value of each feature for a list of enhancers
        #########################################################
        feature_list = self.list_of_features()
        average_feature_dict = {}
        mean_dict = {}
        for i in feature_list:
            average_feature_dict[i] = [[0], [0]]
        enhancer_count = 0
        for i in self.enhancer_list:
            enhancer_count = enhancer_count + 1
            enhancer_feature_info = {}
            for tfs in i.TF_dict:
                feature_type = i.TF_dict[tfs].feature_type
                if enhancer_feature_info.has_key(feature_type):
                    new_val = enhancer_feature_info[feature_type] + (i.TF_dict[keys].score, 1)
                    enhancer_feature_info[feature_type] = new_val
                else:
                    enhancer_feature_info[feature_type] = (i.TF_dict[keys].score, 1)
            for keys in enhancer_feature_info:
                new_val_1 = float(enhancer_feature_info[keys][0])/float(enhancer_feature_info[keys][1])
                new_val_2 = float(enhancer_feature_info[keys][1])
                enhancer_feature_info[keys] = [[new_val_1], [new_val_2]]
                if average_feature_dict[keys] == [[0], [0]]:
                    average_feature_dict[keys] = enhancer_feature_info[keys]
                else:
                    average_feature_dict[keys][0].append(new_val_1)
                    average_feature_dict[keys][1].append(new_val_2)
        for keys in average_feature_dict:
            my_list_1 = average_feature_dict[keys][0]
            my_array_1 = numpy.array(my_list_1)
            my_list_2 = average_feature_dict[keys][1]
            my_array_2 = numpy.array(my_list_2)
            mean_dict[keys] = [numpy.mean(my_array_1),numpy.std(my_array_1)] , \
                              [numpy.mean(my_array_2),numpy.std(my_array_2)]
        print mean_dict

    def print_enhancer_data(self):
        #########################################################
        #for each enhancer in comparison list, prints sum of
        #feature scores, and counts
        #########################################################
        outfilename = "/mnt/home/sonnens2/crm_analysis/output_data/" + \
                       self.list_name + "_" + "summary_info.csv"
        outfile = open(outfilename, "w")
        outfile.write("name,gene,expr,status,distance_from_TSS,phastcons,H3,")
        key_list= []
        feature_list = []
        feature_list_output = []
        header_info = "" 
        for enhancer_items in self.enhancer_list:
            for i in enhancer_items.feature_list:
                if i not in feature_list:
                    feature_list.append(i)
                    ilist = i.split(",")
                    ioutput = "_".join(ilist)
                    feature_list_output.append(ioutput)
            for keys in enhancer_items.chromatin:
                print keys
                if keys not in key_list:
                    key_list.append(keys)
        for tf_keys in sorted(feature_list_output):
            header_info = header_info + tf_keys + "_score,"
        for chromatin_keys in sorted(key_list):
            header_info = header_info + chromatin_keys + "_score,"       
        header_info = header_info[0:-1] + "\r\n"
        outfile.write(header_info)
        for i in self.enhancer_list:
            out_string = (str(i.name) + "," + str(i.gene) + "," + str(i.expr) + "," +   \
                          str(i.status) + "," +  str(i.distance) + "," + 
                          str(i.phast_cons) + "," + str(i.H3) + ",")           
            for j in sorted(feature_list):
                feature_count = 0
                feature_scores = [0]
                for tfs in i.TF_dict:
                    feature_type = i.TF_dict[tfs].feature_type
                    if j == feature_type:
                        feature_count = feature_count + 1
                        feature_scores.append(float(i.TF_dict[tfs].score))
                sum_score = numpy.sum(feature_scores)
                out_string = out_string + (str(sum_score) + "," )
            for keys in key_list:
                if keys in i.chromatin.keys():
                    out_string = out_string + str(i.chromatin[keys]) + ","
                else:
                    out_string = out_string + str(0) + ","
            out_string = out_string[0:-1]  + "\r\n" 
            outfile.write(out_string)
        outfile.close()

class Enhancer_find:
    ###########################################################################
    #checks feature files for putative enhancers associated with a specific 
    #gene, and validates features. includes list of TFs involved in putative 
    #enhancer, list of feature files, coordinates, 
    #and the gene for which enhancers are being found
    ###########################################################################

    def __init__(self, protein_list, coords, gene, TSS, expr, status, input_info = bedfiles):
        """instantiates Enhancer_find set of actions"""
        #########################################################
        #features include: list of proteins that may contribute to this enhancer, 
        #gene enhancer would be associated with, and target coords
        #########################################################
        print "finding enhancers of ", gene
        self.protein_list = protein_list
        self.input_info = input_info
        self.coords = coords
        coord_info = coords.split("\t")
        self.chr = coord_info[0]
        self.coord1 = int(coord_info[1])
        self.coord2 = int(coord_info[2])
        self.gene = gene
        self.TSS = TSS
        self.phast_cons = 0
        self.expr = expr
        self.status = status

    def cluster_TF(self, seed_TF, test_line):
        """uses seed protein and other transcription factors to identify clusters"""
        #########################################################
        #this function adds bedformat transcription factor features
        #if there is any overlap with the seed protein (in this case
        #dorsal), it's added to the list of prospective enhancers
        #The coordinates overlapping features for a putative
        #are averaged, weighting in favor of chip-seq over chip-chip
        #data 3:1 to get the midpoint, and 750 up and down from the 
        #midpoint are defined as the enhancer.         
        #the putative enhancer coords are returned
        #########################################################
        my_feature_list = []
        return_variable = 0
        for each_file in self.input_info:
            infile = open(each_file)
            dataset = ""
            method = ""
            protein = ""
            for line in infile:
                if "browser" in line:
                    pass
                elif "track type" in line:
                    line = line.strip()
                    header_info = line.split("=")
                    filename_info = header_info[2].split(" ")
                    dataset = filename_info[0].strip('\"')
                    method = filename_info[1].strip('\"')
                    protein = filename_info[2].strip('\"')
                else:
                    bed_line = line.split("\t")
                    chrline = seed_TF.chr
                    if bed_line[0] == chrline:
                        TF_coord1 = int(bed_line[1])
                        TF_coord2 = int(bed_line[2])
                        if (seed_TF.coord1 > TF_coord2) or (TF_coord1 > seed_TF.coord2): 
                            #conditions of NOT overlapping features
                            pass
                        elif (seed_TF.coord1 == TF_coord1) and (seed_TF.coord2 == TF_coord2) and \
                             (protein == seed_TF.protein) and (dataset == seed_TF.dataset):
                            pass
                        else:
                            TF_coords = chrline + "\t" + str(TF_coord1) + "\t" + str(TF_coord2)
                            score = bed_line[3].strip()
                            TFactor = TF(protein, TF_coords, score, method, dataset)
                            my_feature_list.append(TFactor)
        if len(my_feature_list) > 1:
            method_weights = {"chip-chip":1, "chip-seq":3}
            my_coord_list = []
            my_weight_list = []
            distance_from_TSS = 0
            for i in my_feature_list:
                my_coord_list.append((i.coord2 - i.coord1)/2 + i.coord1)
                my_weight_list.append(method_weights[i.method])
            midpoint = int(numpy.average((my_coord_list), weights = (my_weight_list)))   
            putative_enhancer_coord1 = midpoint - 1000
            putative_enhancer_coord2 = midpoint + 1000
            putative_enhancer_coords =  self.chr + "\t" + \
                                       str(putative_enhancer_coord1) + "\t" + \
                                       str(putative_enhancer_coord2)
            if self.TSS < midpoint:
                distance_from_TSS = int(math.fabs(putative_enhancer_coord1 - self.TSS))
            else:
                distance_from_TSS = int(math.fabs(putative_enhancer_coord2 - self.TSS))        
            putative_enhancer = Enhancer(self.gene, putative_enhancer_coords, \
                                         distance_from_TSS, self.expr, self.status)
            for i in my_feature_list:
                putative_enhancer.add_TF(i)
            return_variable = putative_enhancer   
        return return_variable 

    def get_TF(self, line, protein, method, dataset):
        #########################################################
        #the coordinates of putative enhancers identified in the 
        #clustering function are compared with the total list of
        #transcription factors.  TFs from the bedfile are added
        #to the putative enhancers 
        #########################################################
        putative_enhancer_list = []
        my_line = line.split("\t")
        chr = my_line[0]
        if chr == self.chr:
            TF_coord1 = int(my_line[1])
            TF_coord2 = int(my_line[2])
            if (TF_coord1 > self.coord2) or (self.coord1 > TF_coord2):
                pass
            else:
                coords = "\t".join(my_line[0:3])
                score = float(my_line[3])
                test_TF = TF(protein, coords, score, method, dataset)
                putative_enhancer = self.cluster_TF(test_TF, line)
                if putative_enhancer != 0:
                    putative_enhancer_list.append(putative_enhancer)
        return(putative_enhancer_list)

    def print_bedfile(self, mylist):
        #########################################################
        #print a list of enhancers in bedgraph format
        #########################################################
        file_info = self.gene + "_" + self.status + "_enhancers"
        outfile_name = "/mnt/home/sonnens2/crm_analysis/output_data/" + \
                       self.gene + "_putative_enhancers.bedgraph"    
        outfile = open(outfile_name, "w")
        outfile.write("browser hide all\r\n")
        outfile.write("browser pack refGene encodeRegions\r\n")
        outfile.write("browser full altGraph\r\n")
        outfile.write("track type=bedGraph name=")
        outfile.write(file_info)
        outfile.write(" description=\"BedGraph format\" visibility=full color=300,400,0 \
                                      altColor=0,400,300 priority=20\r\n")
        mylist.sort(key=lambda x: x.coord1, reverse=False)
        for i in mylist:
            i.print_enhancer_bedfile(outfile)
        outfile.close()

    def merge_enhancers(self, enhancer_list):
        #########################################################
        #
        #########################################################
        midpoint_list = []
        distance_from_TSS = 0
        for i in enhancer_list:
            enhancer_coords = [i.coord2, i.coord1]
            midpoint_list.append(numpy.average(enhancer_coords))
        midpoint =  numpy.average(midpoint_list).astype(int)
        merged_midpoint = numpy.asscalar(midpoint)
        merged_enhancer_coord1 = merged_midpoint - 1000
        merged_enhancer_coord2 = merged_midpoint + 1000
        if self.TSS < merged_midpoint:
            distance_from_TSS = int(math.fabs(merged_enhancer_coord1 - self.TSS))
        else:
            distance_from_TSS = int(math.fabs(merged_enhancer_coord2 - self.TSS)) 
        merged_coords = self.chr + "\t" + str(merged_enhancer_coord1) + \
                        "\t" + str(merged_enhancer_coord2)
        merged_enhancer = Enhancer(self.gene, merged_coords, distance_from_TSS, self.expr, self.status) 
        return merged_enhancer

    def merge_enhancer_list(self, enhancer_list):
        #########################################################
        #this is the worst function ever
        #########################################################
        enhancer_list.sort(key=lambda x: x.coord1, reverse=False)
        print "merging overlapping enhancers"
        small_list = []
        merged_list = []
        for i in enhancer_list[0:-1]:
            small_list.append(i) 
            if enhancer_list.index(i) == (len(enhancer_list) -2):
                if i.coord2 >= enhancer_list[enhancer_list.index(i) + 1].coord1:
                    small_list.append(enhancer_list[enhancer_list.index(i) + 1])
                    merged_enhancer = self.merge_enhancers(small_list)
                    merged_list.append(merged_enhancer)
                else:
                    if len(small_list) > 1:
                        merged_enhancer = self.merge_enhancers(small_list)
                        merged_list.append(merged_enhancer)
                    else:
                        merged_list.append(small_list[0])
                    merged_list.append(enhancer_list[enhancer_list.index(i) + 1])
                    small_list = []
            else:
                if i.coord2 >= enhancer_list[enhancer_list.index(i) + 1].coord1:
                    small_list.append(enhancer_list[enhancer_list.index(i) + 1])
                else:
                    if len(small_list) > 1:
                        merged_enhancer = self.merge_enhancers(small_list)
                        merged_list.append(merged_enhancer)
                    else:
                        merged_list.append(small_list[0])
                    small_list = []
        print "merged"
        return merged_list

    def find_TF(self, choose_protein = "Dorsal"):
        #########################################################
        #
        #########################################################
        putative_enhancer_list = []
        print "finding TFs"
        for each_file in self.input_info:
            infile = open(each_file)
            dataset = ""
            method = ""
            protein = ""
            for line in infile:
                if "browser" in line:
                    pass
                elif "track type" in line:
                    line = line.strip()
                    header_info = line.split("=")
                    filename_info = header_info[2].split(" ")
                    dataset = filename_info[0].strip('\"')
                    method = filename_info[1].strip('\"')
                    protein = filename_info[2].strip('\"') 
                else:
                    if protein == choose_protein:
                        found_values = self.get_TF(line, protein, method, dataset)
                        if len(found_values) > 0:
                            for i in found_values:
                               putative_enhancer_list.append(i)                   
            infile.close()
        if len(putative_enhancer_list) > 0: 
            merged_enhancer_list = self.merge_enhancer_list(putative_enhancer_list)
            self.score_enhancers(merged_enhancer_list)
            self.print_bedfile(merged_enhancer_list)

    def score_enhancers(self, enhancer_list):
        #########################################################
        #
        #########################################################
        list_name = self.gene  + "_" + self.status + "_enhancers"
        print "getting info for ", list_name
        for i in enhancer_list:
            for bedgraph in bedfiles:
                i.read_file(bedgraph)
            #for wigfile in phast_files:
            #    i.read_file(wigfile)
            for wigfile in histone_files:
                i.read_file(wigfile)
        my_compare = Enhancer_compare(enhancer_list, list_name)
        Enhancer_compare.print_enhancer_data(my_compare)
