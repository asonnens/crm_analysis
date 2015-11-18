#chip analysis classes
#this should contain classes for analyzing chip peaks 
#currently just contains a ChipPeak class
#includes default constructor, a create_peak constructor, and a return peak coordinates 


#ChipPeak class 
class ChipPeaks:
	#default constructor (name of peak(1,2,3, etc.), TF, coordinate min/max, and peak value
    def_init_(self):
	    self.name = ""
		self.TF = ""
		self.coords = [0,0]
		self.peak_val = 0
		
	#create a new peak-- takes arguements name, min and max, and peak value
	def_create_peak(self, name, TF, min, max, val):
	    self.name = name
		self.TF = TF
		self.coords = [min, max]
		self.peak_val = val
	
	#gets coordinates given a name I'm not sure if this is correct?
	def_get_coords(self, name):
	    return self.coords[:]

#ChipPeak dict
class ChipPeakDict:

"""
chip peak dictionary class
-put chip peaks in a dictionary organized by 'peak.name', 'peak'
-constructor/initializer
                