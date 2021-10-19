###################################
# script takes in series of
# ECAL hits in order ieta, iphi, hitE,
# and returns .csv file of dimensions
# ieta x iphi, with mapped hitE with 
# converted to # of lego bricks per 
# based on hitE.
#
# to-do:
# -create criteria for # of lego bricks
# based on hitE
# -figure out how to input .txt file
# -take .txt file and put into 2d list
#    of floats (?) 
# -nested for loop outputting commas/
#    # of lego bricks,                 
# -save output .csv file              
######################################

def separate(array):
	output = []
	for i in range(len(array)):
		output.append(array[i].split(","))
	return output

def removeFromEnd(array): #removes last two chars (\n in this case) from end of the 3rd element 
	for i in range(len(array)):                                                   #of each array
		tobesplit = list(array[i][2])
		tobesplit = tobesplit[:-1] #removes last element of list
		switch = ''.join(tobesplit)
		array[i][2] = switch
	return array	

def gather(array, index):
	output = []
	for i in range(len(array)):
		output.append(array[i][index])
	return output

def colorgradiant(hitE, isECAL):
	if isECAL:		
		if(hitE < .24): #hello future physics intern making lego model 4.0, I would suggest making this number .23
			return 1 #so when you are creating seeds for clusters you can go off of their color of where or not
		elif(hitE < .48): #they meet the threshhold or not
			return 2
		elif(hitE < .95):
			return 3
		elif(hitE < 1.43):
			return 4
		elif(hitE < 1.9):
			return 5
		elif(hitE < 2.86):
			return 6
		elif(hitE < 10):
			return 7
		elif(hitE < 40):
			return 8
		else:
			return 9
	else:
		if(hitE < 1.6):
			return 1
		elif(hitE < 7):
			return 2
		else:
			return 3

def maplego(array, isECAL):
	ieta = gather(array, 0)
	iphi = gather(array, 1)
	hitE = gather(array, 2)
	if isECAL:
		eta = 90
		phi = 380
		offset = 40
	else:
		eta = 16
		phi = 72
		offset = 8
	counter = 0
	legos = 0;
	x = [[',' for i in range(eta)] for j in range(phi)]
	for k in range(len(array)):
		#x[int(iphi[k])-1][int(ieta[k])+offset] = str(colorgradiant(float(hitE[k]), isECAL))+","
		x[int(iphi[k])-1][int(ieta[k])+offset] = ieta[k]+";"+iphi[k]+";"+hitE[k]+","
		counter += 1
		legos += colorgradiant(float(hitE[k]), isECAL)
	print "Hits: " + str(counter)
	if isECAL:
		print "Red Legos: " + str(legos)
	else:
		print "Blue Legos: " + str(legos*5)
	return x

###################### MAIN ####################################################################################
		
isECAL = False #Change based on if want to output csv file for hcal or ecal

if isECAL:
	ECAL = open("ehits.txt", "r") #open(txt_file_name, mode) modes: 'r' = read only, 'w' = write, 'a' = append
	output = open("ECAL.csv", "w+") # 'a+' indicates append and create ('+') file if needed
	lines = ECAL.readlines()
	eta = 90
else:
	HCAL = open("hhits.txt", "r")
	output = open("HCAL.csv", "w+") # 'a+' indicates append and create ('+') file if needed
	lines = HCAL.readlines()
	eta = 16 

sep = separate(lines)
clean = removeFromEnd(sep) # clean: array containing arrays in order [ieta, iphi, hitE] 

print(clean)

legomap = maplego(clean, isECAL) # legomap: 2d array with 


for i in range(len(legomap)):
	for ieta in range(eta):
	     output.write(str(legomap[i][ieta]))
	output.write("\n")
output.close() 
