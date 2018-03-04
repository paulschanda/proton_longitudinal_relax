#!/usr/bin/python
    
        
import os             
import sys
#the following line serves just to source the BioPython module, the path where this module is located is defined here. This has to be modified according to the computer used and the correspoding path to biopython
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/biopython-1.44')
from math import cos
import string
from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder
from numarray import *
import random
from math import pi

###########################################################
# Function used to create a directory
def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete        - regular file in the way, raise an exception        - parent directory(ies) does not exist, make them as well    """ 
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        #print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)
###########################################################

def read_input(inputfilename):
	list_specif_deut=[]
	deuteration_list=[]
	deuteration_type_list=[]
        initial_special_list=[]
	initial_special_value_list=[]
	analysis_parameterlist=[]

	file=open(inputfilename,'r')

	field=string.split(file.readline())
	while (field[0].find('#')!=-1):
		field=string.split(file.readline())
	pdbname=str(field[0])
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        freq=float(field[0])
	omega=freq*2.0*3.141592*1e6
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        tauc=float(field[0])*1e-9

	field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        order=float(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        tau=float(field[0])*1e-9
	taui=1/((1/tauc)+(1/tau))
# Deuteration:
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        deut=str(field[0]).lower()
	if deut=='y':  # Case: deuteration
		field=string.split(file.readline())
		while (field[0].find('#')!=-1):
			field=string.split(file.readline())
		frac_non_ex=float(field[0])

                field=string.split(file.readline())
                while (field[0].find('#')!=-1):
                        field=string.split(file.readline())
                frac_ex=float(field[0])

                field=string.split(file.readline())
                while (field[0].find('#')!=-1):
                        field=string.split(file.readline())
                special_deut=str(field[0])

		while True:		
			field=string.split(file.readline())
			if field[0]=='end_deuteration':
				break
			if (field[0].find('#')==-1):
				deut_tupel=(str(field[0]),str(field[1]).upper())
				deuteration_type_list.append(str(field[2]).upper())
				deuteration_list.append(deut_tupel)
		if special_deut=='n':
			deuteration_type_list=[]
			deuteration_list=[]	
	else:
		frac_non_ex=0.0
		frac_ex=0.0
		special_deut='n'
		while True:
			field=string.split(file.readline())
                        if field[0]=='end_deuteration':
				break
# End deuteration
# Initial conditions
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_amide=str(field[0])
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_arom=str(field[0])
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_aliph=str(field[0])
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_methyl=str(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        special_initial=str(field[0]).lower()
	if special_initial=='y':
		while True:
                        field=string.split(file.readline())
                        if field[0]=='end_initial_conditions':
                                break
                        if (field[0].find('#')==-1):
				initial_special=(str(field[0]),str(field[1]).upper())
				initial_special_value_list.append(str(field[2]).lower())
                        	initial_special_list.append(initial_special)
	
        else:
                while True:
                        field=string.split(file.readline())
                        if field[0]=='end_initial_conditions':
                                break
		initial_special_value_list=[]
		initial_special_list=[]

# Simulation details:

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        water_exchange=str(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        water_polarisation=str(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        water_T1=float(field[0])


        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        scantime=float(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        number_scans=float(field[0])
	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        stepsize=float(field[0])
        
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        writeout=str(field[0]).lower()
        
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        write_step=float(field[0])
	############################################
	# Parameters that are analysed
	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        number_of_analyses=int(field[0])

	for counter in range(number_of_analyses):
		parameterlist_local=[]
		nuclei_list=[]

		field=string.split(file.readline())
		while (field[0].find('#')!=-1):
	                field=string.split(file.readline())
        	analysis_type=str(field[0])
		if analysis_type!='e' and analysis_type!='t':
			print "Unknown kind of analysis in input file. The analysis type must be 'e' or 't'. \n\n\n"
		parameterlist_local.append(str(field[0]))
                while True:
                        field=string.split(file.readline())
                        if field[0]=='end_of_list':
                                break
                        if (field[0].find('#')==-1):
				if len(field)>1:
					nucleus=(field[0],field[1])				
				elif (field[0]=='amide' or field[0]=='arom' or field[0]=='methyl' or field[0]=='aliph'):
					nucleus=field[0]
				else:
					print "Unknown nucleus in analysis-list of input file\n\n\n"
				nuclei_list.append(nucleus)	
		parameterlist_local.append(nuclei_list)	
		if analysis_type=='t':
			field=string.split(file.readline())
	                while (field[0].find('#')!=-1):
        	                field=string.split(file.readline())
                	parameterlist_local.append(str(field[0]))
		analysis_parameterlist.append(parameterlist_local)	

	sim_time=number_scans*scantime-stepsize
	return pdbname,omega,tauc,order,taui,deut,frac_non_ex,frac_ex,special_deut,deuteration_list,deuteration_type_list,init_amide,init_arom,init_aliph,init_methyl,special_initial, initial_special_list,initial_special_value_list,water_exchange,water_polarisation,water_T1,sim_time,stepsize,writeout,write_step,analysis_parameterlist,scantime


##############################
def get_pdb_info(pdbname,deut_non_ex,deut_ex,special_deut,deuteration_list,deuteration_type_list):
        # deuteration is the fraction of sites of non-exchangeable hydrogens that are occupied by a deuteron
        # This is the basis for a random rejection of protons to simulate a (partially) deuterated molecule)
	# 
	exchangeable_protons=[('ASN','HD21'),('ASN','HD22'),('GLN','HE21'),('GLN','HE22'),('ARG','HE'),('ARG','HH11'),('ARG','HH12'),('ARG','HH21'),('ARG','HH22'),('CYS','HG'),('HIS','HE2'),('HIS','HD1'),('TYR','HH'),('SER','HG'),('THR','HG1'),('TRP','HE1'),('GLU','HE2'),('ASP','HD2')]


        p=PDBParser(PERMISSIVE=1)
        structure=p.get_structure("example", pdbname)

        counterlist=[]
        sequencelist=[]
        restypelist=[]
        atomnamelist=[]
        x_list=[]
        y_list=[]
        z_list=[]
	exchangeable_protons_indexlist=[]

        counter=-1

        for model in structure.get_list():
                model_id=model.get_id()
                for chain in model.get_list():
                        chain_id=chain.get_id()
                        list_of_residues=chain.get_list()
# LOOP over residues, search all protons and write their coordinates and names

                        for residue in list_of_residues:
                                residue_id=residue.get_id()
                                hetfield, resseq, icode=residue_id
                                current_residue=residue.get_resname()
                                for atom in residue.get_list():
                                        test=atom.get_name()
                                        if test.find('H') !=-1 and test!='NH1' and test!='NH2' and test!='OH':
						typetupel=(current_residue,test)
						numbertupel=(str(resseq),test)
                                        # Deuteration: 
						# first case: specifically deuterated site:
						if (numbertupel in deuteration_list and deuteration_type_list[deuteration_list.index(numbertupel)]=='D') or (typetupel in deuteration_list and deuteration_type_list[deuteration_list.index(typetupel)]=='D'):
							# don't do anything
							print "specifically deuterated site"

						elif (( test!='H' and test!='HN' and typetupel not in exchangeable_protons and random.random()>deut_non_ex ) or ((test=='H' or test=='HN' or typetupel in exchangeable_protons) and random.random()>deut_ex)) or (numbertupel in deuteration_list and deuteration_type_list[deuteration_list.index(numbertupel)]=='H') or (typetupel in deuteration_list and deuteration_type_list[deuteration_list.index(typetupel)]=='H'):

	                                                counter =counter+1
                                                        x,y,z=atom.get_coord()
                                                        restypelist.append(current_residue)
       	                                                counterlist.append(counter)
                                                        sequencelist.append(resseq)
                                                        atomnamelist.append(str(atom.get_name()))
                                                        x_list.append(x)
                                                        y_list.append(y)
                                                        z_list.append(z)

                                                        if (current_residue=='ASN' and test.find('HD') !=-1) or (current_residue=='GLN' and test.find('HE') !=-1) or (current_residue=='ARG' and test.find('HE') !=-1) or (current_residue=='ARG' and test.find('HH') !=-1) or (current_residue=='CYS' and test.find('HG') !=-1) or (current_residue=='HIS' and test=='HE2') or (current_residue=='HIS' and test=='HD1') or (current_residue=='TYR' and test.find('HH') !=-1) or (current_residue=='TYR' and test.find('HG') !=-1) or (current_residue=='SER' and test.find('HG') !=-1) or (current_residue=='THR' and test=='HG1'):
                                                                exchangeable_protons_indexlist.append(counter)


        print " Total number of protons present in the molecule ",len(x_list)
        return counterlist,sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,exchangeable_protons_indexlist

#######################################
def calc_relaxationmatrix(sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,omega,tauc,taui,S2):
	
	####################
	# define constants #
	h=6.62068e-34      #
	pi=3.14159         #
	hbar=h/(pi*2)      #
	gH=2.6751987e8     #
	mu0_4pi=1e-7       #
	####################

	# tauc ... overall tumbling
	# taui ... effective correlation time internal motion
	# S2 ... order parameter methyl groups
	omega2=omega**2  # square of the Larmor frequency
	tauc2=tauc**2    # square of correlation time
	taui2=taui**2	 # square of effective correlation time for internal motion

	# Initialize relaxation matrix
	relax_matrix=resize([0.0],(len(x_list),len(x_list)))

	# Calculate rho*(r**6), i.e. this value has to be multiplied by 1/(r**6) in order to get rho. 
	# This factor is the same irrespective of the distance r
	# Do this for methyl groups and rigid protons, for rho and sigma in a similar way
	rho_factor_ch3=gH**4*hbar**2*mu0_4pi**2*\
	(\
	0.1*tauc*S2+\
	0.1*(1-S2)*taui+\
	0.3*(S2*tauc)/(1+omega2*tauc2)+\
	0.3*((1-S2)*taui)/(1+omega2*taui2)+\
	0.6*(S2*tauc)/(1+4*omega2*tauc2)+\
	0.6*((1-S2)*taui)/(1+4*omega2*taui2))

        rho_factor=gH**4*hbar**2*mu0_4pi**2*\
        (\
        0.1*tauc+\
        0.3*tauc/(1+omega2*tauc2)+\
        0.6*tauc/(1+4*omega2*tauc2))

        sigma_factor_ch3=gH**4*hbar**2*mu0_4pi**2*\
        (\
        -0.1*tauc*S2+\
        -0.1*(1-S2)*taui+\
        0.6*(S2*tauc)/(1+4*omega2*tauc2)+\
        0.6*((1-S2)*taui)/(1+4*omega2*taui2))

        sigma_factor=gH**4*hbar**2*mu0_4pi**2*\
        (\
        -0.1*tauc+\
        0.6*tauc/(1+4*omega2*tauc2))
	
	print " ... Calculating relaxation matrix ... \n"

	for counter in range(len(x_list)):
		sum_rho=0.0
		for counter2 in range(len(x_list)):
			if counter != counter2:
				# check if this is a methyl proton. If so, then use the adapted sigma and rho expressions for the inter-methyl interactions.
				if (((restypelist[counter] == 'MET') and (atomnamelist[counter].find('HE')!=-1)) or ((restypelist[counter] == 'ILE') and (atomnamelist[counter].find('HD1')!=-1)) or ((restypelist[counter] == 'THR') and (atomnamelist[counter].find('HG2')!=-1)) or ( (restypelist[counter] == 'VAL') and (atomnamelist[counter].find('HG')!=-1)) or ((restypelist[counter] == 'ILE') and (atomnamelist[counter].find('HG2')!=-1)) or ((restypelist[counter] == 'LEU') and (atomnamelist[counter].find('HD')!=-1)) or ((restypelist[counter] == 'ALA') and (atomnamelist[counter].find('HB')!=-1))):
					if ((restypelist[counter2] == restypelist[counter]) and ((counter-counter2)**2 <10) and (atomnamelist[counter][1:3]==atomnamelist[counter2][1:3])):
						prefactor_rho=rho_factor_ch3
		                	        prefactor_sigma=sigma_factor_ch3
				# if not intra-methyl
				else:
						prefactor_rho=rho_factor
						prefactor_sigma=sigma_factor
				distance_inv6=1/((1e-10*(((x_list[counter]-x_list[counter2])**2+(y_list[counter]-y_list[counter2])**2+(z_list[counter]-z_list[counter2])**2)**0.5))**6)
				relax_matrix[counter,counter2]=prefactor_sigma*distance_inv6
				sum_rho = sum_rho+prefactor_rho*distance_inv6
		relax_matrix[counter,counter]=sum_rho
	
	print " ... Done ...\n"
		
	return relax_matrix

##########################################
# Calculate the vector of magnetizations reflecting the initial conditions (the Iz - Iz0 terms)
def calc_initial_conditions(sequencelist,restypelist,atomnamelist,init_amide,init_arom,init_aliph,init_methyl,special_initial,initial_special_list,initial_special_value_list):



	unitvector=resize([1.0],(len(atomnamelist),1))		# equilibrium spin magnetization is 1 for each spin
	excitevector=resize([1.0],(len(atomnamelist),1))	
	aromatic_protons=[('PHE','HE1'),('PHE','1HE'),('PHE','HE2'),('PHE','2HE'),('PHE','HD1'),('PHE','1HD'),('PHE','HD2'),('PHE','2HD'),('HIS','HE1'),('HIS','1HE'),('HIS','HD2'),('HIS','2HD'),('TRP','HD1'),('TRP','1HD'),('TRP','HE3'),('TRP','3HE'),('TRP','HZ3'),('TRP','3HZ'),('TRP','HH2'),('TRP','2HH'),('TRP','HZ2'),('TRP','2HZ'),('TYR','HD1'),('TYR','1HD'),('TYR','HE1'),('TYR','1HE'),('TYR','HE2'),('TYR','2HE'),('TYR','HD2'),('TYR','2HD')]
	saturationvector=[]

	for i in range(len(restypelist)):	# create vector reflecting the excitation
						# it has 1 for spins in equilibrium, 0 for 90 deg excitation, -1 for inversion

		
		# Case 1: current atom has individually specified initial conditions
		if ((str(sequencelist[i]),str(atomnamelist[i])) in initial_special_list):
			position=initial_special_list.index((str(sequencelist[i]),str(atomnamelist[i])))
			if str(initial_special_value_list[position])=='s':
				excitevector[i]=0.0
				saturationvector.append(i)
			else:
				excitevector[i]=float(initial_special_value_list[position])
		# Continued Case 1: 
		elif ((str(restypelist[i]),str(atomnamelist[i])) in initial_special_list):
			position=initial_special_list.index((str(restypelist[i]),str(atomnamelist[i])))
                        if str(initial_special_value_list[position])=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(initial_special_value_list[position])

		# Case 2: amide proton
		elif atomnamelist[i] == 'H' or atomnamelist[i] == 'HN':
			if str(init_amide)=='s':
				excitevector[i]=0.0
				saturationvector.append(i)
			else:
				excitevector[i]=float(init_amide)

		# Case 3: aromatic proton
		elif ( (str(restypelist[i]),str(atomnamelist[i])) in aromatic_protons):
			if str(init_arom)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:   
                                excitevector[i]=float(init_arom)
		
		# Case 4: methyl proton
		elif(((restypelist[i] == 'MET') and (atomnamelist[i].find('HE')!=-1)) or ((restypelist[i] == 'ILE') and (atomnamelist[i].find('HD1')!=-1)) or ((restypelist[i] == 'THR') and (atomnamelist[i].find('HG2')!=-1)) or ( (restypelist[i] == 'VAL') and (atomnamelist[i].find('HG')!=-1)) or ((restypelist[i] == 'ILE') and (atomnamelist[i].find('HG2')!=-1)) or ((restypelist[i] == 'LEU') and (atomnamelist[i].find('HD')!=-1)) or ((restypelist[i] == 'ALA') and (atomnamelist[i].find('HB')!=-1))):
			if str(init_methyl)=='s':

				excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_methyl)

		# Last case (5): aliphatic proton
		else:
			if str(init_aliph)=='s':
				excitevector[i]=0.0
                                saturationvector.append(i)
			else:
				excitevector[i] = float(init_aliph)

		#####

	initial_vector=excitevector-unitvector

	return initial_vector,saturationvector
#	return excitevector,saturationvector
def calc_initial_conditions2(sequencelist,restypelist,atomnamelist,init_amide,init_arom,init_aliph,init_methyl,special_initial,initial_special_list,initial_special_value_list):



	unitvector=resize([1.0],(len(atomnamelist),1))		# equilibrium spin magnetization is 1 for each spin
	excitevector=resize([1.0],(len(atomnamelist),1))
	
	aromatic_protons=[('PHE','HE1'),('PHE','1HE'),('PHE','HE2'),('PHE','2HE'),('PHE','HD1'),('PHE','1HD'),('PHE','HD2'),('PHE','2HD'),('HIS','HE1'),('HIS','1HE'),('HIS','HD2'),('HIS','2HD'),('TRP','HD1'),('TRP','1HD'),('TRP','HE3'),('TRP','3HE'),('TRP','HZ3'),('TRP','3HZ'),('TRP','HH2'),('TRP','2HH'),('TRP','HZ2'),('TRP','2HZ'),('TYR','HD1'),('TYR','1HD'),('TYR','HE1'),('TYR','1HE'),('TYR','HE2'),('TYR','2HE'),('TYR','HD2'),('TYR','2HD')]
	saturationvector=[]

	for i in range(len(restypelist)):	# create vector reflecting the excitation
						# it has 1 for spins in equilibrium, 0 for 90 deg excitation, -1 for inversion

		
		# Case 1: current atom has individually specified initial conditions
		if ((str(sequencelist[i]),str(atomnamelist[i])) in initial_special_list):
			position=initial_special_list.index((str(sequencelist[i]),str(atomnamelist[i])))
			if str(initial_special_value_list[position])=='s':
				excitevector[i]=0.0
				saturationvector.append(i)
			else:
				excitevector[i]=cos(float(initial_special_value_list[position])*pi/180.0)
		# Continued Case 1: 
		elif ((str(restypelist[i]),str(atomnamelist[i])) in initial_special_list):
			position=initial_special_list.index((str(restypelist[i]),str(atomnamelist[i])))
                        if str(initial_special_value_list[position])=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=cos(float(initial_special_value_list[position])*pi/180.0)

		# Case 2: amide proton
		elif atomnamelist[i] == 'H' or atomnamelist[i] == 'HN':
			if str(init_amide)=='s':
				excitevector[i]=0.0
				saturationvector.append(i)
			else:
				excitevector[i]=cos(float(init_amide)*pi/180.0)

		# Case 3: aromatic proton
		elif ( (str(restypelist[i]),str(atomnamelist[i])) in aromatic_protons):
			if str(init_arom)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:   
                                excitevector[i]=cos(float(init_arom)*pi/180.0)
		
		# Case 4: methyl proton
		elif(((restypelist[i] == 'MET') and (atomnamelist[i].find('HE')!=-1)) or ((restypelist[i] == 'ILE') and (atomnamelist[i].find('HD1')!=-1)) or ((restypelist[i] == 'THR') and (atomnamelist[i].find('HG2')!=-1)) or ( (restypelist[i] == 'VAL') and (atomnamelist[i].find('HG')!=-1)) or ((restypelist[i] == 'ILE') and (atomnamelist[i].find('HG2')!=-1)) or ((restypelist[i] == 'LEU') and (atomnamelist[i].find('HD')!=-1)) or ((restypelist[i] == 'ALA') and (atomnamelist[i].find('HB')!=-1))):
			if str(init_methyl)=='s':

				excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=cos(float(init_methyl)*pi/180.0)

		# Last case (5): aliphatic proton
		else:
			if str(init_aliph)=='s':
				excitevector[i]=0.0
                                saturationvector.append(i)
			else:
				excitevector[i] = cos(float(init_aliph)*pi/180.0)

		#####

	#initial_vector=excitevector-unitvector

	return excitevector,saturationvector,unitvector
#	return excitevector,saturationvector
#######################################
##############################################################################
# This function calculates the stdev of the values of a list of values, that 
# is given as an input. Additionally the mean value is also an input parameter
def calc_stdev(value_list,avg):
  diff_sqr=0.0
  for j in range (len(value_list)):
   diff_sqr = diff_sqr + (float(value_list[j]) - avg)**2
  stdev=(diff_sqr*(1.0/j))**0.5
  return stdev
##############################################################################

##############################################################################
def calc_avg(value_list):
  sum=0.0
  for j in range (len(value_list)):
   	sum=sum+ (value_list[j])
  average=sum/len(value_list)
  return average 
####################################################################

##############################################################################
# MAIN
##############################################################################


#try:
pdbname,omega,tauc,order,taui,deut,frac_non_ex,frac_ex,special_deut,deuteration_list,deuteration_type_list,init_amide,init_arom,init_aliph,init_methyl,special_initial, initial_special_list,initial_special_value_list,water_exchange,water_polarisation,water_T1,maxtime,timestep,writeout,write_step,analysis_parameterlist,recycledelay=read_input(sys.argv[1])




timesteplist=[i*timestep for i in range((int(maxtime/timestep+0.5))+1)] # list of time steps used in simulation

numberscans=int(maxtime/recycledelay)+1

pulsestep=[k * int(recycledelay/timestep) for k in range(numberscans)]


counterlist,sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,exchangeable_protons_indexlist=get_pdb_info(pdbname,frac_non_ex,frac_ex,special_deut,deuteration_list,deuteration_type_list)

############################################
# Now I have these lists:
# counterlist a counter of all protons
# sequencelist contains sequence position of the atom
# restypelist contains the residue type
# atomnamelist contains atomnames
# initial_special_list: contains tupels of either ('resnumber','atomname') or ('retype','atomname') of atoms which have special initial conditions
# initial_special_value_list: contains the initial conditions of these nuclei
# x_list contains x-coord
# y_list contains y coord
# z_list contains z-coord
# tauc,taui,order ... correlation times, and S2 for CH3 groups
#
# a list of lists containing the desired output-analysis is called analysis_parameterlist [['t', ['amide', ('ARG', 'HA')], '0.63'], ['e', ['methyl', ('588', 'H')]]]
############################################

relaxationmatrix=calc_relaxationmatrix(sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,omega,tauc,taui,order)

#initial_magnetization_vector,saturationvector,magvector=calc_initial_conditions2(sequencelist,restypelist,atomnamelist,init_amide,init_arom,init_aliph,init_methyl,special_initial,initial_special_list,initial_special_value_list)
excitevector,saturationvector,unitvector=calc_initial_conditions2(sequencelist,restypelist,atomnamelist,init_amide,init_arom,init_aliph,init_methyl,special_initial,initial_special_list,initial_special_value_list)

initial_magnetization_vector=excitevector-unitvector
#for i in range(len(x_list)):
#       print counterlist[i],sequencelist[i],restypelist[i],atomnamelist[i],initial_magnetization_vector[i]

#print exchangeable_protons_indexlist," exchangeable"
#print
#print saturationvector," saturated protons"

############################################
# Now I have:
# The relaxation matrix
# the vector reflecting the initial conditions (i.e. the Iz - Iz0 terms for each proton spin)
# saturationvector: contains the number (ie from the counterlist) of the atoms that are kept saturated throughout the simulation
# exchangeable_protons_vector: similar as the saturationvector, this list contains the number (ie from the counterlist) of the atoms that are in fast exchange with water
############################################

vector_of_vector_of_magnetizations=[] # The results of each iterative step are written into this list of lists

# initialize the magnetization vector that represents the state of the spin system
magnetization_vector=initial_magnetization_vector
# and put this initial state as a first list into the list of lists "vector_of_vector_of_magnetizations"
vector_of_vector_of_magnetizations.append(magnetization_vector)



print "\n ... Calculating the time course of the magnetization ...\n"

##########################################
# Calculate the time evolution of the magnetization
##########################################
# The calculation is splitted up in different cases, so that the conditions have to be checked only once, and
# not at every single loop iteration.

counter=0

# Case 1: water exchange ignored, no saturation of protein proton spins
if (water_exchange=='n' and saturationvector==[]):
	
	for time in timesteplist[1:]:   # UNDER CONSTRUCTION HERE IN THIS IF- CASE
		
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)
		else:
			derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
			magnetization_vector=magnetization_vector+derivative*timestep
			vector_of_vector_of_magnetizations.append(magnetization_vector)

# Case 2: water exchange ignored, saturated protein proton spins
elif (water_exchange=='n' and saturationvector!=[]):
	for index in saturationvector:
		magnetization_vector[index]=-1 # polarization minus unit vector, makes -1
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:
			for index in saturationvector:  # set saturated proton polarization to 0
				magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
                	derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)

# Case 3: water is saturated and exchanges with labile protons, no direct saturation of protein proton spins
elif (water_exchange=='y' and water_polarisation=='s' and  saturationvector==[]):
	for index in exchangeable_protons_indexlist:
		magnetization_vector[index]=-1 # polarization minus unit vector, makes -1
	for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)		
		else:	
			for index in exchangeable_protons_indexlist:  # set exchangeable protons to 0 all the time
				magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
                	derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)

# Case 4: water saturated and exchanging, protein protons directly saturated
elif (water_exchange=='y' and water_polarisation=='s' and  saturationvector!=[]):
	for index in exchangeable_protons_indexlist:
                magnetization_vector[index]=-1 # polarization minus unit vector, makes -1
	for index in saturationvector:  # set saturated proton polarization to 0
		magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:		
			for index in exchangeable_protons_indexlist:
	         	       magnetization_vector[index]=-1 # polarization minus unit vector, makes -1
			for index in saturationvector:  # set saturated proton polarization to 0
	        	        magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
			derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)


# Case 5: water exchange of labile protons, water in equilibrium, no direct saturation of protein protons
elif (water_exchange=='y' and float(water_polarisation)==1.0 and  saturationvector==[]): 	
	for index in exchangeable_protons_indexlist:
		magnetization_vector[index]=0
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:		
			for index in exchangeable_protons_indexlist:  # set exchangeable proton polarization to +z
				magnetization_vector[index]=0
                	derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)

# Case 6: water exchange of labile protons, water in equilibrium, direct saturation of protein protons
elif (water_exchange=='y' and float(water_polarisation)==1.0 and  saturationvector!=[]):
	for index in exchangeable_protons_indexlist:
                magnetization_vector[index]=0
	for index in saturationvector:
		magnetization_vector[index]=-1 # polarization minus unit vector, makes -1
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:		
                	for index in saturationvector:  # set saturated proton polarization to 0
                        	magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
			for index in exchangeable_protons_indexlist:  # set exchangeable proton polarization to +z
                        	magnetization_vector[index]=0
			derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)

# Case 7: some general initial condition of the water polarization, no direct rf saturation of protein protons
elif (water_exchange=='y' and saturationvector==[]):
	water_deviation_initial=-1+float(water_polarisation)
	for index in exchangeable_protons_indexlist:
		magnetization_vector[index]=float(water_polarisation)-1.0  
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:		
			water_magn=water_deviation_initial*exp(-time/float(water_T1))
			for index in exchangeable_protons_indexlist:  # set exchangeable proton polarization to the current water-magnetization 
                        	magnetization_vector[index]=water_magn
                	derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
                	magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)


# Case 8: some general initial condition of the water polarization, direct rf saturation of protein protons
elif (water_exchange=='y' and saturationvector!=[]):
        water_deviation_initial=-1+float(water_polarisation)
        for index in exchangeable_protons_indexlist:
                magnetization_vector[index]=float(water_polarisation)-1.0  
	for index in saturationvector:  # set saturated proton polarization to 0
                        magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
        for time in timesteplist[1:]:
		counter=counter+1	
		if counter in pulsestep:
			for element in range(len(magnetization_vector)):
				magnetization_vector[element]=((magnetization_vector[element]+1.0)*excitevector[element])-1.0
			
			vector_of_vector_of_magnetizations.append(magnetization_vector)	
		else:
				
                	water_magn=water_deviation_initial*exp(-time/float(water_T1))
	                for index in exchangeable_protons_indexlist:  # set exchangeable proton polarization to the current water-magnetization
        	                magnetization_vector[index]=water_magn
			for index in saturationvector:  # set saturated proton polarization to 0
                        	magnetization_vector[index]=-1  # polarization minus unit vector, makes -1
	                derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
        	        magnetization_vector=magnetization_vector+derivative*timestep
                	vector_of_vector_of_magnetizations.append(magnetization_vector)

else:
	print "There is a problem. There is no case covering the current settings\n of initialwater-polarization and proton-saturation.\n Check program or input file\n"
 
print "Simulation finished. Writing output files.\n\n"

#######################################################
# GENERATING OUTPUT FILES CONTAINING THE MAGNETIZATION FOR EACH ATOM
if writeout=='y':
	_mkdir('data_per_residue')
	# remove all files which may be present from previous runs
	for datei in os.listdir('data_per_residue'):
		os.remove('data_per_residue/'+str(datei))

# The vector_of_vector_of_magnetizations contains for each timepoint the vector of the magnetization
# The first index refers to the time point, the second to the atom number

	for k in range(len(atomnamelist)):
		filename=file('data_per_residue/'+str(sequencelist[k])+'_'+str(restypelist[k])+'_'+str(atomnamelist[k])+'.dat','w')
		filename.write('#atom number '+str(counterlist[k])+' res '+str(restypelist[k])+' '+str(sequencelist[k])+' '+str(atomnamelist[k])+'\n')

		for t in range(len(timesteplist)):
			if t%(write_step/timestep)==0:
				filename.write(str(timesteplist[t])+'\t'+str((vector_of_vector_of_magnetizations[t][k][0])+1)+'\n')
		filename.close()
#######################################################

###############################################################
# ANALYSIS OF DATA AND OUTPUT

aromatic_protons=[('PHE','HE1'),('PHE','1HE'),('PHE','HE2'),('PHE','2HE'),('PHE','HD1'),('PHE','1HD'),('PHE','HD2'),('PHE','2HD'),('HIS','HE1'),('HIS','1HE'),('HIS','HD2'),('HIS','2HD'),('TRP','HD1'),('TRP','1HD'),('TRP','HE3'),('TRP','3HE'),('TRP','HZ3'),('TRP','3HZ'),('TRP','HH2'),('TRP','2HH'),('TRP','HZ2'),('TRP','2HZ'),('TYR','HD1'),('TYR','1HD'),('TYR','HE1'),('TYR','1HE'),('TYR','HE2'),('TYR','2HE'),('TYR','HD2'),('TYR','2HD')]

methyl_protons=[('ALA','HB1'),('ALA','HB2'),('ALA','HB3'),('ALA','1HB'),('ALA','2HB'),('ALA','3HB'),('MET','HE1'),('MET','HE2'),('MET','HE3'),('MET','1HE'),('MET','2HE'),('MET','3HE'),('ILE','HD11'),('ILE','HD12'),('ILE','HD13'),('ILE','1HD1'),('ILE','2HD1'),('ILE','3HD1'),('THR','HG21'),('THR','HG22'),('THR','HG23'),('THR','1HG2'),('THR','2HG2'),('THR','3HG2'),('VAL','1HG1'),('VAL','2HG1'),('VAL','3HG1'),('VAL','HG11'),('VAL','HG12'),('VAL','HG13'),('ILE','HG21'),('ILE','HG22'),('ILE','HG23'),('ILE','1HG2'),('ILE','2HG2'),('ILE','3HG2'),('LEU','HD11'),('LEU','HD12'),('LEU','HD13'),('LEU','1HD1'),('LEU','2HD1'),('LEU','3HD1'),('LEU','HD21'),('LEU','HD22'),('LEU','HD23'),('LEU','1HD2'),('LEU','2HD2'),('LEU','3HD2'),('VAL','1HG2'),('VAL','2HG2'),('VAL','3HG2'),('VAL','HG21'),('VAL','HG22'),('VAL','HG23')]


if len(analysis_parameterlist) > 0:
	try:
		os.remove('simulation_results.dat')
	except:
		print
	filename_analyses=file('simulation_results.dat','w')
	filename_analyses.write('# SIMULATION RESULTS, BASED ON INPUT PARAMETERS SPECIFIED IN FILE: '+str(sys.argv[1]))

# Iteration over all the different analyses
for i in range(len(analysis_parameterlist)):
	
	filename_analyses.write('\n\n# Results of Analysis number '+str(i+1)+'\n\n')

	########################
	# Case 1: This is an analysis at which timepoint a certain magnetization value is reached, on average
	if analysis_parameterlist[i][0]=='t':
		timepoint_list=[]
		cutofflevel=float(analysis_parameterlist[i][2])
		filename_analyses.write('# Determining the time points where the magnetization reaches a level of '+str(cutofflevel)+'\n\n')

		# Check all atomic nuclei
		for k in range(len(atomnamelist)):
			counter=0

			# Check if the current nucleus is of interest (ie if it is in the specified list of nuclei to be analysed)
			if ((restypelist[k],atomnamelist[k]) in analysis_parameterlist[i][1] or (str(sequencelist[k]),atomnamelist[k]) in analysis_parameterlist[i][1] or (('amide' in analysis_parameterlist[i][1]) and (atomnamelist[k] =='H' or atomnamelist[k] =='HN')) or (('arom' in analysis_parameterlist[i][1]) and (restypelist[k],atomnamelist[k]) in aromatic_protons) or (('methyl' in analysis_parameterlist[i][1]) and (restypelist[k],atomnamelist[k]) in methyl_protons) or (('aliph' in analysis_parameterlist[i][1]) and (atomnamelist[k] !='H' or atomnamelist[k] !='HN') and ((restypelist[k],atomnamelist[k]) not in aromatic_protons) and ((restypelist[k],atomnamelist[k]) not in methyl_protons ))  ):

				for t in range(len(timesteplist)):
				
					if ((vector_of_vector_of_magnetizations[t][k][0])+1)>cutofflevel:
						timepoint_list.append(timesteplist[t])	# Now I found the effective T1 and put it into list.
						filename_analyses.write(str(restypelist[k])+'\t'+str(sequencelist[k])+'\t'+str(atomnamelist[k])+'\t'+str(timesteplist[t])+'\n')
						break
					if counter==len(timesteplist)-1 :
						filename_analyses.write(str(restypelist[k])+'\t'+str(sequencelist[k])+'\t'+str(atomnamelist[k])+'\t DID NOT REACH TARGET LEVEL, NOT IN AVERAGE!!\n')
						break
					counter = counter + 1					
		try:
			average_time=calc_avg(timepoint_list)
			stdev=calc_stdev(timepoint_list,average_time)
			filename_analyses.write('\n\t Average time point '+str(average_time)+'\n\t Standard deviation '+str(stdev))	
		except:
			filename_analyses.write('\n\t Averaging not possible because there are not enough data points')

	
	########################
	# Case 2: This is an average of the magnetization at the end of the simulation
	if analysis_parameterlist[i][0]=='e':
		filename_analyses.write('\n# Determining the magnetization present at the end of the simulation\n\n')
		
		magnetization_end_list=[]
                # Check all atomic nuclei
                for k in range(len(atomnamelist)):	
			# Check if the current nucleus is of interest (ie if it is in the specified list of nuclei to be analysed)
			if ((restypelist[k],atomnamelist[k]) in analysis_parameterlist[i][1] or (str(sequencelist[k]),atomnamelist[k]) in analysis_parameterlist[i][1] or (('amide' in analysis_parameterlist[i][1]) and (atomnamelist[k] =='H' or atomnamelist[k] =='HN')) or (('arom' in analysis_parameterlist[i][1]) and (restypelist[k],atomnamelist[k]) in aromatic_protons) or (('methyl' in analysis_parameterlist[i][1]) and (restypelist[k],atomnamelist[k]) in methyl_protons) or (('aliph' in analysis_parameterlist[i][1]) and not (atomnamelist[k] =='H' or atomnamelist[k] =='HN' or ((restypelist[k],atomnamelist[k]) in aromatic_protons) or ((restypelist[k],atomnamelist[k])  in methyl_protons ))))  :


				magnetization_end_list.append(vector_of_vector_of_magnetizations[-1][k][0]+1)				
				filename_analyses.write(str(restypelist[k])+'\t'+str(sequencelist[k])+'\t'+str(atomnamelist[k])+'\t'+str(vector_of_vector_of_magnetizations[-1][k][0]+1)+'\n')
		average_time=calc_avg(magnetization_end_list)
		stdev=calc_stdev(magnetization_end_list,average_time)
		filename_analyses.write('\n\t Average end magnetization '+str(average_time)+'\n\t Standard deviation '+str(stdev))

print "\n\n   Done.\n\n"
