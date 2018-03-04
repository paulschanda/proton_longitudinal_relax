#!/usr/bin/python

import os       
import sys
#the following line serves just to source the BioPython module, the path where this module is located is defined here. This has to be modified according to the computer used and the correspoding path to biopython
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/biopython-1.44')
        
import string
from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder
from numarray import *
import random
import shutil

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
############################################################################################
# This function reads the input file
def read_input(inputfilename):

        initial_special_list=[]
	initial_special_value_list=[]
	analysis_parameterlist=[]
	exchange_list=[]	# This list takes tupels of (res-type,atom-type) which exchange with solvent (this means that this atom-type behaves the same in all residues of the same type) or it can also take tupels like (res-number,atom-type) in case a specific atom in a specific residue (eg Uracil 3, H1) does something different from the other Uracil-H1 atoms.
	exchange_rate_list=[]	# This list takes the exchange rates of the protons in the above list; in the same order
	file=open(inputfilename,'r')

	field=string.split(file.readline())
# Read pdb name
	while (field[0].find('#')!=-1):
		field=string.split(file.readline())
	pdbname=str(field[0])
# Read spectrometer frequency and calculate Larmor frequency (specfreq times 2pi)
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        freq=float(field[0])
	omega=freq*2.0*3.141592*1e6
# Read the overall tumblling correlation time
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        tauc=float(field[0])*1e-9
##
# Exchange with water
# Take into account solvent exchange ?
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        exch=str(field[0]).lower()
	if exch=='y':  # Case:  there is solvent exchange of labile protons
		# Go through the whole list of protons that can exchange:
		while True:		
			field=string.split(file.readline())
			if field[0]=='end_exchange':
				break
			if (field[0].find('#')==-1):
				exch_tupel=(str(field[0]),str(field[1]).upper())
				exchange_rate_list.append(float(field[2]))
				exchange_list.append(exch_tupel)
	else:
		while True:
			field=string.split(file.readline())
                        if field[0]=='end_exchange':
				break
# End exchange

# Initial conditions
	# H5'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H5prime=str(field[0])
	# H5''	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H5dprime=str(field[0])
	# H4'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H4prime=str(field[0])
	# H3'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H3prime=str(field[0])
	# H2'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H2prime=str(field[0])
	# H1'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_H1prime=str(field[0])
	# HO2'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_HO2prime=str(field[0])
	# HO3'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_HO3prime=str(field[0])
	# HO5'	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_HO5prime=str(field[0])
	# A H8	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_A_H8=str(field[0])
	# A H61	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_A_H61=str(field[0])
	# A H62	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_A_H62=str(field[0])
	# A H2	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_A_H2=str(field[0])
	# G H8	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_G_H8=str(field[0])
	# G H1	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_G_H1=str(field[0])
	# G H21	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_G_H21=str(field[0])
	# G H22
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_G_H22=str(field[0])
	# C H41
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_C_H41=str(field[0])
	# C H42
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_C_H42=str(field[0])
	# C H5	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_C_H5=str(field[0])
	# C H6	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_C_H6=str(field[0])
	# U H3 
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_U_H3=str(field[0])
	# U H5	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_U_H5=str(field[0])
	# U H6	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        init_U_H6=str(field[0])

# special treatment for selected atoms ?
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
        number_of_sims=int(field[0])
	
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        water_polarisation=float(field[0])

        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        water_T1=float(field[0])
# Simulation duration in seconds
        field=string.split(file.readline())
        while (field[0].find('#')!=-1):
                field=string.split(file.readline())
        sim_time=float(field[0])

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
				elif (field[0]=='imino' or field[0]=='aromatic' or field[0]=='amino' or field[0]=='sugar'):
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


	return pdbname,omega,tauc,exch,exchange_list,exchange_rate_list,init_H5prime,init_H5dprime,init_H4prime,init_H3prime,init_H2prime,init_H1prime,init_HO2prime,init_HO3prime,init_HO5prime,init_A_H8,init_A_H61,init_A_H62,init_A_H2,init_G_H8,init_G_H1,init_G_H21,init_G_H22,init_C_H41,init_C_H42,init_C_H5,init_C_H6,init_U_H3,init_U_H5,init_U_H6,special_initial, initial_special_list,initial_special_value_list,number_of_sims,water_polarisation,water_T1,sim_time,stepsize,writeout,write_step,analysis_parameterlist


############################################################################################
# This function reads the pdb and returns a couple of lists containing different properties of the H-atoms (their coordinates, types, etc)
# It makes use of the BioPython package (actually, this is the only function where BioPython is required)
def get_pdb_info(pdbname):
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

			for residue in list_of_residues:	#residue is for each res an expression like this: <Residue   A het=  resseq=1 icode= >
                                residue_id=residue.get_id()	# residue_id is an expression that looks like this: (' ', 1, ' ')
                                hetfield, resseq, icode=residue_id	# resseq is the sequence position of the residue
                                current_residue=residue.get_resname()	# current_residue is the name of the res (eg U or A)
                                for atom in residue.get_list():		# go through all the atoms
                                        test=atom.get_name()		# write the atom name (eg H4' into the variable 'test'
					if test.find('H') !=-1:		# go only through the H-atoms
						counter =counter+1
                                                x,y,z=atom.get_coord()
                                                restypelist.append(current_residue.lstrip())
       	                                        counterlist.append(counter)
                                                sequencelist.append(resseq)
                                                atomnamelist.append(str(atom.get_name()))
                                                x_list.append(x)
                                                y_list.append(y)
                                                z_list.append(z)
						# write the exchangeable protons into a list. See here which protons I assume exchangeable:
						if (current_residue.find('U')!=-1 and test=='H3') or (current_residue.find('C')!=-1 and test=='H41') or (current_residue.find('C')!=-1 and test=='H42') or (current_residue.find('G')!=-1 and test=='H1') or (current_residue.find('G')!=-1 and test=='H21') or (current_residue.find('G')!=-1 and test=='H22') or (current_residue.find('A')!=-1 and test=='H61') or (current_residue.find('A')!=-1 and test=='H62') or test.find('HO2')!=-1 or test.find('HO5')!=-1 or test.find('HO3')!=-1:
                             				exchangeable_protons_indexlist.append(counter)
	return counterlist,sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,exchangeable_protons_indexlist
#######################################################
def calc_relaxationmatrix(sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,omega,tauc):
	
	####################
	# define constants #
	h=6.62068e-34      #
	pi=3.14159         #
	hbar=h/(pi*2)      #
	gH=2.6751987e8     #
	mu0_4pi=1e-7       #
	####################

	# tauc ... overall tumbling
	omega2=omega**2  # square of the Larmor frequency
	tauc2=tauc**2    # square of correlation time

	# Initialize relaxation matrix
	relax_matrix=resize([0.0],(len(x_list),len(x_list)))

	# Calculate rho*(r**6), i.e. this value has to be multiplied by 1/(r**6) in order to get rho. 
	# This factor is the same irrespective of the distance r
	# Note that there is no internal dynamics at all. spectral densities are simply the rigid ones with tumbling.

        prefactor_rho=gH**4*hbar**2*mu0_4pi**2*\
        (\
        0.1*tauc+\
        0.3*tauc/(1+omega2*tauc2)+\
        0.6*tauc/(1+4*omega2*tauc2))

        prefactor_sigma=gH**4*hbar**2*mu0_4pi**2*\
        (\
        -0.1*tauc+\
        0.6*tauc/(1+4*omega2*tauc2))
	
	print " ... Calculating relaxation matrix ... \n"

	for counter in range(len(x_list)):
		sum_rho=0.0
		for counter2 in range(len(x_list)):
			if counter != counter2:
				distance_inv6=1/((1e-10*(((x_list[counter]-x_list[counter2])**2+(y_list[counter]-y_list[counter2])**2+(z_list[counter]-z_list[counter2])**2)**0.5))**6)
				relax_matrix[counter,counter2]=prefactor_sigma*distance_inv6
				sum_rho = sum_rho+prefactor_rho*distance_inv6
		relax_matrix[counter,counter]=sum_rho
	
	print " ... Done ...\n"
		
	return relax_matrix

#######################################################
# Calculate the vector of magnetizations reflecting the initial conditions (the Iz - Iz0 terms)
def calc_initial_conditions(sequencelist,restypelist,atomnamelist,init_H5prime,init_H5dprime,init_H4prime,init_H3prime,init_H2prime,init_H1prime,init_HO2prime,init_HO3prime,init_HO5prime,init_A_H8,init_A_H61,init_A_H62,init_A_H2,init_G_H8,init_G_H1,init_G_H21,init_G_H22,init_C_H41,init_C_H42,init_C_H5,init_C_H6,init_U_H3,init_U_H5,init_U_H6,special_initial,initial_special_list,initial_special_value_list):
	unitvector=resize([1.0],(len(atomnamelist),1))		# equilibrium spin magnetization is 1 for each spin
	excitevector=resize([1.0],(len(atomnamelist),1))	
	saturationvector=[]

	for i in range(len(restypelist)):	# create vector reflecting the excitation
						# it has 1 for spins in equilibrium, 0 for 90 deg excitation, -1 for inversion

		# CASE 1: ATOM THAT HAS INDIVIDUALLY SPECIFIED INITIAL CONDITIONS		
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

		# Case 2: H5'
		elif atomnamelist[i] == "H5'":
			if str(init_H5prime)=='s':
				excitevector[i]=0.0
				saturationvector.append(i)
			else:
				excitevector[i]=float(init_H5prime)
                # Case 3: H5''
                elif atomnamelist[i] == "H5''":
                        if str(init_H5dprime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_H5dprime)

                # Case 4: H4'
                elif atomnamelist[i] == "H4'":
                        if str(init_H4prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_H4prime)

                # Case 5: H3'
                elif atomnamelist[i] == "H3'":
                        if str(init_H3prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_H3prime)

                # Case 6: H2'
                elif atomnamelist[i] == "H2'":
                        if str(init_H2prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_H2prime)

                # Case 7: H1'
                elif atomnamelist[i] == "H1'":
                        if str(init_H1prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_H1prime)

                # Case 8: HO2'
                elif atomnamelist[i] == "HO2'":
                        if str(init_HO2prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_HO2prime)

                # Case 9: HO3'
                elif atomnamelist[i] == "HO3'":
                        if str(init_HO3prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_HO3prime)

                # Case 10: HO5'
                elif atomnamelist[i] == "HO5'":
                        if str(init_HO5prime)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_HO5prime)

                # Case 11: A H8
                elif atomnamelist[i] == 'H8' and restypelist[i].find('A')!=-1:
                        if str(init_A_H8)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_A_H8)

                # Case 12: A H61
                elif atomnamelist[i] == 'H61' and restypelist[i].find('A')!=-1:
                        if str(init_A_H61)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_A_H61)

                # Case 13: A H62
                elif atomnamelist[i] == 'H62' and restypelist[i].find('A')!=-1:
                        if str(init_A_H62)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_A_H62)

                # Case 14: A H2
                elif atomnamelist[i] == 'H2' and restypelist[i].find('A')!=-1:
                        if str(init_A_H2)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_A_H2)

                # Case 15: G H8
                elif atomnamelist[i] == 'H8' and restypelist[i].find('G')!=-1:
                        if str(init_G_H8)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_G_H8)

                # Case 16: G H1
                elif atomnamelist[i] == 'H1' and restypelist[i].find('G')!=-1:
                        if str(init_G_H1)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_G_H1)

                # Case 17: G H21
                elif atomnamelist[i] == 'H21' and restypelist[i].find('G')!=-1:
                        if str(init_G_H21)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_G_H21)

                # Case 18: G H22
                elif atomnamelist[i] == 'H22' and restypelist[i].find('G')!=-1:
                        if str(init_G_H22)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_G_H22)

                # Case 19: C H41
                elif atomnamelist[i] == 'H41' and restypelist[i].find('C')!=-1:
                        if str(init_C_H41)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_C_H41)

                # Case 20: C H42
                elif atomnamelist[i] == 'H42' and restypelist[i].find('C')!=-1:
                        if str(init_C_H42)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_C_H42)
                # Case 21: C H5
                elif atomnamelist[i] == 'H5' and restypelist[i].find('C')!=-1:
                        if str(init_C_H5)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_C_H5)

                # Case 22: C H6
                elif atomnamelist[i] == 'H6' and restypelist[i].find('C')!=-1:
                        if str(init_C_H6)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_C_H6)

                # Case 23: U H3
                elif atomnamelist[i] == 'H3' and restypelist[i].find('U')!=-1:
                        if str(init_U_H3)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_U_H3)

                # Case 24: U H5
                elif atomnamelist[i] == 'H5' and restypelist[i].find('U')!=-1:
                        if str(init_U_H5)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_U_H5)

                # Case 25: U H6
                elif atomnamelist[i] == 'H6' and restypelist[i].find('U')!=-1:
                        if str(init_U_H6)=='s':
                                excitevector[i]=0.0
                                saturationvector.append(i)
                        else:
                                excitevector[i]=float(init_U_H6)

		else:
			print "atom type not found\n"

		#####

	initial_vector=excitevector-unitvector

	return initial_vector,saturationvector
#	return excitevector,saturationvector
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


#######################################################

# MAIN
try:
	inputfile=sys.argv[1]
	outdir=sys.argv[2]
except:
	print " Usage:\n Use the input file as argument one\n and the name of the directory where\n the data should be stored as argument two"
	sys.exit(0)

# Write the input file into the output directory
_mkdir(outdir)
shutil.copyfile(inputfile,outdir+'/inputfile')
	
pdbname,omega,tauc,exch,exchange_list,exchange_rate_list,init_H5prime,init_H5dprime,init_H4prime,init_H3prime,init_H2prime,init_H1prime,init_HO2prime,init_HO3prime,init_HO5prime,init_A_H8,init_A_H61,init_A_H62,init_A_H2,init_G_H8,init_G_H1,init_G_H21,init_G_H22,init_C_H41,init_C_H42,init_C_H5,init_C_H6,init_U_H3,init_U_H5,init_U_H6,special_initial, initial_special_list,initial_special_value_list,number_of_sims,water_polarisation,water_T1,sim_time,stepsize,writeout,write_step,analysis_parameterlist=read_input(sys.argv[1])

# Treatment of water exchange: each exchangeable proton will during the simulation cross-relax with the others, and in addition, after each simulation step (relaxation)
# it has a certain probability to exchange with water. If it exchanges, then the new polarisation of this site is simply the water polarisation, otherwise it will
# carry the polarisation it had at the end of the simulation step on to the next simulation step. The probability of exchange is given by the exchange rate (in seconds ^-1)
# normalized to the size of the simulation steps.
# Exchange rates: Calculate for each exchanging atom (in the exchange_list) the probability at each simulation step to exchange with water.
# Have to read in the exchange rates (in sec^-1) from the exchange_rate_list, and multiply with stepsize
exch_probab_list=[]
for i in range(len(exchange_rate_list)):
	prob=exchange_rate_list[i]*stepsize
	if prob>0.9999:		# if exchange rate is faster than the simulation step size
		prob=0.9999
	exch_probab_list.append(prob)
	
# LIST OF TIMESTEPS OF THE SIMULATION
timesteplist=[i*stepsize for i in range((int(sim_time/stepsize+0.5))+1)] # list of time steps used in simulation



#print pdbname,omega,tauc,exch,exchange_list,exchange_rate_list,init_H5prime,init_H5dprime,init_H4prime,init_H3prime,init_H2prime,init_H1prime,init_HO2prime,init_HO3prime,init_HO5prime,init_HO3prime,init_A_H8,init_A_H61,init_A_H62,init_A_H2,init_G_H8,init_G_H1,init_G_H21,init_G_H22,init_C_H41,init_C_H42,init_C_H5,init_C_H6,init_U_H3,init_U_H5,init_U_H6,special_initial, initial_special_list,initial_special_value_list,water_polarisation,water_T1,sim_time,stepsize,writeout,write_step,analysis_parameterlist

#### READ PDB FILE AND GET ATOM COORDINATES
print " Reading pdb file\n"
counterlist,sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,exchangeable_protons_indexlist = get_pdb_info(pdbname)
print " Finished. Number of protons found: ",len(atomnamelist)

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
# tauc ... correlation time overall tumbling
# a list of lists containing the desired output-analysis is called analysis_parameterlist [['t', [('U', 'H3'), ('G', 'H1')], '0.63']]
# exch_probab_list: list containing for all the exchangeable protons their per-stepsiez-probability to get exchange with water
############################################

print " Calculating the relaxation matrix\n\n"

# CALCULATE THE RELAXATION MATRIX
relaxationmatrix=calc_relaxationmatrix(sequencelist,restypelist,atomnamelist,x_list,y_list,z_list,omega,tauc)

# MAKE THE VECTOR OF INITIAL Iz VALUES
initial_magnetization_vector,saturationvector=calc_initial_conditions(sequencelist,restypelist,atomnamelist,init_H5prime,init_H5dprime,init_H4prime,init_H3prime,init_H2prime,init_H1prime,init_HO2prime,init_HO3prime,init_HO5prime,init_A_H8,init_A_H61,init_A_H62,init_A_H2,init_G_H8,init_G_H1,init_G_H21,init_G_H22,init_C_H41,init_C_H42,init_C_H5,init_C_H6,init_U_H3,init_U_H5,init_U_H6,special_initial,initial_special_list,initial_special_value_list)


# EXCHANGE WITH SOLVENT: SETTING UP THE LISTS THAT WILL BE USED DURING THE SIMULATION
# Set up a list that contains a number for each proton. This number will then later be set to the
# per-simulation-step-time probability that this atom exchanges with water
# This list is as long as the list that contains the atoms, of course, and is ordered in the same way
exchange_probability_vector=[0.0*i for i in range(len(initial_magnetization_vector))]

for i in (exchangeable_protons_indexlist):
	tupelres=(restypelist[i].rstrip(),atomnamelist[i])
	for count in range(len(exchange_list)):
		if exchange_list[count] == tupelres:	# The current atom type is generally exchangeable
			prob = exch_probab_list[count]
			exchange_probability_vector[i]=prob
# Now go again through the exchangeable_protons_indexlist and search for specific assignments of ex rates
# which differ from the global exchange rates for this proton type. (eg if HO2' of res 3 behaves differently
# than the HO2' assumed generally
for i in (exchangeable_protons_indexlist):
	tupelseq=(str(sequencelist[i]),atomnamelist[i])
	for count in range(len(exchange_list)):
		if exchange_list[count] == tupelseq:	# The current atom type is specifically exchangeable
			prob = exch_probab_list[count]
			exchange_probability_vector[i]=prob
#for i in range(len(exchange_probability_vector)):
#	print sequencelist[i],restypelist[i],atomnamelist[i],exchange_probability_vector[i]
#### END OF SETUP OF SOLVENT EXCHANGE LIST

############################################
# Now I have:
# The relaxation matrix
# the vector reflecting the initial conditions (i.e. the Iz - Iz0 terms for each proton spin)
# saturationvector: contains the number (ie from the counterlist) of the atoms that are kept saturated throughout the simulation
# exchangeable_protons_indexlist: similar as the saturationvector, this list contains the number (ie from the counterlist) of the atoms that are in fast exchange with water
# exchange_probability_vector: tells for each exchangeable atom what is the probab that this atom exchanges in a simulation time step
############################################


print "\n ... Calculating the time course of the magnetization ...\n"

superlist=[]	# This list will be filled with the simulation results from the individual simulations

##########################################
# Calculate the time evolution of the magnetization
##########################################

# Do the calculations number_of_sims times
for simulationscounter in range(number_of_sims):
	print " Running simulation number", simulationscounter
	# Initialization of magnetization vectors
	vector_of_vector_of_magnetizations=[] # The results of each iterative step are written into this list of lists
	# initialize the magnetization vector that represents the state of the spin system
	magnetization_vector=initial_magnetization_vector
	# and put this initial state as a first list into the list of lists "vector_of_vector_of_magnetizations"
	vector_of_vector_of_magnetizations.append(magnetization_vector)

	if (exch=='y'  and  saturationvector==[]): 	# Maybe I later want to implement different schemes, with saturation ...
		for time in timesteplist[1:]:
			current_water=(water_polarisation-1.0)*exp(-time/water_T1)
			for index in exchangeable_protons_indexlist:  # set exchangeable proton polarization to the water polarization
				test=random.random()
				if random.random() < exchange_probability_vector[index]:	# Event: proton gets exchanged with relaxed water
					magnetization_vector[index]=current_water
			derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
			magnetization_vector=magnetization_vector+derivative*stepsize
			vector_of_vector_of_magnetizations.append(magnetization_vector)
	
	elif (exch=='n'  and  saturationvector==[]):
		for time in timesteplist[1:]:
			derivative=(-1)*matrixmultiply(relaxationmatrix,magnetization_vector)
			magnetization_vector=magnetization_vector+derivative*stepsize
			vector_of_vector_of_magnetizations.append(magnetization_vector)
	
		
	else:					# Maybe I use this later
		print " Case not implmented yet\n"
		sys.exit(0)
		
	superlist.append(vector_of_vector_of_magnetizations)
print "  Simulation finished. Averaging over the simulations.\n\n"


averaged_vector_of_vector_of_magnetizations=[]


# Average over the simulations
for timestep in range(len(superlist[0])):
	proton_list=[]
	for atomnumber in range(len(superlist[0][0])):		# Go through the different atoms
		simaveragelist=[]
		for simnumber in range(len(superlist)):		# Go through the different simulations
			simaveragelist.append(superlist[simnumber][timestep][atomnumber][0])
			#print superlist[0][timestep][atomnumber],superlist[1][timestep][atomnumber]
		#print simaveragelist
		atomaverage=calc_avg(simaveragelist)
		#print simaveragelist,atomaverage
		proton_list.append(atomaverage)
	#print len(proton_list)
	averaged_vector_of_vector_of_magnetizations.append(proton_list)

print " Creating output files"
#print len(averaged_vector_of_vector_of_magnetizations),len(averaged_vector_of_vector_of_magnetizations[0])
# 
#######################################################
# GENERATING OUTPUT FILES CONTAINING THE MAGNETIZATION FOR EACH ATOM
if writeout=='y':
	
	# remove all files which may be present from previous runs
	for datei in os.listdir(outdir):
		os.remove(str(outdir)+'/'+str(datei))

# The vector_of_vector_of_magnetizations contains for each timepoint the vector of the magnetization
# The first index refers to the time point, the second to the atom number

	for k in range(len(atomnamelist)):
		filename=file(str(outdir)+'/'+str(sequencelist[k])+'_'+str(restypelist[k])+'_'+str(atomnamelist[k])+'.dat','w')
		filename.write('#atom number '+str(counterlist[k])+' res '+str(restypelist[k])+' '+str(sequencelist[k])+' '+str(atomnamelist[k])+'\n')

		for t in range(len(timesteplist)):
			if t%(write_step/stepsize)==0:
				filename.write(str(timesteplist[t])+'\t'+str((averaged_vector_of_vector_of_magnetizations[t][k])+1)+'\n')
		filename.close()
#######################################################

###############################################################
# ANALYSIS OF DATA AND OUTPUT

aromatic_protons=[('U','H5'),('U','H6'),('C','H5'),('C','H6'),('A','H8'),('A','H2'),('G','H8')]
imino_protons=[('G','H1'),('U','H3')]
amino_protons=[('A','H61'),('A','H62'),('C','H41'),('C','H42')]
sugar_protons=[('A',"H1'"),('A',"H2'"),('A',"H3'"),('A',"H4'"),('A',"H5'"),('A',"H5''"),('A',"HO2'"),('A',"HO3'"),('A',"HO5'"), ('G',"H1'"),('G',"H2'"),('G',"H3'"),('G',"H4'"),('G',"H5'"),('G',"H5''"),('G',"HO2'"),('G',"HO3'"),('G',"HO5'"), ('C',"H1'"),('C',"H2'"),('C',"H3'"),('C',"H4'"),('C',"H5'"),('C',"H5''"),('C',"HO2'"),('C',"HO3'"),('C',"HO5'"), ('U',"H1'"),('U',"H2'"),('U',"H3'"),('U',"H4'"),('U',"H5'"),('U',"H5''"),('U',"HO2'"),('U',"HO3'"),('U',"HO5'")]

if len(analysis_parameterlist) > 0:
	try:
		os.remove(str(outdir)+'simulation_results.dat')
	except:
		print
	filename_analyses=file(str(outdir)+'/simulation_results.dat','w')
	filename_analyses.write('# SIMULATION RESULTS, BASED ON INPUT PARAMETERS SPECIFIED IN FILE: '+str(sys.argv[1]))

shutil.copyfile(inputfile,outdir+'/inputfile')
# Iteration over all the different analyses
for i in range(len(analysis_parameterlist)):
	
	filename_analyses.write('\n\n# Results of Analysis number '+str(i+1)+'\n\n')

	########################
	# Case 1: This is an analysis at which timepoint a certain magnetization value is reached, on average
	if analysis_parameterlist[i][0]=='t':
		timepoint_list=[]
		cutofflevel=float(analysis_parameterlist[i][2])
		filename_analyses.write('# Determining the time points where the magnetization reaches a level of '+str(cutofflevel)+'\n\n')
		print analysis_parameterlist[i][1]
		# Check all atomic nuclei
		for k in range(len(atomnamelist)):
			counter=0
			# Check if the current nucleus is of interest (ie if it is in the specified list of nuclei to be analysed)
			if (  ( (restypelist[k],atomnamelist[k]) in analysis_parameterlist[i][1]) or ((str(sequencelist[k]),atomnamelist[k]) in analysis_parameterlist[i][1]) or ( ('imino' in analysis_parameterlist[i][1]) and ( (restypelist[k],atomnamelist[k]) in imino_protons  )) or ( ('aromatic' in analysis_parameterlist[i][1]) and ( (restypelist[k],atomnamelist[k]) in aromatic_protons  )) or ( ('amino' in analysis_parameterlist[i][1]) and ( (restypelist[k],atomnamelist[k]) in amino_protons  )) or ( ('sugar' in analysis_parameterlist[i][1]) and ( (restypelist[k],atomnamelist[k]) in sugar_protons  ))  ):
				for t in range(len(timesteplist)):
				
					if ((averaged_vector_of_vector_of_magnetizations[t][k])+1)>cutofflevel:
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
