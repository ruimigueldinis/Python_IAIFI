#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for file organization for hypothesis test code in C++ (load ROOT 6 before running the hypothesis testing code)
# Run this code with ROOT 6 by going inside the file .bashrc and reloading using source ~/.profile
# This code prepares the hypothesis files to be tested in the other codes


# This code I created lets me construct observables and also plot 1D histograms of the same observable from ROOT files. 
# In this case the output file came from the analysis performed in MadAnalysis5.
# I must create new functions or variables to check! Just create a new function and adapt the Running_Spin_Observable 
# function! It's different from the other codes I created, because in this case it plots all the spins in one figure!
# For a specific spin just change the input string to the desired case! I'm trying to make all my python codes in the same format to be better to adapt it accordingly.

import uproot3
import ROOT
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
import math

from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator

from ROOT import Math
from ROOT import TLorentzVector

#----------------------------------------------#
#--------------   Formatting    ---------------#
#----------------------------------------------#

def set_style(font,size):
    rc('text', usetex=True)
    rc('font',**{'family':'serif','serif':[font]})
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.size'] = size
    return True

#----------------------------------------------#
#------------ Read Data from ROOT -------------#
#----------------------------------------------#

def extract_data(input_ROOT, ttree, var):
    # Open the ROOT file and get the TTree
    file_data = uproot3.open(input_ROOT)
    # Extract the variable data
    data = file_data[ttree].array(var)

    return data

spin_title_dict = {
    "spin_0_plus": r"$0^{+}$"+ " ",
    "spin_0_minus": r"$0^{-}$"+ " ",
    "spin_1_plus": r"$1^{+}$"+ " ",
    "spin_1_minus": r"$1^{-}$"+ " ",
    "spin_2_plus": r"$2^{+}$"+ " ",
}

spin_case_dict = {
    "spin_0_plus": "0_plus",
    "spin_0_minus":"0_minus",
    "spin_1_plus": "1_plus",
    "spin_1_minus": "1_minus",
    "spin_2_plus": "2_plus",
}

#----------------------------------------------#
#------------ Plot DATA from ROOT -------------#
#----------------------------------------------#

def get_PtEtaPhiM(input_ROOT,var_Pt,var_Eta,var_Phi,var_M):	
	# Save the Exp variables into an array (including the cases without -10000)
	Exp_Pt  = extract_data(input_ROOT,'WithCuts;1','Exp'+var_Pt)
	Exp_Eta = extract_data(input_ROOT,'WithCuts;1','Exp'+var_Eta)	
	Exp_Phi = extract_data(input_ROOT,'WithCuts;1','Exp'+var_Phi)
	Exp_M   = extract_data(input_ROOT,'WithCuts;1','Exp'+var_M)
	# Save the Parton-level variables into an array
	Pt      = extract_data(input_ROOT,'WithCuts;1',var_Pt)
	Eta     = extract_data(input_ROOT,'WithCuts;1',var_Eta)	
	Phi     = extract_data(input_ROOT,'WithCuts;1',var_Phi)
	M       = extract_data(input_ROOT,'WithCuts;1',var_M)

	# Filter the Exp and Parton variables into an array without '-1000' default values
	# that the analysis creates on the .root file
	Exp_Pt  = Exp_Pt[Exp_Pt > -1000]
	Exp_Eta = Exp_Eta[Exp_Eta > -1000]
	Exp_Phi = Exp_Phi[Exp_Phi > -1000]
	Exp_M   = Exp_M[Exp_M > -1000]

	Pt      = Pt[Pt > -1000]
	Eta     = Eta[Eta > -1000]
	Phi     = Phi[Phi > -1000]
	M       = M[M > -1000]

	return Exp_Pt, Exp_Eta, Exp_Phi, Exp_M, Pt, Eta, Phi, M
	

def set_PtEtaPhiM(Exp_Pt, Exp_Eta, Exp_Phi, Exp_M, Pt, Eta, Phi, M):
    # Create empty lists to hold the resulting four momentum
    Exp_Momentum = []
    Momentum     = []

    # Loop over all the values
    for i in range(len(Exp_Pt)):
        
        # Create a TLorentzVector with the Exp Pt, Eta, Phi, and M
        Exp    = TLorentzVector()
        Exp.SetPtEtaPhiM(Exp_Pt[i], Exp_Eta[i], Exp_Phi[i], Exp_M[i])
        
        # Create a TLorentzVector with the Parton-level Pt, Eta, Phi, and M
        Parton = TLorentzVector()
        Parton.SetPtEtaPhiM(Pt[i], Eta[i], Phi[i], M[i])

        # Append the particle's Px, Py, Pz, and Energy to each list
        Exp_Momentum.append([Exp.Px(), Exp.Py(), Exp.Pz(), Exp.E()])
        Momentum.append([Parton.Px(), Parton.Py(), Parton.Pz(), Parton.E()])

    # Transpose the resulting arrays: The way numpy goes through lists
    Exp_Momentum = np.transpose(np.array(Exp_Momentum))
    Momentum     = np.transpose(np.array(Momentum))

    return Exp_Momentum, Momentum

# MAKE NOTICE THAT get_P..., the output is not a numpy array. I can create a new function if it is important in the future.

#----------------------------------------------#
#-----------  Important Functions  ------------#
#----------------------------------------------#

def triple_product(vec_A,vec_B,vec_C):
	# This function calculates the triple product: (A x B) . C = A . (B x C)
	# First: B x C = (b2*c3 - b3*c2)i + (b3*c1 - b1*c3)j + (b1*c2 - b2*c1)k
	i = vec_B[:, 1] * vec_C[:, 2] - vec_B[:, 2] * vec_C[:, 1]
	j = vec_B[:, 2] * vec_C[:, 0] - vec_B[:, 0] * vec_C[:, 2]
	k = vec_B[:, 0] * vec_C[:, 1] - vec_B[:, 1] * vec_C[:, 0]

	norm = np.sqrt(vec_A[:, 0]**2 + vec_A[:, 1]**2 + vec_A[:, 2]**2 ) * np.sqrt(i**2 + j**2 + k**2)
	cross_product = vec_A[:, 0] * i + vec_A[:, 1] * j + vec_A[:, 2] * k

	return cross_product, norm

#----------------------------------------------#
#-----------  Kinematic Variables  ------------#
#----------------------------------------------#

def observable_b4_like(input_ROOT,var_A,var_B):
	
	Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A = get_PtEtaPhiM(input_ROOT,var_A[0],var_A[1],var_A[2],var_A[3])
	Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B = get_PtEtaPhiM(input_ROOT,var_B[0],var_B[1],var_B[2],var_B[3])
	
	four_p, four_p_exp = set_PtEtaPhiM(Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A)
	four_n, four_n_exp = set_PtEtaPhiM(Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B)
	
	var_positive=four_p_exp[2]                                                   # Gets pz from TopQ
	var_negative=four_n_exp[2]                                                   # Gets pz from TbarQ 
	mod_positive=np.sqrt(four_p_exp[0]**2 + four_p_exp[1]**2 + four_p_exp[2]**2) # Calculates the |p|^2 of the TopQ
	mod_negative=np.sqrt(four_n_exp[0]**2 + four_n_exp[1]**2 + four_n_exp[2]**2) # Calculates the |p|^2 of the TbarQ
	
	b4_like= (var_positive)*(var_negative)/mod_positive/mod_negative             # (pz_t * pz_tbar)/|p_t|^2 / |p_tbar|^2 
	return b4_like


def Running_Spin_Observable(input_directory,mass,spin,output_name,num_bins,var_positive,var_negative):
	hep.style.use('ATLAS')
	for i, spin_case in enumerate(spin):
		# Create input file name for current spin case and mass
		input_file = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
		sb4 = observable_b4_like(input_file,var_positive,var_negative)
		# Create 1D histogram for current spin case and mass
		plt.hist(sb4, bins=num_bins, histtype=u'step', label="Spin " + spin_title_dict[spin_case],density=True)
		# Set labels and title for current plot        
		plt.xlabel("$b_4$")
		plt.ylabel("$\\frac{1}{N}\\frac{dN}{d b_4} $")
		hep.atlas.label("$b_4$ with Top Quarks",loc=2,data=True,lumi=100,com=14,lumi_format="{0:.1f}")
		plt.legend(loc='best',frameon=False)

 	# Save plot as pdf and show plot
	plt.savefig(output_name+".pdf")
	plt.show()
	return True


#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#

place="/Users/ruimiguelsilva/documents/madanalysis5/bin/ttH_dilep_REC5_September/Output/Output_ttY_Spin0Plus_0GeV/"

Top_Quark       = ["PtTopQ","EtaTopQ","PhiTopQ","mTopQ"]
AntiTop_Quark   = ["PtTbarQ","EtaTbarQ","PhiTbarQ","mTbarQ"]

Positive_Lepton = ["PtLepP","EtaLepP","PhiLepP","mLepP"]
Negative_Lepton = ["PtLepN","EtaLepN","PhiLepN","mLepN"]

B_Top_Quark_Sys = ["PtBt","EtaBt","PhiBt","mBt"]
AntiB_AntiTop_Quark_Sys = ["PtBbtb","EtaBbtb","PhiBbtb","mBbtb"]

Running_Spin_Observable(place,0,["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"],"b4_with_tops_scalar_0gev",100,Top_Quark,AntiTop_Quark)
Running_Spin_Observable(place,1,["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"],"b4_with_tops_scalar_1gev",100,Top_Quark,AntiTop_Quark)