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
        Exp_Momentum.append([Exp.Px(), Exp.Py(), Exp.Pz()])
        Momentum.append([Parton.Px(), Parton.Py(), Parton.Pz()])

    # Transpose the resulting arrays: The way numpy goes through lists
    Exp_Momentum = np.transpose(np.array(Exp_Momentum))
    Momentum     = np.transpose(np.array(Momentum))

    return Exp_Momentum, Momentum

# MAKE NOTICE THAT REMOVED THE ENERGY FROM THE 4-VECTORS ------->>> Exp_Momentum.append([Exp.Px(), Exp.Py(), Exp.Pz(), Exp.E()])


'''
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

#----------------------------------------------#
#-----------  Kinematic Variables  ------------#
#----------------------------------------------#


def observable_N9(input_ROOT, var_A, var_B, var_C, var_D):

    Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A = get_PtEtaPhiM(input_ROOT,var_A[0],var_A[1],var_A[2],var_A[3]) #TOP
    Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B = get_PtEtaPhiM(input_ROOT, var_B[0], var_B[1], var_B[2], var_B[3])  # l+
    Exp_Pt_C, Exp_Eta_C, Exp_Phi_C, Exp_M_C, Pt_C, Eta_C, Phi_C, M_C = get_PtEtaPhiM(input_ROOT, var_C[0], var_C[1], var_C[2], var_C[3])  # b quark/l-
    Exp_Pt_D, Exp_Eta_D, Exp_Phi_D, Exp_M_D, Pt_D, Eta_D, Phi_D, M_D = get_PtEtaPhiM(input_ROOT, var_D[0], var_D[1], var_D[2], var_D[3])  # Tbar

    four_topQ, four_topQ_exp = set_PtEtaPhiM(Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A) #T(Anti T)-QUARK
    four_Lep, four_Lep_exp = set_PtEtaPhiM(Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B)  # ± LEPTON
    four_Bt, four_Bt_exp = set_PtEtaPhiM(Exp_Pt_C, Exp_Eta_C, Exp_Phi_C, Exp_M_C, Pt_C, Eta_C, Phi_C, M_C)  # B(Anti B)-QUARK / lep negativo
    four_tBar, four_tBar_exp = set_PtEtaPhiM(Exp_Pt_D, Exp_Eta_D, Exp_Phi_D, Exp_M_D, Pt_D, Eta_D, Phi_D, M_D) #Tbar

    # Dot Product of (p_t . p_l+) = (E_t * E_l+) - (px_t*px_l+ + py_t*py_l+ + py_t*py_l+)
    product_positive = np.sum(four_topQ_exp[3,:] * four_Lep_exp[3,:]) - np.sum(four_topQ_exp[:3,:] * four_Lep_exp[:3,:], axis=0)
    product_negative = np.sum(four_tBar_exp[3,:] * four_Bt_exp[3,:]) - np.sum(four_tBar_exp[:3,:] * four_Bt_exp[:3,:], axis=0)

    # Removing the Energy component
    PxPyPz_topQ	= four_topQ_exp[:3, :]
    PxPyPz_LepP	= four_Lep_exp[:3, :]
    PxPyPz_Bt	= four_Bt_exp[:3, :]
    PxPyPz_tBar	= four_tBar_exp[:3, :]

    nt= (-PxPyPz_topQ/172.5) + (172.5/(product_positive))*PxPyPz_LepP
    ntbar= (PxPyPz_tBar/172.5) - (172.5/(product_negative))*PxPyPz_Bt

    mod_nt=np.sqrt(nt[0]**2 + nt[1]**2 + nt[2]**2)  # Calculates the |p|^2 of the TopQ
    mod_ntbar = np.sqrt(ntbar[0]**2 + ntbar[1]**2 + ntbar[2]**2)  # Calculates the |p|^2 of the Lepton +

    print(mod_ntbar)

    ntx=nt[2]
    ntbarx=ntbar[2]

    print(np.shape(ntx))


    #kZ = np.array([1/(np.sqrt(2)),1/(np.sqrt(2)), 0])
    k = np.array([0,0,1])

    #dot_product = np.dot(nt.T,k)
    #dot_product2 = np.dot(ntbar.T,k)
    #print(np.shape(dot_product))
    # Calculate the cross product manually
    #nt_Z=np.cross(nt.T,kZ).T
    #ntbar_Z=np.cross(ntbar.T,kZ).T
    #print("nt_Z Shape:", nt_Z.shape)
    #print("nt_Z:", nt_Z)
    #my_obs= np.sign(nt_Z*ntbar_Z)*np.arcsin(nt_Z*ntbar_Z)

    #my_obs= nt*ntbar
    my_obs= (ntx)*(ntbarx)/mod_nt/mod_ntbar

    print(np.shape(my_obs))
    return my_obs

'''

def observable_b2(input_ROOT,var_A,var_B):
    
    Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A = get_PtEtaPhiM(input_ROOT,var_A[0],var_A[1],var_A[2],var_A[3])
    Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B = get_PtEtaPhiM(input_ROOT,var_B[0],var_B[1],var_B[2],var_B[3])
    
    four_topQ, four_topQ_exp = set_PtEtaPhiM(Exp_Pt_A, Exp_Eta_A, Exp_Phi_A, Exp_M_A, Pt_A, Eta_A, Phi_A, M_A) #TOP QUARK
    four_Tbar, four_Tbar_exp = set_PtEtaPhiM(Exp_Pt_B, Exp_Eta_B, Exp_Phi_B, Exp_M_B, Pt_B, Eta_B, Phi_B, M_B) #ANTI TOP
    
    mod_topQ=np.sqrt(four_topQ[0]**2 + four_topQ[1]**2 + four_topQ[2]**2)   # Calculates the |p|^2 of the TopQ
    mod_Tbar=np.sqrt(four_Tbar[0]**2 + four_Tbar[1]**2 + four_Tbar[2]**2)   # Calculates the |p|^2 of the TBar

    normaTopQ=four_topQ/mod_topQ
    normaTbar=four_Tbar/mod_Tbar

    # Define the z-axis unit vector

    kZ = np.array([0,1/(np.sqrt(2)),1/(np.sqrt(2))])
    #kZ = np.array([0,1,0])

    # Calculate the cross product manually
    TopQ_Z=np.cross(normaTopQ.T,kZ).T
    TBar_Z=np.cross(normaTbar.T,kZ).T

    b2= np.sum(TopQ_Z[:3,:] * TBar_Z[:3,:], axis=0)
    return b2

def Running_Spin_Observable(input_directory,mass,spin,output_name,num_bins,var_A,var_B):
    hep.style.use('ATLAS')
    for i, spin_case in enumerate(spin):
        # Create input file name for current spin case and mass
        input_file = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
        b2 = observable_b2(input_file,var_A,var_B)
        # Create 1D histogram for current spin case and mass
        plt.hist(b2, bins=num_bins, histtype=u'stepfilled',alpha=0.5,edgecolor='black', linewidth=1, label="Spin " + spin_title_dict[spin_case],density=True)
        # Set labels and title for current plot        
        plt.xlabel("$\\widetilde{b_{2}}$")
        plt.ylabel("$\\frac{1}{N}\\frac{dN}{d \\widetilde{b_{2}}} $")

        #plt.xlabel("$b_{2}$")
        #plt.ylabel("$\\frac{1}{N}\\frac{dN}{d b_{2}} $")
        hep.atlas.text("$m_{Y_0}=$" + str(mass) +" GeV",loc=0)
        hep.atlas.text("$\\widetilde{b_{2}} = (\\overrightarrow{p_t} \\times \hat{k_d}) \cdot (\\overrightarrow{p_{\overline{t}}} \\times \hat{k_d}); \\qquad \hat{k_d}=(0,\\frac{1}{\sqrt{2}},\\frac{1}{\sqrt{2}})$",loc=3)
        #hep.atlas.text("$\\widetilde{b_{2}} = (\\overrightarrow{p_t} \\times \hat{k_d}) \cdot (\\overrightarrow{p_{\overline{t}}} \\times \hat{k_d}); \\qquad \hat{k_d}=(0,1,0)$",loc=3)
        #hep.atlas.text("$b_{2} = (\\overrightarrow{p_t} \\times \hat{k_d}) \cdot (\\overrightarrow{p_{\overline{t}}} \\times \hat{k_d}); \\qquad \hat{k_d}=(0,0,1)$",loc=3)
        plt.legend(loc='best',frameon=False)

    # Save plot as pdf and show plot
    plt.savefig(output_name + str(mass) + "GeV" + ".png",dpi=300)
    plt.show()

    return True
#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#

input_directory="/Users/ruimiguelsilva/documents/madanalysis5/bin/ttH_dilep_REC5_September/Output/"

Top_Quark       = ["PtTopQ","EtaTopQ","PhiTopQ","mTopQ"]
AntiTop_Quark   = ["PtTbarQ","EtaTbarQ","PhiTbarQ","mTbarQ"]

Top_Quark_ttH       = ["PtTopQ_ttH","EtaTopQ_ttH","PhiTopQ_ttH","mTopQ_ttH"]
AntiTop_Quark_ttH   = ["PtTbarQ_ttH","EtaTbarQ_ttH","PhiTbarQ_ttH","mTbarQ_ttH"]

Positive_Lepton = ["PtLepP","EtaLepP","PhiLepP","mLepP"]
Negative_Lepton = ["PtLepN","EtaLepN","PhiLepN","mLepN"]

B_Top_Quark_Sys = ["PtBt","EtaBt","PhiBt","mBt"]
AntiB_AntiTop_Quark_Sys = ["PtBbtb","EtaBbtb","PhiBbtb","mBbtb"]


Running_Spin_Observable(input_directory,125,["spin_0_plus", "spin_0_minus"],"B2_Gunion_Parton_sqrt",100,Top_Quark,AntiTop_Quark)

