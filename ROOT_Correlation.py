#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for file organization for hypothesis test code in C++ (load ROOT 6 before running the hypothesis testing code)
# Run this code with ROOT 6 by going inside the file .bashrc and reloading using source ~/.profile
# This code prepares the hypothesis files to be tested in the other codes


# This code I created lets me plot "Correlation Plots" from the same ROOT file. 
# In this case the output file came from the analysis performed in MadAnalysis5.
# For a specific spin just change the input string to the desired case! CHANCHE THE SPECIFIC YLABEL TO ADJUST TO THE SPECIFIC CASE!

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

def extract_data(input_ROOT,ttree,varA,varB):
	
	file_data=uproot3.open(input_ROOT)
	x_data=file_data[ttree].array(varA)
	y_data=file_data[ttree].array(varB)
	
	return x_data, y_data 

#----------------------------------------------#
#------------ Plot DATA from ROOT -------------#
#----------------------------------------------#

def Running_pT(input_directory, varA, varB, ttree, mass_list, spin_case):
    hep.style.use('ATLAS')
    for i, mass in enumerate(mass_list):
    # Construct the file name
        direc = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
        x, y = extract_data(direc, ttree, varA, varB)
        # Plot the data on the current axis
        plt.hist2d(x, y, bins=[150, 150], range=[[0, 300], [0, 300]],cmap='twilight',density=True)
        cbar = plt.colorbar()
        cbar.set_label("$\\frac{1}{N} \\frac{d^2 N}{(d p_T)^2}$")
        # Set the labels for the current axis
        plt.xlabel("Parton Level")
        plt.ylabel("Rec. with Truth-Match")
        #plt.ylabel("Rec. without Truth-Match")
    
        hep.atlas.text(dictionary_spin[spin_case] + str(mass) +" GeV",loc=2)
        hep.atlas.text(dictionary_pT[varB],loc=4)
        plt.text(0.25, 0.82, "MadGraph5_aMC@NLO", transform=plt.gca().transAxes, ha="center", va="center")
        plt.tight_layout()
    
        plt.savefig(str(varA) + "_" +str(spin_case) + "_" + str(mass) + "GeV" + ".png",dpi=300)
        plt.show()
    return True

def Running_eta(input_directory, varA, varB, ttree, mass_list, spin_case):
    hep.style.use('ATLAS')
    for i, mass in enumerate(mass_list):
    # Construct the file name
        direc = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
        x, y = extract_data(direc, ttree, varA, varB)
        # Plot the data on the current axis
        plt.hist2d(x, y, bins=[150, 150], range=[[-5, 5], [-5, 5]],cmap='twilight',density=True)
        cbar = plt.colorbar()
        cbar.set_label("$\\frac{1}{N} \\frac{d^2 N}{(d \eta)^2}$")
        # Set the labels for the current axis
        plt.xlabel("Parton Level")
        plt.ylabel("Rec. with Truth-Match")
        #plt.ylabel("Rec. without Truth-Match")
    
        hep.atlas.text(dictionary_spin[spin_case] + str(mass) +" GeV",loc=2)
        hep.atlas.text(dictionary_eta[varB],loc=4)
        plt.text(0.25, 0.82, "MadGraph5_aMC@NLO", transform=plt.gca().transAxes, ha="center", va="center")
        plt.tight_layout()
    
        plt.savefig(str(varA) + "_" +str(spin_case) + "_" + str(mass) + "GeV" + ".png",dpi=300)
        plt.show()
    return True


#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#

dictionary_pT={
"PtTTbar":r"$p_{T_{(t \overline{t})}}$ [GeV]",
"PtTopQ":r"$p_{T_{(t)}}$ [GeV]",
"PtWp":r"$p_{T_{(W^+)}}$ [GeV]",
"PtNeu":r"$p_{T_{(\nu)}}$ [GeV]",
"PtHiggs":r"$p_{T_{(Y_0)}}$ [GeV]",
}

dictionary_eta={
"EtaTTbar":r"$\eta_{(t \overline{t})}$ [GeV]",
"EtaTopQ":r"$\eta_{(t)}$ [GeV]",
"EtaWp":r"$\eta_{(W^+)}$ [GeV]",
"EtaNeu":r"$\eta_{(\nu)}$ [GeV]",
"EtaHiggs":r"$\eta_{(Y_0)}$ [GeV]",
}

dictionary_spin = {
"spin_0_plus": r"$m_{Y_{0^+}}=$",
"spin_0_minus": r"$m_{Y_{0^-}}=$",
"spin_1_plus": r"$m_{Y_{1^+}}=$",
"spin_1_minus": r"$m_{Y_{1^-}}=$",
"spin_2_plus": r"$m_{Y_{2^+}}=$",
}

input_1=""

#------------     REC/PARTON     -------------#

#Running_pT(input_1,"RecPtTTbar","PtTTbar","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_pT(input_1,"RecPtTTbar","PtTTbar","WithCuts;1",[1,10,100,125], "spin_0_minus")
#Running_eta(input_1,"RecEtaTTbar","EtaTTbar","WithCuts;1",[1,10,100,125], "spin_0_plus")

# TOP QUARK
#Running_pT(input_1,"RecPtTopQ","PtTopQ","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_eta(input_1,"RecEtaTopQ","EtaTopQ","WithCuts;1",[1,10,100,125], "spin_0_plus")

# W+
#Running_pT(input_1,"RecPtWp","PtWp","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_eta(input_1,"RecEtaWp","EtaWp","WithCuts;1",[1,10,100,125], "spin_0_plus")

# NEUTRINO
#Running_pT(input_1,"RecPtNeu","PtNeu","WithCuts;1",[1,10,100,125], "spin_0_plus")


#------------     EXP/PARTON     -------------#

#Running_pT(input_1,"ExpPtTTbar","PtTTbar","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_pT(input_1,"ExpPtTTbar","PtTTbar","WithCuts;1",[1,10,100,125], "spin_0_minus")
#Running_eta(input_1,"ExpEtaTTbar","EtaTTbar","WithCuts;1",[1,10,100,125], "spin_0_plus")

# TOP QUARK
#Running_pT(input_1,"ExpPtTopQ","PtTopQ","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_eta(input_1,"ExpEtaTopQ","EtaTopQ","WithCuts;1",[1,10,100,125], "spin_0_plus")

# W+
#Running_pT(input_1,"ExpPtWp","PtWp","WithCuts;1",[1,10,100,125], "spin_0_plus")
#Running_eta(input_1,"ExpEtaWp","EtaWp","WithCuts;1",[1,10,100,125], "spin_0_plus")

# NEUTRINO
#Running_pT(input_1,"ExpPtNeu","PtNeu","WithCuts;1",[1,10,100,125], "spin_0_plus")
