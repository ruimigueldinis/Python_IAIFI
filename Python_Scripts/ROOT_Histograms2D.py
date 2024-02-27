#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for file organization for hypothesis test code in C++ (load ROOT 6 before running the hypothesis testing code)
# Run this code with ROOT 6 by going inside the file .bashrc and reloading using source ~/.profile
# This code prepares the hypothesis files to be tested in the other codes


# This code I created lets me plot two observables from the same ROOT file. 
# In this case the output file came from the analysis performed in MadAnalysis5.
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

def extract_data(input_ROOT,ttree,varA,varB):
	
	file_data=uproot3.open(input_ROOT)
	x_data=file_data[ttree].array(varA)
	y_data=file_data[ttree].array(varB)
	
	return x_data, y_data 

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

def Running_Mass(input_directory, varA, varB, ttree, mass_list, spin_case, output_name):
    # Set the size of the figure
    a4=(11.69,8.27)
    num_masses = len(mass_list)
    # Calculate the number of rows needed for the subplots
    num_rows = num_masses // 3 + (num_masses % 3 > 0)
    # Create the subplots
    fig, axs = plt.subplots(num_rows, 3, figsize=a4)
    # Flatten the axes array to make it easier to loop over
    axs = axs.ravel()
    # Loop over the masses
    for i, mass in enumerate(mass_list):
        # Construct the file name
        direc = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
        x, y = extract_data(direc, ttree, varA, varB)
        # Plot the data on the current axis
        axs[i].hist2d(x, y, bins=[100, 100], cmap=plt.cm.Reds)
        # Set the labels for the current axis
        axs[i].set_xlabel(dictionary_inX[varA])
        axs[i].set_ylabel(dictionary_inY[varB])
        # Set the title for the current axis
        axs[i].set_title("Spin " + spin_title_dict[spin_case] + str(mass) + " GeV")

    # Remove any unused axes
    for j in range(i+1, num_rows*3):
        fig.delaxes(axs[j])

    # Adjust the layout of the figure
    fig.tight_layout()
    # Save plot as pdf and show plot
    plt.savefig(output_name + ".pdf")
    plt.show()
    return True

def Running_Spin(input_directory, varA, varB, ttree, mass, spin, output_name):
    # Set the size of the figure
    a4=(11.69,8.27)
    num_spin = len(spin)
    # Calculate the number of rows needed for the subplots
    num_rows = num_spin // 3 + (num_spin % 3 > 0)
    # Create the subplots
    fig, axs = plt.subplots(num_rows, 3, figsize=a4)
    # Flatten the axes array to make it easier to loop over
    axs = axs.ravel()
    for i, spin_case in enumerate(spin):
        # Create input file name for current spin case and mass
        input_file = input_directory + "SET01_dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic.root"
        # Extract data for current spin case and mass
        x, y = extract_data(input_file, ttree, varA, varB) 
        # Create 2D histogram for current spin case and mass
        axs[i].hist2d(x, y, bins=[100, 100], cmap=plt.cm.Reds)
        # Set labels and title for current plot
        axs[i].set_xlabel(dictionary_inX[varA])
        axs[i].set_ylabel(dictionary_inY[varB])
        axs[i].set_title("Spin " + spin_title_dict[spin_case] + str(mass) + " GeV")
        
    # Remove any unused axes
    for j in range(i+1, num_rows*3):
        fig.delaxes(axs[j])

    # Adjust layout of subplots
    fig.tight_layout()
    # Save plot as pdf and show plot
    plt.savefig(output_name+".pdf")
    plt.show()
    return True

#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#

activated = set_style('Latin Modern Roman',10)

input_1=""

# Define dictionary for input variable names
dictionary_inX = {
    "AngHiggsTopQ":r"$\Delta \phi_{Y, t}$",
    "AngHiggsTbarQ":r"$\Delta \phi_{Y, \bar{t}}$"
}

dictionary_inY = {
    "AngHiggsTopQ":r"$\Delta \phi_{Y, t}$",
    "AngHiggsTbarQ":r"$\Delta \phi_{Y, \bar{t}}$"
}

input_files=""

Running_Mass(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", [0,1,10,100,125,1000], "spin_0_plus", "Scalar_Mediator_Angle_HiggsTopQ_vs_HiggsTbarQ")
Running_Mass(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", [0,1,10,100,125,1000], "spin_0_minus", "Pseudoscalar_Mediator_Angle_HiggsTopQ_vs_HiggsTbarQ")
Running_Mass(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", [0,1,10,100,1000], "spin_1_minus", "Axial_Mediator_Angle_HiggsTopQ_vs_HiggsTbarQ")
Running_Mass(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", [0,1,10,100,1000], "spin_1_plus", "Vector_Mediator_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Mass(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", [0,1,10,100,1000], "spin_2_plus", "Tensor_Mediator_Angle_HiggsTopQ_vs_HiggsTbarQ")

Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 0, ["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"], "0GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 1, ["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"], "1GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 10, ["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"], "10GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 100, ["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"], "100GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 125, ["spin_0_plus", "spin_0_minus"], "125GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
#Running_Spin(input_1, "AngHiggsTopQ", "AngHiggsTbarQ", "NoCuts;1", 1000, ["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"], "1000GeV_Angle_HiggsTopQ_vs_HiggsTbarQ")
