#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for file organization for hypothesis test code in C++ (load ROOT 6 before running the hypothesis testing code)
# Run this code with ROOT 6 by going inside the file .bashrc and reloading using source ~/.profile
# This code prepares the hypothesis files to be tested in the other codes


# This code plots 1D histograms from ROOT files. In this case the Output file came from the combinatorial analysis performed in MadAnalysis5
# It's important to know that you can run Running_Mass by itself, but it will only generate a single output pdf. For a greater number of variables,
# is better to run_Running_Mass_chunks, as it divides the input variables in chunks of 6 and saves the plots in separate files.
# For a small number of variables, I advice to comment the line where it deletes axes, and adjust manually the number of rows in Running_Mass
# This is because, I defined a subplot as being in the format of an A4. When a single plot is made maybe is best to adjust the figsize.
# Command density is for normalizing the plots , and alpha parameter is for the oppacity of the colors.

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

def Running_Mass(input_directory, ttree, mass_list, spin_case, output_name, variables, num_bins,case_signal_back):
    # Set the size of the figure
    a4=(8.27,11.69)
    num_variables = len(variables)
    # Calculate the number of rows needed for the subplots
    num_rows = num_variables // 2 + (num_variables % 2 > 0)
    # Create the subplots
    fig, axs = plt.subplots(num_rows,2, figsize=a4, sharey=False)
    fig.suptitle("Spin " + spin_title_dict[spin_case])
    # Flatten the axes array to make it easier to loop over
    axs = axs.ravel()
    # Loop over the masses
    for i, mass in enumerate(mass_list):
        # Construct the file name
        direc = input_directory + "dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic/Output/Combinatorial_ttH_"+str(case_signal_back)+"_m"+str(mass)+"gev_s"+spin_case_dict[spin_case]+".root"
        k=0
        # Loop over the variables
        for j, var in enumerate(variables):
            data = extract_data(direc, ttree, var)
            # Plot the data on the current axis
            axs[k].hist(data, bins=num_bins, histtype=u'step',label=str(mass) + "GeV",density=True)
            # Set the labels for the current axis
            axs[k].set_xlabel(dictionary_inX[var])
            axs[k].set_ylabel(dictionary_inY[var])
            axs[k].legend(loc='best', fontsize=5,frameon=False)
            k+=1
    
    # Remove any unused axes
    if (num_variables % 2) != 0:
    	fig.delaxes(axs[num_variables])

    # Adjust the layout of the figure
    fig.tight_layout()
    # Save plot as pdf and show plot
    plt.savefig(output_name + ".pdf")
    plt.show()
    return True

def Running_Spin(input_directory, ttree, mass, spin, output_name, variables, num_bins,case_signal_back):
    # Set the size of the figure
    a4=(8.27,11.69)
    num_variables = len(variables)
    # Calculate the number of rows needed for the subplots
    num_rows = num_variables // 2 + (num_variables % 2 > 0)
    # Create the subplots
    fig, axs = plt.subplots(num_rows,2, figsize=a4, sharey=False)
    fig.suptitle("$m_Y = $ " + str(mass) + " GeV")
    # Flatten the axes array to make it easier to loop over
    axs = axs.ravel()
    # Loop over the masses
    for i, spin_case in enumerate(spin):
        # Construct the file name
        direc = input_directory + "dm_simp_mass" + str(mass) + "gev_" + spin_case + "_dileptonic/Output/Combinatorial_ttH_"+str(case_signal_back)+"_m"+str(mass)+"gev_s"+ spin_case_dict[spin_case] +".root"
        k=0
        # Loop over the variables
        for j, var in enumerate(variables):
            data = extract_data(direc, ttree, var)
            # Plot the data on the current axis
            axs[k].hist(data, bins=num_bins,histtype=u'step',label="Spin " + spin_title_dict[spin_case],density=True)
            # Set the labels for the current axis
            axs[k].set_xlabel(dictionary_inX[var])
            axs[k].set_ylabel(dictionary_inY[var])
            axs[k].legend(loc='best', fontsize=8,frameon=False)
            k+=1
    
    # Remove any unused axes
    if (num_variables % 2) != 0:
        fig.delaxes(axs[num_variables])

    # Adjust the layout of the figure
    fig.tight_layout()
    # Save plot as pdf and show plot
    plt.savefig(output_name + ".pdf")
    plt.show()
    return True
	
#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#
	
activated = set_style('Latin Modern Roman',10)

input_1=""

variables_1=["CombDeltaR_lp_bT","CombDeltaPhi_lp_bT","CombDeltaThe_lp_bT","CombMass_lp_bT","CombDeltaR_ln_bbarTbar","CombDeltaPhi_ln_bbarTbar","CombDeltaThe_ln_bbarTbar","CombDeltaR_lp_bbarTbar","CombDeltaPhi_lp_bbarTbar","CombDeltaThe_lp_bbarTbar","CombDeltaR_ln_bT","CombDeltaPhi_ln_bT","CombDeltaThe_ln_bT"]

# Define dictionary for input variable names
dictionary_inX={
"CombDeltaR_lp_bT":r"$\Delta R_{\ell^{+}, b_{T}}$",
"CombDeltaPhi_lp_bT":r"$\Delta \Phi_{\ell^{+}, b_{T}}$",
"CombDeltaThe_lp_bT":r"$\Delta \Theta_{\ell^{+}, b_{T}}$",  
"CombMass_lp_bT":r"$M_{\ell^{+}, b_{T}}$",  

"CombDeltaR_ln_bbarTbar":r"$\Delta R_{\ell^{-}, \bar{b}_{\bar{T}}}$",   
"CombDeltaPhi_ln_bbarTbar":r"$\Delta \Phi_{\ell^{-}, \bar{b}_{\bar{T}}}$",      
"CombDeltaThe_ln_bbarTbar":r"$\Delta \Theta_{\ell^{-}, \bar{b}_{\bar{T}}}$",    
    
"CombDeltaR_lp_bbarTbar":r"$\Delta R_{\ell^{+}, \bar{b}_{\bar{T}}}$",       
"CombDeltaPhi_lp_bbarTbar":r"$\Delta \Phi_{\ell^{+}, \bar{b}_{\bar{T}}}$",          
"CombDeltaThe_lp_bbarTbar":r"$\Delta \Theta_{\ell^{+}, \bar{b}_{\bar{T}}}$",

"CombDeltaR_ln_bT":r"$\Delta R_{\ell^{-}, b_{T}}$",
"CombDeltaPhi_ln_bT":r"$\Delta \Phi_{\ell^{-}, b_{T}}$",    
"CombDeltaThe_ln_bT":r"$\Delta \Theta_{\ell^{-}, b_{T}}$",  
}

dictionary_inY={
"CombDeltaR_lp_bT":r"$\frac{1}{N}\frac{dN}{d\Delta R_{\ell^{+}, b_{T}}}$",
"CombDeltaPhi_lp_bT":r"$\frac{1}{N}\frac{dN}{d  \Delta \Phi _{\ell^{+}, b_{T}}  }$",
"CombDeltaThe_lp_bT":r"$\frac{1}{N}\frac{dN}{d  \Delta \Theta_{\ell^{+}, b_{T}}   }$",  
"CombMass_lp_bT":r"$\frac{1}{N}\frac{dN}{d  M_{\ell^{+}, b_{T}}    }$", 
    
"CombDeltaR_ln_bbarTbar":r"$\frac{1}{N}\frac{dN}{d   \Delta R_{\ell^{-}, \bar{b}_{\bar{T}}}  }$",       
"CombDeltaPhi_ln_bbarTbar":r"$\frac{1}{N}\frac{dN}{d \Delta \Phi_{\ell^{-}, \bar{b}_{\bar{T}}}    }$",      
"CombDeltaThe_ln_bbarTbar":r"$\frac{1}{N}\frac{dN}{d \Delta \Theta_{\ell^{-}, \bar{b}_{\bar{T}}}    }$",
        
"CombDeltaR_lp_bbarTbar":r"$\frac{1}{N}\frac{dN}{d  \Delta R_{\ell^{+}, \bar{b}_{\bar{T}}}   }$",       
"CombDeltaPhi_lp_bbarTbar":r"$\frac{1}{N}\frac{dN}{d  \Delta \Phi_{\ell^{+}, \bar{b}_{\bar{T}}}   }$",      
"CombDeltaThe_lp_bbarTbar":r"$\frac{1}{N}\frac{dN}{d \Delta \Theta_{\ell^{+}, \bar{b}_{\bar{T}}}    }$",
        
"CombDeltaR_ln_bT":r"$\frac{1}{N}\frac{dN}{d  \Delta R_{\ell^{-}, b_{T}}   }$",     
"CombDeltaPhi_ln_bT":r"$\frac{1}{N}\frac{dN}{d   \Delta \Phi_{\ell^{-}, b_{T}}  }$",        
"CombDeltaThe_ln_bT":r"$\frac{1}{N}\frac{dN}{d \Delta \Theta_{\ell^{-}, b_{T}}   }$",       
}

#Running_Mass(input_1,"WithCutsComb",[0,1,10,100,125,1000],"spin_0_plus","Scalar_signal",variables_1,100,"signal")
#Running_Mass(input_1,"WithCutsComb",[0,1,10,100,125,1000],"spin_0_plus","Scalar_back_allComb",variables_1,100,"back_allComb")
#Running_Mass(input_1,"WithCutsComb",[0,1,10,100,125,1000],"spin_0_plus","Scalar_back_noOver",variables_1,100,"back_noOver")

Running_Spin(input_1,"WithCutsComb",0,["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"],"0GeV_signal",variables_1,100,"signal")
#Running_Spin(input_1,"WithCutsComb",1,["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"],"1GeV_signal",variables_1,100,"signal")
#Running_Spin(input_1,"WithCutsComb",10,["spin_0_plus", "spin_0_minus", "spin_1_plus", "spin_1_minus"],"10GeV_signal",variables_1,100,"signal")
