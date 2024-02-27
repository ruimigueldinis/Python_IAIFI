#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for file organization for hypothesis test code in C++ (load ROOT 6 before running the hypothesis testing code)
# Run this code with ROOT 6 by going inside the file .bashrc and reloading using source ~/.profile
# This code prepares the hypothesis files to be tested in the other codes


# This code performs a binary classification task on a dataset using the K-nearest neighbors (KNN) algorithm.
# In this case the output file came from the analysis performed in MadAnalysis5.


from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler 
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report, confusion_matrix

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
#------------ Read Data from ROOT -------------#
#----------------------------------------------#

def extract_data(input_ROOT, ttree, var):
    # Open the ROOT file and get the TTree
    file_data = uproot3.open(input_ROOT)
    # Extract the variable data
    data = file_data[ttree].array(var)

    return data

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

#----------------------------------------------#
#-----------  Important Functions  ------------#
#----------------------------------------------#

def QuadVector_Extractor(input_directory,input_ROOT,var):
	# It gets this values from the variables Pt,Eta,Phi and M 
	Exp_Pt, Exp_Eta, Exp_Phi, Exp_M, Pt, Eta, Phi, M 	= get_PtEtaPhiM(input_directory+input_ROOT,var[0],var[1],var[2],var[3])
	QuadVector_p, QuadVector_p_exp 						= set_PtEtaPhiM(Exp_Pt, Exp_Eta, Exp_Phi, Exp_M, Pt, Eta, Phi, M)

	# Transpose the arrays to have the same number of features. Necessary for the concatenation of the datasets
	QuadVector_p 										= np.transpose(QuadVector_p)
	QuadVector_p_exp 									= np.transpose(QuadVector_p_exp)
	
	return QuadVector_p, QuadVector_p_exp

#----------------------------------------------#
#------------     Directories     -------------#
#----------------------------------------------#

place=""

Top_Quark       = ["PtTopQ","EtaTopQ","PhiTopQ","mTopQ"]
AntiTop_Quark   = ["PtTbarQ","EtaTbarQ","PhiTbarQ","mTbarQ"]

Positive_Lepton = ["PtLepP","EtaLepP","PhiLepP","mLepP"]
Negative_Lepton = ["PtLepN","EtaLepN","PhiLepN","mLepN"]

B_Top_Quark_Sys = ["PtBt","EtaBt","PhiBt","mBt"]
AntiB_AntiTop_Quark_Sys = ["PtBbtb","EtaBbtb","PhiBbtb","mBbtb"]


QuadVector_ttbar, QuadVector_exp_ttbarvjets			= QuadVector_Extractor(place, "SET02_ttbar_dileptonic_NLO.root", Top_Quark)	
QuadVector_ttbarvjets, QuadVector_exp_ttbarvjets	= QuadVector_Extractor(place, "SET02_ttbar_v_jets_NLO_dilep.root", Top_Quark)	
QuadVector_Higgs, QuadVector_exp_Higgs				= QuadVector_Extractor(place, "SET02_h_NLO.root", Top_Quark)	
QuadVector_zvjets, QuadVector_exp_zvjets			= QuadVector_Extractor(place, "SET02_z_v_jets_NLO.root", Top_Quark)	

#print(QuadVector_ttbar.shape)

#----------------------------------------------#
#------------   Machine Learning    -----------#
#----------------------------------------------#


# Concatenate the four datasets
x = np.concatenate((QuadVector_ttbar, QuadVector_ttbarvjets, QuadVector_Higgs, QuadVector_zvjets), axis=0)

# Create labels for the datasets (0 for ttbar, 1 for ttbarvjets, 2 for Higgs, and 3 for zvjets)
y_ttbar = np.zeros(QuadVector_ttbar.shape[0])
y_ttbarvjets = np.ones(QuadVector_ttbarvjets.shape[0])
y_Higgs = np.full(QuadVector_Higgs.shape[0], 2)
y_zvjets = np.full(QuadVector_zvjets.shape[0], 3)
y = np.concatenate((y_ttbar, y_ttbarvjets, y_Higgs, y_zvjets), axis=0)

# Split the data into training and test sets
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3)

# Fit the KNN model
knn = KNeighborsClassifier(n_neighbors=5)
knn.fit(x_train, y_train)

# Predict the labels for the test set
y_pred = knn.predict(x_test)

# Calculate the confusion matrix
cm = confusion_matrix(y_test, y_pred)
#print(cm)
#print(classification_report(y_test, y_pred))

hep.style.use('ATLAS')
cmap=plt.cm.Greens
classes=[r"Dileptonic $t\bar{t}$", r"$t\bar{t}+$jets", "Higgs", r"Z+$\nu$"]

plt.imshow(cm, interpolation='nearest', cmap=cmap)
plt.title("KNeighborsClassifier on Dominant Backgrounds",fontsize=12)

plt.colorbar()
tick_marks = np.arange(len(classes))
plt.xticks(tick_marks, classes)
plt.yticks(tick_marks, classes)

fmt = 'd'
thresh = cm.max() / 2.
for i, j in np.ndindex(cm.shape):
	plt.text(j, i, format(cm[i, j], fmt), horizontalalignment="center",color="white" if cm[i, j] > thresh else "black")
                 
plt.savefig("Confusion_Matrix_KNN_Background.pdf")                
plt.show()
