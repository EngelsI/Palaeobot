#!/usr/bin/env python3
"""
Created on Fri Oct  7 16:00:55 2022

@author: Ian
"""
##################### PACKAGES #######################################
import os
import math
import pandas as pd
from zipfile import ZipFile
from pyteomics import mgf
import matplotlib.pyplot as plt
import numpy as np
from ms2pip.ms2pipC import MS2PIP
from deeplc import DeepLC
from ms2pip.single_prediction import SinglePrediction
import csv
import warnings
from Bio import SeqIO
import re
import time
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings('ignore')
#####################################################################

################### Functions #######################################
def insilico_peptide(pep, new_pep, temp, MZ):
    new_pep_list = {}
    for AA in pep:
        for AA2 in new_pep:
            if AA2[-1].isdigit()==True and AA[-1].isdigit()==True:
                continue
            if (MZ-0.01< (AA_code[AA]+new_pep[AA2]) < MZ+0.01):
                add = AA2+'|'+AA
                temp.add(add)
            elif (AA_code[AA]+new_pep[AA2]) < MZ:
                new_pep_list[AA2+'|'+AA] = AA_code[AA]+new_pep[AA2]
            elif (AA_code[AA]+new_pep[AA2]) > MZ:
                temp.add(AA2)
    return new_pep_list, temp

def jaccard(list1, list2, pep):
    mass = find_mass(pep)
    if list(list1).count(0)>1 or len(list1)==0:
        return 0
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    j = ((float(intersection) / union))
    return j*mass/len(pep)

def insilico_spectrum(peptides, spectrums, intensity, ion, to_add): 
    if len(spectrums)==0:
        return [], []
    if ('k' in peptides.split('|') or 'r' in peptides.split('|') or 'Carbamidomethyl10' in peptides.split('|')) and ('K' in peptides.split('|') or 'R' in peptides.split('|')  or 'Carbamidomethyl2' in peptides.split('|')):
        final_y_array = [0]*len(spectrums)
        return spectrums, np.array(final_y_array)
    if ion == 'Y' and ('k' in peptides[1:] or 'r' in peptides[1:]):
        final_y_array = [0]*len(spectrums)
        return spectrums, np.array(final_y_array)
    if ion == 'Y' and 'Gln->pyro-Glu' in peptides or 'Glu->pyro-Glu' in peptides:
        return [], []
    if ion == 'B' and ('k' in peptides.split('|') or 'r' in peptides.split('|')):
        return [], []
    x_array = []
    x_point = ion_types[ion]+to_add #19 for Y ions, 1 for b ions
    for AA in peptides.split('|'):
        x_array.append(AA_code[AA]+x_point)
        x_point += AA_code[AA] 
  
    x_array = sorted(x_array)
    if min(spectrums)>max(x_array):
        final_y_array = [0]*len(spectrums)
        return spectrums, np.array(final_y_array)
    pep_spectrum = {}
    for loc,num in enumerate(spectrums):
        former=0
        for value in x_array:
            if value-dev <= num <= value+dev:
                if intensity[loc]>former:
                    pep_spectrum[value]=num

    # if 3>len(peptides.split('|'))>6:
    #     double_missing_peak = False
    #     zero_test = ion_types[ion]
    #     zero_count = 0
    #     first = 0
    #     much_miss = False
    #     for nummer, amino in enumerate(peptides.split('|')):
    #         zero_test += AA_code[amino]
    #         if zero_test not in pep_spectrum.keys() and double_missing_peak == True:
    #             final_y_array = [0]*len(spectrums)
    #             return spectrums, np.array(final_y_array)
    #         elif zero_test not in pep_spectrum.keys() and nummer <= 1:
    #             first += 1
    #         elif zero_test not in pep_spectrum.keys():
    #             double_missing_peak = True
    #             if much_miss == False:
    #                 much_miss = True
    #                 double_missing_peak = False
    #             # if first == 1 and nummer ==2:
    #             #     final_y_array = [0]*len(spectrums)
    #             #     return spectrums, np.array(final_y_array)
    #             zero_count += 1
    #         else:
    #             double_missing_peak = False
    #             much_miss = False
            
    #     if zero_count > (abs(math.ceil(len(peptides.split('|'))/3))):
    #         final_y_array = [0]*len(spectrums)
    #         return spectrums, np.array(final_y_array)
    final_y_array = []
    for y_loc, point in enumerate(spectrums):
        if point in pep_spectrum.values():
            final_y_array.append(intensity[y_loc])
        else:
            final_y_array.append(1e-15)
    return spectrums, np.array(final_y_array)

def pearson_correlation(real_spectra, insilico_spectra, overlap = False, i=0):
    # if overlap == 'Overlap' and list(insilico_spectra).count(0)>2:
    #     corr_matrix = np.corrcoef(real_spectra, insilico_spectra)
    #     return corr_matrix
    if len(real_spectra)<=2:
        return[[0,0]]
    if list(real_spectra) == list(insilico_spectra):
        return [[1,1]]
    if list(insilico_spectra).count(0)>1: #or sum(insilico_spectra)/len(insilico_spectra) <= min(insilico_spectra)
        return [[0,0]]
    #corr_matrix = [[0,np.sum(insilico_spectra)/np.sum(real_spectra)]]
    corr_matrix = np.corrcoef(real_spectra, insilico_spectra)
    if math.isnan(corr_matrix[0][1]) == True:
        return [[0,0]]
    return corr_matrix

def find_treshold(pearson_list,MZ=250):
    Y = [num[-1] for num in pearson_list if num[-1]>0]
    T = Y
    Y = sorted(Y)[::-1]
    if len(Y)<=500:
        print('no threshold')
        return 1e-15
    much = len(Y)*0.1
    Y = Y[int(len(Y)*0.1):]
    
    z=np.polyfit([x for x, a in enumerate(Y)], Y, 1)
    
    if much + len([num for num in Y if num >=z[1]]) <min(MZ*3,1500):
        z = [0,z[1]*0.85]
    while len([num for num in T if num >= z[1]])>min(MZ*(MZ/420)*3,2500):
        z = [0, z[1]*1.01]
    print('threshold is', z[1], 'with maximum ', min(MZ*(MZ/420)*3,2500))
    return z[1]
  
def find_treshold_overlap(pearson_list):
    #cluster with pearson jaccard score to find best group. Only keep points from top group
    points = [[num[1], num[2]*num[3]] for num in pearson_list if num[3]>0.45 and num[2]>0.45]
    if len(points)==0:
        return []
    a = [num[0] for num in points]
    pearson_list = [num for num in pearson_list if num[3]>0.45 and num[2]>0.45 and num[1]>=max(a)*0.4]
    points = [[num[0],num[1]] for num in points if num[0]>=max(a)*0.4]
    if len(points)==0:
        return []
    
    a = [num[0] for num in points]
    b = [num[1] for num in points]
    kmeans_pr = [1 if (num[1]>=max(b)*0.9 and num[0]>=max(a)*0.7) or (num[0]>=max(a)*0.6 and num[1]>=max(b)*0.8) or num[1]>=0.9 or num[0]>=0.85*max(a) else 0 for num in points]
    add = 0
    while kmeans_pr.count(1)>250 and add<0.1:
        to_much_results = int(kmeans_pr.count(1))
        print(to_much_results, ' is to much, reducing! ...')
        add += 0.01
        kmeans_pr = [1 if (num[1]>=max(b)*(0.9+add) and num[0]>=max(a)*(0.75+add)) or (num[0]>=max(a)*(0.75+add) and num[1]>=max(b)*(0.8+add)) or num[1]>=(0.98+add) or num[0]>=(max(a)*(0.9+add)) else 0 for num in points]
    plt.scatter(a,b, c=kmeans_pr)
    plt.title('Selection of potential peptides \n match score to spectrum')
    plt.ylabel('pearson correlation * coverage')
    plt.xlabel('Peak counting score')
    plt.show()
    
    pearson_list = [pearson_list[i] for i, u in enumerate(kmeans_pr) if u == 1]
    pearson_list = [pearson_list[i] for i, u in enumerate(pearson_list) if u[3]>0.45 and u[2]>0.45 and u[1]>=max(a)*0.4]
    return pearson_list

def make_plots(xpoints, ypoints, peptide, PC_list, ion):
    cut = find_treshold(PC_list)
    ylist = []
    xlist = []
    for a, AA in enumerate(PC_list):
        ylist.append(AA[1])
        xlist.append(a)
    plt.title('Pearson correlation plot')
    plt.bar(xlist, ylist)
    plt.axhline(cut, c='r', alpha=0.5)
    plt.show()

def find_mass(peptide):
    mass = 0
    for AA in peptide:
        mass += AA_code[AA]
    return mass            

def select_top_scores(PC_list, treshold):
    x = {}
    for item in PC_list:
        if item[-1] >= treshold and item[-1] != 0:
            to_add = item[0].split('|')
            if int(item[1]) != 0:
                x[str(item[1])+'*'+'|'.join(to_add)] = find_mass(to_add)+item[1]
            else:
                x['|'.join(to_add)] = find_mass(to_add)
    # plt.bar([len(PC_list),len(x)], ['possible fingerprints', 'remaining fingerprints'])
    # plt.title('Remaining potential TAGs \n after peak counting')
    # plt.show()
    return x

def find_overlap(y_ions, b_ions, yextra, bextra): #adapt for downstream
    y_ions = y_ions.split('|')[::-1]
    b_ions = b_ions.split('|')
    to = 1+int(len(y_ions)/3)
    to_b = 1+int(len(b_ions)/3)
    if spectrum['params']['charge'][0]==1:
        print('overlap of 1+')
        to = len(y_ions)
        to_b = 0
    if (''.join(b_ions[-1]) not in ''.join(y_ions)[:to]) or (''.join(y_ions[0]) not in ''.join(b_ions)[-to_b:]):
        return[]
    if 'Methyl' in y_ions or'Gln->pyro-Glu' in y_ions or 'Glu->pyro-Glu' in y_ions:
        return []
    if 'Methyl' in b_ions[1:] or 'Gln->pyro-Glu' in b_ions[1:] or 'Glu->pyro-Glu' in b_ions[1:] or 'k' in b_ions or 'r' in b_ions or 'Carbamidomethyl10' in b_ions:
        return []
    for i in range(to_b, len(b_ions)):
        b_ions_overlap = b_ions[i:len(b_ions)]
        if b_ions_overlap ==y_ions[0:len(b_ions_overlap)]:
            massa = find_mass(b_ions[0:i]+y_ions)+yextra+bextra
            if pp_mass-1-math.floor(21/spectrum['params']['charge'][0]) <= massa/spectrum['params']['charge'][0] <= pp_mass+1-math.floor(21/spectrum['params']['charge'][0]): #include h3o+ and h+
                adding = [b_ions[0:i] + y_ions, yextra, bextra]
                result = []
                ptm = []
                if bextra != 0:
                    to_ptm = '0|'+N_terms[bextra]
                    ptm.append(to_ptm)
                for loc, element in enumerate(adding[0]):
                    if element in PTM:
                        result.append(PTM[element])
                        to_ptm = str(loc+1)+'|'+''.join(c for c in element if c.isdigit()==False)
                        ptm.append(to_ptm)
                    else:
                        result.append(element)
                if ('k' in result or 'r' in result) and ('K' in result or 'R' in result):
                    return []
                return ("".join(result), "|".join(ptm), adding)
    return [] 

def make_overlap_spectrum(y_peptide, b_peptide, spectrum, intensity, yextra, bextra):
    original_peak = []
    x_array_Y = []
    x_point_Y = ion_types['Y'] + yextra
    H2O = False
    NH3 = False
    heavy_kr = False
    for AA in y_peptide:
        x_array_Y.append(AA_code[AA]+x_point_Y)
        original_peak.append(AA_code[AA]+x_point_Y)
        if AA in mods['-H2O'] or H2O == True:
            H2O = True
            x_array_Y.append(AA_code[AA]+x_point_Y-18.010565)
        if AA in mods['-NH3'] or NH3 == True:
            NH3 = True
            if heavy == True and AA in ['k', 'r']:
                x_array_Y.append(AA_code[AA]+x_point_Y-18.026549)
                heavy_kr = True
            elif heavy != True:
                x_array_Y.append(AA_code[AA]+x_point_Y-17.026549)
            else:
                if heavy_kr == True:
                    x_array_Y.append(AA_code[AA]+x_point_Y-18.026549)
                x_array_Y.append(AA_code[AA]+x_point_Y-17.026549)  
        x_point_Y += AA_code[AA]    
    x_array_B = []
    x_point_B = ion_types['B'] + bextra
    H2O = False
    NH3 = False
    heavy_kr = False
    for AA in b_peptide:
        x_array_B.append(AA_code[AA]+x_point_B)
        original_peak.append(AA_code[AA]+x_point_B)
        if AA in mods['-H2O'] or H2O == True:
            H2O = True
            x_array_B.append(AA_code[AA]+x_point_B-18.010565)
        if AA in mods['-NH3'] or NH3 == True:
            NH3 = True
            if heavy == True and AA.upper() in ['k', 'r']:
                x_array_B.append(AA_code[AA]+x_point_B-18.026549)
                heavy_kr=True
            elif heavy != True:
                x_array_B.append(AA_code[AA]+x_point_B-17.026549)
            else:
                if heavy_kr==True:
                    x_array_B.append(AA_code[AA]+x_point_B-18.026549)
                x_array_B.append(AA_code[AA]+x_point_B-17.026549)
        x_point_B += AA_code[AA] 
    x_array = set(x_array_B + x_array_Y) #takes unique points, otherwise maybe bias
    x_array = sorted(x_array)#[num for num in x_array if max(spectrum)+19 >= num >= min(spectrum)])
    original_peak = sorted(original_peak)#[num for num in original_peak if max(spectrum)+19 >= num >= min(spectrum)])
    to_zero = 0
    pep_spectrum ={}
    for num in spectrum:
        for value in x_array:
            if value-dev <= num <= value+dev:
                pep_spectrum[value]=num
                if value in original_peak:
                    to_zero += 1
    # if to_zero<len(original_peak)/3:
    #     final_y_array = [0]*len(spectrum)
    #     return spectrum, np.array(final_y_array)
    final_y_array = []
    for y_loc, point in enumerate(spectrum):
        if point in pep_spectrum.values():
            final_y_array.append(intensity[y_loc])
        else:
            final_y_array.append(1e-15)
    return spectrum, np.array(final_y_array)

def make_final_plots(xpoints, ypoints, y_pep, b_pep, PC_list):
    d = ypoints
    a, c = make_overlap_spectrum(y_pep, b_pep, xpoints, d, 0, 0)
    #fig1 = plt.figure()
    #ax1 = fig1.add_subplot(111)
    # fig1, ax1 = plt.subplots()
    #ax2 = ax1.twinx()
    plt.stem(xpoints, d, 'b', label='True signal', markerfmt=' ')
    plt.stem(a, -c, 'g', label='in-silico', markerfmt=' ')
    plt.title('Cleaned spectrum VS theoretical spectrum')
    plt.ylabel('Intensity')
    plt.xlabel('m/z')
    #fig1.tight_layout()
    plt.show()
        
    ylist = []
    xlist = []
    for a, AA in enumerate(PC_list):
        ylist.append(AA[1])
        xlist.append(a)
    
    plt.plot(xlist, ylist)
    plt.show()

def tic_normalize(intensity):
        
        intensity = np.array(intensity)
        return intensity / np.sum(intensity)

def find_intensity(item, spectrum, intensity, ms2pip, ptm): #if 0 use ms2pip intensity?
    AA_code = extra_AA   
    true_intensity_B = []
    true_intensity_Y = []
    loc_ptm = []
    which = []
    bextra = 0
    do = False
    for el in ptm.split('|'):
        try:
            if int(el)==0:
                do = True
            else:
                loc_ptm.append(int(el))
        except ValueError:
            if do == True:
                do = False
                for ex, name in N_terms.items():
                    if name == el:
                        bextra = ex
            else:
                which.append(el)
            pass
    item = item.replace('I', 'L')
    item = list(item)
    for i, number in enumerate(loc_ptm):
        for k,v in PTM.items():
            if v==item[number-1] and which[i] == ''.join(c for c in k if c.isdigit()==False):
                item[number-1]= k
                break
    y_peptide = item[::-1][:-1]
    b_peptide = item[:-1]

    x_array_Y = []
    x_point_Y = ion_types['Y'] 
    for AA in y_peptide:
        x_array_Y.append(AA_code[AA]+x_point_Y)
        x_point_Y += AA_code[AA]
    x_array_B = []
    x_point_B = ion_types['B'] +bextra
    for AA in b_peptide:
        x_array_B.append(AA_code[AA]+x_point_B)
        x_point_B += AA_code[AA]
     
    pep_spectrum = {}
    for loc, value in enumerate(x_array_Y):
        done = False
        for num in spectrum:
            if value-dev <= num <= value+dev:
                pep_spectrum[value]=num
                done = True
        if done == False:
            if loc<= 2 or loc>= len(y_peptide)-1:#value < min(spectrum):
                pep_spectrum[value] = 'NP'+str(loc)
            else:
                pep_spectrum[value] = 0
    final_y_array = []
    for point in sorted(pep_spectrum.keys()):
        if pep_spectrum[point] ==0:
            final_y_array.append(1e-15)
            true_intensity_Y.append(True)
        elif 'NP' in str(pep_spectrum[point]):
            true_intensity_Y.append(False)
            #final_y_array.append(ms2pip[int(pep_spectrum[point][2:])+len(y_peptide)])
        else:
            itemindex = np.where(spectrum == pep_spectrum[point])[0][0]
            final_y_array.append(intensity[itemindex])
            true_intensity_Y.append(True)
    pep_spectrum = {}
    for loc, value in enumerate(x_array_B):
        done = False
        for num in spectrum:
            if value-dev <= num <= value+dev:
                pep_spectrum[value]=num
                done = True
        if done == False:
            if loc <= 2 or loc>= len(b_peptide)-1:#value < min(spectrum):
                pep_spectrum[value] = 'NP'+str(loc)
            else:
                pep_spectrum[value] = 0
    final_b_array = []
    for point in sorted(pep_spectrum.keys()):
        if pep_spectrum[point] ==0:
            final_b_array.append(1e-15)
            true_intensity_B.append(True)
        elif 'NP' in str(pep_spectrum[point]):
            true_intensity_B.append(False)
            #final_b_array.append(ms2pip[int(pep_spectrum[point][2:])])
        else:
            itemindex = np.where(spectrum == pep_spectrum[point])[0][0]
            final_b_array.append(intensity[itemindex])
            true_intensity_B.append(True)
    final_int = true_intensity_B + true_intensity_Y
    ms2pip = np.array(ms2pip)
    ms2pip = ms2pip[final_int]
    return np.array(final_b_array + final_y_array), ms2pip #tic_normalize(np.array(final_b_array + final_y_array))

def catch_asterix(temp):
    new = []
    for i in temp:
        if '*' in i and len(i)>0:
            element = i.split('*')
            if element[1][0]=='|':
                new.append((element[1][1:], float(element[0])))
            else: 
                new.append((element[1], float(element[0])))
        else:
            new.append((i,0))
    return new
        
def make_list_pep(MZ, ion, AAs):
    temp =set()  
    pep = list(AAs.keys())
    if ion=='B':
        new_pep = {}
        for i,u in AAs.items():
            new_pep[i]=u
        for item in N_terms.keys():
            new_pep[str(item)+'*']=item
    elif ion=='Y':
        new_pep = {}
        for i,u in AAs.items():
            new_pep[i]=u
        for item in C_terms.keys():
            new_pep[str(item)+'*']=item
    print('making the in-silico peptides')
    while len(new_pep)!=0:
        new_pep, temp = insilico_peptide(pep, new_pep, temp, MZ)
    temp = catch_asterix(temp)
    return temp

def find_pep_ions_part1(xpoints_data, ypoints_data, MZ, ion, temp):
    # MZ_range = xpoints_data#spectrum['m/z array']
    # ind = np.argsort(MZ_range)
    # MZ_range = MZ_range[ind]
    print('Start Pearson Corr')
    print(len(temp))
    xpoints = xpoints_data[xpoints_data<=(MZ+19)]
    ypoints = ypoints_data[0:len(xpoints)]
    PC_list =[]
    for peptide, extra in temp:
        a, c = insilico_spectrum(peptide, xpoints, ypoints, ion, extra)
        pc = jaccard(c,ypoints,peptide.split('|'))
        # found_intensity = np.sum(c)
        # all_intensity = np.sum(ypoints_data)
        # pc = pc#*found_intensity/all_intensity
        #pc = pearson_correlation(ypoints, c, i=i)[0][1]
        add = tuple([peptide, extra, pc])
        PC_list.append(add)
       
    PC_list = sorted(PC_list, key=lambda x: x[-1])[::-1]
    
    #make_plots(xpoints, ypoints, PC_list[0][0], PC_list, ion)
    return PC_list

def find_pep_ions_iterations(PC_list, MZ, ion):
    temp = set()
    pep = list(AA_code.keys())
    new_pep = select_top_scores(PC_list, find_treshold(PC_list,MZ))
    while len(new_pep)!=0:
        new_pep, temp = insilico_peptide(pep, new_pep, temp, MZ)      
    print('Start Pearson Corr')
    temp = catch_asterix(temp)
    print(len(temp))
    if len(temp)>=35e4:
        return []
    ind = xpoints_data<=(MZ+19)
    xpoints = xpoints_data[ind]
    ypoints = ypoints_data[ind]
    
    PC_list =[]
    for peptide, extra in temp:
        a, c = insilico_spectrum(peptide, xpoints, ypoints, ion, extra)
        pc = jaccard(c,ypoints, peptide.split('|'))
        # found_intensity = np.sum(c)
        # all_intensity = np.sum(ypoints_data)
        # pc = pc#*found_intensity/all_intensity
        #pc = pearson_correlation(ypoints, c, i=i)[0][1]
        add = tuple([peptide,extra, pc])
        PC_list.append(add)
       
    PC_list = sorted(PC_list, key=lambda x: x[1])[::-1]
    #make_plots(xpoints, ypoints, PC_list[0][0], PC_list, ion)
    return PC_list

def find_overlap_sequence(Y_list, B_list):
    print('start making overlap')
    overlap_output = []
    for y_ion in Y_list:
        for b_ion in B_list:
            overlap = find_overlap(y_ion[0], b_ion[0], y_ion[1], b_ion[1]) 
            if len(overlap)>0:
                overlap_output.append(overlap)
    overlap_output_final = []
    for element in overlap_output:
        if element not in overlap_output_final:
            overlap_output_final.append(element)
    # overlap_output = sorted(overlap_output, key=lambda x: x[2])[::-1]
    print('found ', len(overlap_output)-len(overlap_output_final), ' double sequences')
    print('amount of possible sequences: ',len(overlap_output_final))
    plt.bar(['Y-ion fragments', 'B-ion fragments', 'Resulting peptides'],[len(Y_list), len(B_list),len(overlap_output_final)])
    plt.title('Amount of overlapping ion fragment to peptide')
    plt.show()
    return overlap_output_final

def score_overlap(overlap_output, x, y):
    
    print('final pearson correlation')
    
    PC_list =[]
    for item in overlap_output:
        a, c = make_overlap_spectrum(item[2][0][::-1], item[2][0], x, y, item[2][1], item[2][2])
        j = jaccard(y, c, item[2][0])
        pc = pearson_correlation(y, c)[0][1]
        found_intensity = np.sum(c)
        all_intensity = np.sum(y)
        add = tuple([item, j, found_intensity/all_intensity, pc]) #j*found_intensity/all_intensity
        PC_list.append(add)
        # if 'GAYVEVTAk' in item[0]:
        #     print(make_overlap_spectrum(item[2][::-1], item[2], xpoints_data, ypoints_data))
        #     print(jaccard(ypoints_data, c, item[0]))
        #     print(pearson_correlation(ypoints_data, c)[0][1])
        #     print(np.sum(c))
    
    PC_list = sorted(PC_list, key=lambda x: x[3])[::-1]

    if len(PC_list)==0:
        return []
    
    begin = len(PC_list)
    PC_list = find_treshold_overlap(PC_list)
    plt.bar(['Before', 'After'],[begin, len(PC_list)])
    plt.title('remaining peptides after peptide to spectrum match scoring')
    plt.show()
    print('kept ', len(PC_list), 'of ', begin)
    
    combined_peptide_list = set()
    for item, j_score, cov, pc in PC_list:
        combined_peptide_list = combined_peptide_list| {(item[0], item[1], j_score, cov, pc)}
    return combined_peptide_list

def score_pip(combined_peptide_list, xpoints_data, ypoints_data):
    pip_score = []
    for item, ptm, pc_score, cov in combined_peptide_list:
        ms2pip_sp = SinglePrediction()
        mz, intensity, annotation = ms2pip_sp.predict(item, '-', spectrum['params']['charge'][0])
        y = item[::-1]
        b = item
        to_test = find_intensity(y[:-1], b[:-1],  xpoints_data, ypoints_data, intensity)
        if len(to_test)==len(intensity): #If to_test shorter than intensity, this means that the total found peptide cannot be explained by the spectrum as it is larger than measured. 
            pc = pearson_correlation(to_test, intensity, 'Overlap')[0][1]
            pip_score.append([item, ptm, pc, pc_score, cov])
    
    pip_score = sorted(pip_score, key=lambda x: x[2])[::-1]
    #plot_pip_score(pip_score)  
    return pip_score
    
def plot_pip_score(pip_score):
    ylist = []
    xlist = []
    for a, AA in enumerate(pip_score):
        ylist.append(AA[1])
        xlist.append(a)
    
    plt.plot(xlist, ylist)
    plt.show()
    return 

def deeplc(title, contamination_file):
    print('Start DeepLC')
    calibration_file = "./"+contamination_file
    peptide_file = "./records_"+title+".csv"
    pep_df = pd.read_csv(peptide_file, sep=",")
    cal_df = pd.read_csv(calibration_file, sep=",")
    cal_df['modifications'] = cal_df['modifications'].fillna("")
    pep_df['modifications'] = pep_df['modifications'].fillna("")
    dlc = DeepLC(verbose=True)
    dlc.calibrate_preds(seq_df=cal_df)
    preds = dlc.make_preds(seq_df=pep_df)
    print('End  of deeplc')
    return preds

def crap_f():
    print('making database')
    crap = {}
    files_own = []
    crap_concat = {}
    
    for fastafile in os.walk('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/palaeo_db'):#fasta_db
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta'):
                files_own.append('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/palaeo_db/'+i) #fasta_db
    for i in files_own:
        concat = ''
        print(i.split('/')[-1])
        if 'other' not in i:
            continue
        for record in SeqIO.parse(i, "fasta"):
            if 'X' in record.seq or 'B' in record.seq or 'U' in record.seq or 'J' in record.seq or 'O' in record.seq or 'Z' in record.seq:
                continue
            crap[record.seq]=record.description
            concat+=str(record.seq)+'X'
        crap_concat[concat]=i.split('/')[-1]
    concat=''
    for record in SeqIO.parse("C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/crap.fasta.txt", "fasta"):
        crap[record.seq]=record.id
        concat+=str(record.seq)+'X'
    for record in SeqIO.parse('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/contaminants.fasta', "fasta"):
        crap[record.seq]=record.id
        concat+=str(record.seq)+'X'
    lysyl = 'MHKRTYLNACLVLALAAGASQALAAPGASEMAGDVAVLQASPASTGHARFANPNAAISAAGIHFAAPPARRVARAAPLAPKPGTPLQVGVGLKTATPEIDLTTLEWIDTPDGRHTARFPISAAGAASLRAAIRLETHSGSLPDDVLLHFAGAGKEIFEASGKDLSVNRPYWSPVIEGDTLTVELVLPANLQPGDLRLSVPQVSYFADSLYKAGYRDGFGASGSCEVDAVCATQSGTRAYDNATAAVAKMVFTSSADGGSYICTGTLLNNGNSPKRQLFWSAAHCIEDQATAATLQTIWFYNTTQCYGDASTINQSVTVLTGGANILHRDAKRDTLLLELKRTPPAGVFYQGWSATPIANGSLGHDIHHPRGDAKKYSQGNVSAVGVTYDGHTALTRVDWPSAVVEGGSSGSGLLTVAGDGSYQLRGGLYGGPSYCGAPTSQRNDYFSDFSGVYSQISRYFAP'
    crap[lysyl]='lysyl'
    concat+=lysyl+'X'
    keratin1 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSEGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1] = 'keratin1'
    concat+=keratin1+'X'
    keratin1_2 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSGGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1_2] = 'keratin1_2'
    concat+=keratin1_2+'X'
    crap_concat[concat]='contamination'
    concat=''
    for i in ['AETSELHTSLk', 'GAYVEVTAk', 'LGNEQGVSr', 'LVGTPAEEr', 'LDSTSLPVAk', 'AGLLVAEGVTk',
                  'LGLDFDSFr', 'GFTAYYLPr', 'SGGLLWQLVr', 'AVGANPEQLTr', 'SAEGLDASASLr',
                  'VFTPELVDVAk','VGNELQYVALr', 'YLELAPGVDNSk', 'DGTFAVDGPGVLAk', 'YDSLNNTEVSGLr',
                  'SPYVLTGPGVVEYk', 'ALENDLGVPSDATVk', 'AVYFYAPQLPLYANk', 'TVESLFPEEAETPGSAVr',  'STQAALDQLNGK', 'ALLVASGHLK']:
        concat += i.upper()+'X'
    crap_concat[concat]='pepcal'
    
    return crap_concat, crap

def personal_input(path):
    files_own = []
    for fastafile in os.walk(path+'/palaeo_db'):
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta'):
                files_own.append(path+'/palaeo_db/'+i)
    db = {}
    for i in files_own:
        for record in SeqIO.parse(i, "fasta"):
            db[record.seq]=record.description
    return db
            
def filters(title, df):
    already = []
    spectrum_calibration = {}
    print('save in csv file')
    with open('records_Calibration.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['seq', 'modifications', 'tr'])
        ms2pip_sp = SinglePrediction()
        df = df.sort_values(by=['rt_time'], ascending=True, ignore_index=True)
        for record, ptm, rt, x, y, cs, nr, fn, cov, pc_sc in df[['peptide', 'PTM', 'rt_time', 'xpoints_data', 'ypoints_data', 'charge_state', 'spectrum_nr', 'file_name', 'coverage', 'pc_score']].values:
            for sequence in crap.keys():
                if record.upper() in sequence and nr != 152:#◘hier wegdoen
                    mz, intensity, annotation = ms2pip_sp.predict(record.upper(), '-', cs)
                    intensity = [num if num >0 else 1e-15 for num in intensity]
                    to_test, new_intensity = find_intensity(record,  x, y, intensity, ptm)
                    if len(to_test)==len(new_intensity): 
                        pc = pearson_correlation(to_test, new_intensity)[0][1]
                        if (pc < 0.45 and (pc_sc<0.98 or len(record)<6 or pc <0)) or (len(record)<6 and pc_sc != 1):
                            print(crap[sequence], ' is not a good fit', pc)
                            break
                        print(crap[sequence])
                        if record.upper() not in already:
                            writer.writerow([record, ptm, rt])
                            already.append(record)
                            if nr not in spectrum_calibration.keys():
                                spectrum_calibration[nr]=[(record, ptm, crap[sequence], pc, fn, cov, nr)]
                            else:
                                spectrum_calibration[nr] = spectrum_calibration[nr]+[(record, ptm, crap[sequence], pc, fn, cov, nr)]
                        elif pc >=0.45:
                            if nr not in spectrum_calibration.keys():
                                spectrum_calibration[nr]=[(record, ptm, crap[sequence], pc, fn, cov, nr)]
                            else:
                                spectrum_calibration[nr] = spectrum_calibration[nr]+[(record, ptm, crap[sequence], pc, fn, cov, nr)]
                        break    
    if len(spectrum_calibration)>50:
        print('df length before:',len(df))
        df = df.loc[~df['spectrum_nr'].isin(list(spectrum_calibration.keys()))]
        to_csv = [[seq.upper(), ptm] for seq, ptm in df[['peptide', 'PTM']].values]
        peptides_to_test = write_csv(title, to_csv)
        print('kept', len(df), 'comming from spectra ', set(df['spectrum_nr'].values))
    
        preds=deeplc(title, 'records_Calibration.csv')
        df['deeplc_preds']=np.array(preds)
    
        difference = []
        for real, predicted in df[['rt_time', 'deeplc_preds']].values:
            difference.append(abs(real-predicted))
        df['difference']=np.array(difference)
        print('end of difference')
        
        X = np.array([[i, u] for i, u in df[['deeplc_preds', 'rt_time']].values])
        mini = min(df['rt_time'].values)*0.05
        Y2 = np.array([1 if i<=mini else 0 for i in df[['difference']].values])
    
        plt.title('plot of all predictions after DeepLC')
        plt.scatter([num[0] for i, num in enumerate(X)], [num[1] for i, num in enumerate(X)], c='b', alpha=0.01)
        plt.scatter([num[0] for i, num in enumerate(X) if Y2[i]==1], [num[1] for i, num in enumerate(X) if Y2[i]==1], c='g', alpha=0.1)
        plt.show()
        df['error_interval'] = np.array([mini]*len(df))
        df = df[abs(df['difference'])<=mini] 
        df = df.set_index('peptide', drop=False)
    else:
        print('not enough for deeplc')
        df = df.set_index('peptide', drop=False)
        df['difference']=np.array([0]*len(df))
    print('save results for ms2pip')
    nr = 0
    with open('peprec.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['spec_id','modifications','peptide','charge'])
        for seq, ch, ptm in df[['peptide', 'charge_state', 'PTM']].values:
            if ptm == '':
                ptm = '-'
            ch = int(''.join(c for c in str(ch) if c.isdigit()==True))
            writer.writerow([nr, '-', seq, ch])
            nr += 1
    print('start ms2pip')
    params = {
        "ms2pip": {
            "ptm": [
                "Oxidation,15.994915,opt,M",
                # "Phospho,79.966331,opt,S",
                # "Oxidation,15.994915,opt,P",
                # "Carbamyl,43.005814,opt,T",
                # "Carbamyl,43.005814,opt,R",
                # "kynurenin,3.994915,opt,W",
                # "dihydroxy,31.989829,opt,W",
                # "Gln->pyro-Glu,-18.010565,opt,Q",
                # "Phospho,79.966331,opt,Y",
                # "Phospho,79.966331,opt,T",
                # "dihydroxy,31.989829,opt,M",
                # "Carbonyl,13.979265,opt,K",
                # "Carbonyl,13.979265,opt,R",
                # "Deamidation,0.984016,opt,N",
                # "Deamidation,0.984016,opt,Q"
                ],
                "frag_method": "HCD",
                "frag_error": 0.02,
                "out": "csv",
                "sptm": [], "gptm": [],
                }}

    ms2pip = MS2PIP("peprec.csv",params = params, return_results=True)
    predictions = ms2pip.run()
    actual = []
    for pred in predictions['prediction']:
        actual.append((2**pred)-0.001)
    predictions['prediction']=np.array(actual)
    
    
    print('Start MS2PIP pearson')
    print(len(df))
    pip_score = []
    # ms2pip_sp = SinglePrediction()
    
    nr = 0
    for item, ptm, x, y, cs in df[['peptide', 'PTM', 'xpoints_data', 'ypoints_data', 'charge_state']].values:
        print(nr+1, ' of the', len(df))
        intensity = [num for num in predictions['prediction'][predictions['spec_id']==str(nr)].values]#mz, intensity, annotation = ms2pip_sp.predict(item.upper(), '-', cs)
        intensity = [num if num >0 else 1e-15 for num in intensity]
        to_test, new_intensity = find_intensity(item,  x, y, intensity, ptm)
        if len(to_test)==len(new_intensity): #If to_test shorter than intensity, this means that the total found peptide cannot be explained by the spectrum as it is larger than measured. 
            pc = pearson_correlation(to_test, new_intensity)[0][1]
            pip_score.append(pc)
        else:
            pip_score.append(0)
        nr += 1
    df['ms2pip_score'] = np.array(pip_score)
    plt.plot(sorted(pip_score)[::-1])
    plt.title('Peptide to MS2PIP scoring')
    plt.xlabel('Peptide number')
    plt.ylabel('MS2PIP vs peptide score')
    plt.show()
    return df, spectrum_calibration
         
def write_csv(title, combined_peptide_list):
    print('save in csv file')
    with open('records_'+title+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['seq', 'modifications'])
        for record, ptm in combined_peptide_list:
            writer.writerow([record, ptm])
    return 0
     
def set_data(xpoints_data, ypoints_data, pp_mass, charge, spectrum, tag=False):
          
    ind = np.argsort(xpoints_data)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    ypoints_data = tic_normalize(ypoints_data)
        
    if 'charge' in spectrum['params'].keys():
        if spectrum['params']['charge'][0] > 2:
            print(pp_mass)
            pp_mass = pp_mass-((spectrum['params']['charge'][0]-2)**2)/spectrum['params']['charge'][0]#((pp_mass*3)/2)-0.5
            print('checking a non-double charged!!', pp_mass)
    ind = []    
    for del_y, i in enumerate(xpoints_data):
        if spectrum['params']['pepmass'][0]-0.02 <i< spectrum['params']['pepmass'][0]+0.02:#precursor mass
            ind.append(False)
        elif spectrum['params']['pepmass'][0]-(18.010565/spectrum['params']['charge'][0])-0.02 <i< spectrum['params']['pepmass'][0]-(18.010565/spectrum['params']['charge'][0])+0.02: #precursor mass - H2O
            ind.append(False)
        elif spectrum['params']['pepmass'][0]-(17.026549/spectrum['params']['charge'][0])-0.02 <i< spectrum['params']['pepmass'][0]-(17.026549/spectrum['params']['charge'][0])+0.02: #precursor mass - NH4
            ind.append(False)
        else:
            ind.append(True)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    ypoints_data = tic_normalize(ypoints_data)
    
    if len(ypoints_data)>0:
        ind = ypoints_data >1e-5
        xpoints_data = xpoints_data[ind]
        ypoints_data = ypoints_data[ind]
        charge = charge[ind]
    
    if tag==True:
        dev = 20*(spectrum['params']['charge'][0]*pp_mass)/1e6 #20ppm
        if dev <= 0.01:
              dev = 0.01
        # if dev >0.05:#
        #     dev = 0.05
    else:
        dev = 15*(spectrum['params']['charge'][0]*pp_mass)/1e6 #15ppm
        if dev <= 0.01:
              dev = 0.01
        if dev >0.0185:#
            dev = 0.0185
    if len(xpoints_data)> 1 and tag==False:
        plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
        plt.title(title+ ' before mirror')
        plt.show()
    if len(xpoints_data)<=(spectrum['params']['charge'][0]*pp_mass)/15+25 and tag==False:
        print('only ghost performed')
        neutral_lossX, neutral_lossY, mdev = mirror_B(xpoints_data, ypoints_data, pp_mass, 1)
    elif tag==False:
        print('mirroring')
        neutral_lossX, neutral_lossY, mdev = mirror_B(xpoints_data, ypoints_data, pp_mass)
        print(len(neutral_lossX))
    if (spectrum['params']['charge'][0]*pp_mass)/20+30 >= len(xpoints_data):
        print('only ghost performed 2')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass,spectrum, 1, tag=tag)
    else:
        print('mirroring')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass, spectrum, tag=tag)
        print(len(xpoints_data))
    if tag==False and len(neutral_lossY)>0:
        neutral_lossY = tic_normalize(neutral_lossY)
        ind = neutral_lossY >1e-4
        neutral_lossX = neutral_lossX[ind]
        neutral_lossY = neutral_lossY[ind]
    
    ypoints_data = tic_normalize(ypoints_data)
    if len(ypoints_data)>0:
        ind = ypoints_data>1e-4
        xpoints_data = xpoints_data[ind]
        ypoints_data = ypoints_data[ind]
    if len(xpoints_data)>1 and tag==False:
        plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
        plt.title(title+ ' after')
        plt.show()
    if len(set([int(xs) for xs in xpoints_data]))<=pp_mass*spectrum['params']['charge'][0]/100:
        xpoints_data = []
        ypoints_data = []
    if tag==True:
        return xpoints_data, ypoints_data, pp_mass, dev
    spectrum_new = spectrum.copy()
    spectrum_new['intensity array'] = neutral_lossY
    spectrum_new['m/z array'] = neutral_lossX
    spectrum_new['charge array'] = np.ma.MaskedArray([1]*len(neutral_lossX))
    to_mgf = spectrum_new
    return xpoints_data, ypoints_data, pp_mass, dev, to_mgf, neutral_lossX, neutral_lossY
        
def pip_data(xpoints_data, ypoints_data, charge, pp_mass):
    ind = np.argsort(xpoints_data)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    
    xpoints_data = np.array([num for num in xpoints_data if num <= (((pp_mass)*2))-20]) 
    ind = np.argsort(xpoints_data)
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    
    for del_y, i in enumerate(xpoints_data):
        if charge[del_y] >1 and (round(i,2) in [round(num,2) for num in xpoints_data[charge==1]]):
            ypoints_data[del_y] = min(ypoints_data)
        if pp_mass-0.1 <i< pp_mass+0.1:#precursor mass
            ypoints_data[del_y] = min(ypoints_data)
        elif pp_mass-(18.010565/2)-0.01 <i< pp_mass-(18.010565/2)+0.01: #precursor mass - H2O
            ypoints_data[del_y] = min(ypoints_data)
        elif pp_mass-(17.026549/2)-0.01 <i< pp_mass-(17.026549/2)+0.01: #precursor mass - NH4
            ypoints_data[del_y] = min(ypoints_data) 
    ypoints_data = tic_normalize(ypoints_data)
    
    return xpoints_data, ypoints_data

def w_mgf(files):
    mgf.write(files, output='./MGF_mirror_test')

def rank_results(pip_results):
    print('start ranking')
    rank_score = [0 if d < 0.45 else a*b*c*d*(err-abs(e)) for a, b, c, d, e, err in pip_results[['jaccard_score', 'coverage', 'pc_score', 'ms2pip_score', 'difference', 'error_interval']].values]
    rank_score = [0 if (score<0) else score for score in rank_score] #score>1 or 
    pip_results['rank_score'] = np.array(rank_score)
    drop_peps = []
    column_values = ['peptide', 'PTM', 'spectrum_nr', 'charge_state', 'jaccard_score', 'coverage', 'pc_score', 'ms2pip_score', 'rank_score', 'file_name']
    final_df = pd.DataFrame(columns = column_values) 
    for spectra in set(pip_results['spectrum_nr'].values):
        ranking = pip_results['rank_score'][pip_results['spectrum_nr']==spectra].values
        ranking = sorted(ranking)[::-1][:5]
        ranking = ranking[-1]
        peps = pip_results['peptide'][(pip_results['spectrum_nr']==spectra) & (((pip_results['rank_score'] >= ranking)|(pip_results['peptide'][-1]=='K')|(pip_results['peptide'][-1]=='R'))&(pip_results['rank_score'] != 0))].values
        
        for pep in peps:
            array = np.array([[pep, ptm, spc, chs, js, cov, pc,pip, rs, fn] for pep, ptm, spc, chs, js, cov, pc,pip, rs, fn in pip_results[['peptide', 'PTM', 'spectrum_nr', 'charge_state', 'jaccard_score', 'coverage', 'pc_score','ms2pip_score', 'rank_score', 'file_name']][(pip_results['spectrum_nr']==spectra) & (pip_results['peptide']==pep)].values])
            df_add = pd.DataFrame(data = array, columns = column_values)
            final_df = pd.concat([final_df, df_add], ignore_index = True)
            # drop_peps.append(pep)
            # pip_results = pip_results.loc[~((pip_results['peptide'] == pep) & (pip_results['spectrum_nr'] ==spectra))]
    in_db = []
    for seq in final_df['peptide']:
        to_add = []
        for sequence in personal_db.keys():
            if seq in sequence:
                to_add.append(personal_db[sequence])
        for sequence in crap_individual.keys():
            if seq in sequence:
                to_add.append(crap[sequence])
        if len(to_add)>0:
            annot = '|'.join(to_add[0:10])
            in_db.append('Found '+str(len(to_add))+(' hit(s):')+annot)
        else:
            in_db.append('unkown')
    final_df['ID'] = np.array(in_db)
    final_df= final_df.set_index('peptide')
    #pip_results = pip_results.drop(drop_peps)
    return final_df

def find_consensus(ranked_results): #dit nog aanpassen zodat chimeren 2 uitkomen en niet XXXXXXX
    print('making consensus')
    fin_pep = [[num[0][0]]*num[1] for num in ranked_results[0:20]]
    fin_pep_all = []
    for l in fin_pep:
        fin_pep_all = fin_pep_all+l
    
    final_best = []
    for length in range(len(min(fin_pep_all, key=len)), len(max(fin_pep_all, key=len))+1):
        bestseqs = [[]]
        fin_pep = [num for num in fin_pep_all if len(num)==length]
        if len(fin_pep)>0:
            n = len(fin_pep[0])
            profile = {N:[0]*n for N in set(''.join(fin_pep))}
            for seq in fin_pep:
                for i, char in enumerate(seq):
                    profile[char][i] += 1
            for i in range(n):
                d = {N:profile[N][i] for N in set(''.join(fin_pep))}
                if max(d.values())/sum(d.values())>=0.5:
                    m = max(d.values())/2
                    l = [N for N in set(''.join(fin_pep)) if d[N] >= m]
                    bestseqs = [ s+[N] for N in l for s in bestseqs ]
                else:
                    bestseqs = [ s+['X'] for s in bestseqs ]
            final_best.append([''.join(num) for num in bestseqs])
    return final_best 

def filter_peaks(x, y, dev):
    print('length of xpoints was ', len(x))
    keep_iterating = True
    while keep_iterating == True:
        matrix = []
        matrix_row = []
        for element1 in x:
            for element2 in x:
                calc = abs(element1-element2)
                done = False
                for key, val in list(AA_code.items()):
                    if val-dev/2 <= calc <= val+dev/2:
                        matrix_row.append(key)
                        done = True
                        break
                if done==False:
                    matrix_row.append(0)
            matrix.append(matrix_row)
            matrix_row=[]
        matrix= pd.DataFrame(matrix)
        ind = []
        adjust = []
        for col in range(0,len(matrix)):
            if len(set(matrix.loc[col]))>1:
                ind.append(True)
                adjust.append(len(set(matrix.loc[col]))-1)
            else:
                ind.append(False)
        former_x_length = len(x)
        x = x[ind]
        y= y[ind]
        if former_x_length == len(x):
            keep_iterating = False
        print('length of xpoints is now ', len(x))
    # y = np.array([u*adjust[n] for n, u in enumerate(y)])
    return x, y, matrix

def mirror_A(x, y, peptide_mass, spectrum, mdev = 0, takex=[], takey=[], tag=False):
    peptide_mass = peptide_mass*spectrum['params']['charge'][0]
    if peptide_mass>2200 and tag==False:
        return [], [], 0
    if mdev == 1:
        dev = 0.08
        takex = x
        takey = y
    elif mdev == 0:
        takex = x
        takey = y
        dev = 20*peptide_mass/1e6 #20ppm
        if dev <= 0.02:
              dev = 0.02
        if dev >= 0.05:
            dev = 0.05
    else:
        dev = mdev*0.97
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[])#1.0073,19.0226
        mirror_x = sorted([(peptide_mass - u,i) for i,u in enumerate(x_withh2o)])
        new_x =[]
        for locx,element in enumerate(x):
            added = False
            for value,loc in mirror_x:
                if element-dev<=value <= element+dev and ((y[loc]+y[locx])>sorted(y)[int(len(y)*0.6)] or mdev not in [0,1]):
                    new_x.append(True)
                    added = True
                    break
            if added == False:
                new_x.append(False)
        iteration += 1
        xnew = x[new_x]
        ynew = y[new_x]
        if iteration == 1 and (mdev == 0 or mdev == 1):
            x_it1 = xnew
            y_it1 = ynew
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew, spectrum) #peaks that were not accounted for can be double charged

    xnew = list(x_it1) + list(xnew)
    ynew = list(y_it1) + list(ynew)
    xnew= np.array(xnew)
    ynew=np.array(ynew)
    ind = np.argsort(xnew)
    xnew = xnew[ind]
    ynew = ynew[ind]
    if len(xnew)==0:
        return xnew, ynew, mdev
    
    iterating = True
    finalx = []
    finaly = []
    while iterating == True:
        keepx = []
        keepy = []
        already = []
        for i in range(0,len(xnew)-1):
            if xnew[i+1]-xnew[i] >dev and xnew[i] not in already:#0.01 terug zetten????
                keepx.append(xnew[i])
                keepy.append(ynew[i])
            elif xnew[i] not in already:
                keepx.append(xnew[i])
                keepy.append(ynew[i])
                already.append(xnew[i+1])
        keepx.append(xnew[-1])
        keepy.append(ynew[-1])
        if len(finalx)==len(keepx):
            iterating = False
        finalx = keepx
        finaly = keepy
        xnew = keepx
        ynew = keepy
    already = []
    finalx = np.array(finalx)
    finaly = np.array(finaly)
    if len(finalx)>min(peptide_mass/20+30,175) and tag==False:
        finalx, finaly, mdev = mirror_A(finalx, finaly,peptide_mass/spectrum['params']['charge'][0], spectrum, dev, takex, takey, tag=tag)
    elif len(finalx)>peptide_mass/10 and tag==True:
        finalx, finaly, mdev = mirror_A(finalx, finaly,peptide_mass/spectrum['params']['charge'][0], spectrum, dev, takex, takey, tag=tag)
    elif peptide_mass/80 >= len(finalx) and peptide_mass>900:
        return [], [], mdev
    elif peptide_mass/80+1000 <= len(finalx)<=max(peptide_mass/20+30,175) and tag==False:
        finalx = list(finalx)
        finaly = list(finaly)
        for loc, i in enumerate(takex):
            if peptide_mass >1200 and min(extra_AA.values())<= i <= max(extra_AA.values())+20 and i not in finalx:
                for to_print, aa in extra_AA.items():
                    for io in [19.0226, 1.0073]:
                        aa = aa+io
                        if i-0.018 <= aa <= i+0.018 and loc not in already:
                            adding = True
                            for a in finalx:
                                for XXX in [18.010565, 17.026549]:
                                    if a+XXX-0.018 <= i <= a+XXX+0.018 and io==1.0073:
                                        adding = False
                            if adding == True:
                                already.append(loc)
                                finalx.append(takex[loc])
                                finaly.append(takey[loc])
                                print('hohoho eentje gemist', to_print, '+', io)
                                break
        finalx = np.array(finalx)
        finaly = np.array(finaly)
    ind = np.argsort(finalx)
    finalx = finalx[ind]
    finaly = finaly[ind]
    return finalx, finaly, mdev

def mirror_B(x,y,peptide_mass, mdev = 0, takex=[], takey=[]):
    peptide_mass = peptide_mass*spectrum['params']['charge'][0]
    if peptide_mass>2200:
        return [], [], 0
    if mdev == 1:
        dev = 0.1
        takex = x
        takey = y
    elif mdev == 0:
        takex = x
        takey = y
        dev = 20*peptide_mass/1e6 #20ppm
        if dev <= 0.025:
              dev = 0.025
        if dev >= 0.05:
            dev = 0.05
    else:
        dev = mdev*0.97
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[])#1.0073,19.0226
        mirror_x = sorted([(peptide_mass - u,i) for i,u in enumerate(x_withh2o)])
        new_x =[]
        for locx,element in enumerate(x):
            added = False
            for value,loc in mirror_x:
                if element-dev<=value <= element+dev and ((y[loc]+y[locx])>sorted(y)[int(len(y)*0.6)] or mdev not in [0,1]):
                    new_x.append(True)
                    added = True
                    for XXX in [18.010565, 17.026549]:
                        if element-XXX-dev <= element-XXX <= value-XXX+dev:
                            for element2 in x:
                                if element -XXX-dev <= element2 <= element-XXX+dev:
                                    adj = list(x).index(element2)
                                    new_x[adj]=True
                    break
            if added == False:
                new_x.append(False)
        iteration += 1
        xnew = x[new_x]
        ynew = y[new_x]
        if iteration == 1 and (mdev == 0 or mdev == 1):
            x_it1 = xnew
            y_it1 = ynew
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew, spectrum) #peaks that were not accounted for can be double charged

    xnew = list(x_it1) + list(xnew)
    ynew = list(y_it1) + list(ynew)
    xnew= np.array(xnew)
    ynew=np.array(ynew)
    ind = np.argsort(xnew)
    xnew = xnew[ind]
    ynew = ynew[ind]
    if len(xnew)==0:
        return xnew, ynew, mdev
    finalx = xnew
    finaly = ynew
    if len(finalx)>max(peptide_mass/15+25,200):
        finalx, finaly, mdev = mirror_B(np.array(finalx), np.array(finaly),peptide_mass/spectrum['params']['charge'][0], dev, takex, takey)
    xnew = list(finalx)
    ynew = list(finaly)
    
    iterating = True
    finalx = []
    finaly = []
    while iterating == True:
        keepx = []
        keepy = []
        already = []
        for i in range(0,len(xnew)-1):
            if xnew[i+1]-xnew[i] >dev and xnew[i] not in already:#0.01
                keepx.append(xnew[i])
                keepy.append(ynew[i])
            elif xnew[i] not in already:
                keepx.append(xnew[i])
                keepy.append(ynew[i])
                already.append(xnew[i+1])
        keepx.append(xnew[-1])
        keepy.append(ynew[-1])
        if len(finalx)==len(keepx):
            iterating = False
        finalx = keepx
        finaly = keepy
        xnew = keepx
        ynew = keepy
    finalx = np.array(finalx)
    finaly = np.array(finaly)
    ind = np.argsort(finalx)
    finalx = finalx[ind]
    finaly = finaly[ind]
    return finalx, finaly, mdev

def ghost_peaks(x, y, peptide_mass, xnew, ynew, spectrum):
    peptide_mass = peptide_mass/spectrum['params']['charge'][0]
    x = list(x)
    y=list(y)
    xnew = list(xnew)
    ynew = list(ynew)
    ghost = [num*2-1 for i, num in enumerate(x) if max(extra_AA.values())<=num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew] #+(4*(spectrum['params']['charge'][0]-2))
    xghost = [num for i, num in enumerate(x) if max(extra_AA.values())<=num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew]
    ghost_intensity = [num for i, num in enumerate(y) if max(extra_AA.values())<=x[i]<=peptide_mass*(spectrum['params']['charge'][0]/2) and x[i] not in xnew]
    ghost = list(np.array(ghost))
    x = xghost+ghost
    
    y=ghost_intensity+ghost_intensity
    x= np.array(x)
    y=np.array(y)
    
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    return x, y

def IL_adjusted(cbl):
    cbl_fin=[]
    for element2 in cbl:
        element = element2[0]
        I_list = [element]
        while 'L' in element:
            x = element.find('L')
            temp_Ilist = []
            for i in I_list:
                new_element = i[:x]+'I'+i[x+1:]
                temp_Ilist.append(new_element)
            I_list = I_list + temp_Ilist
            element = element[:x]+'I'+element[x+1:]
        cbl_fin = cbl_fin+ [(y, element2[1], element2[2], element2[3], element2[4]) for y in I_list]
    return cbl_fin

def deamidation(cpl):
    #Deamidation NQ
    cpl_DN=[]
    for element2 in cpl:
        element = element2[0]
        ptm = element2[1]
        I_list = [(element, ptm)]
        while 'D' in element:
            x = element.find('D')
            temp_Ilist = []
            for i, ptm in I_list:
                new_element = i[:x]+'N'+i[x+1:]
                if len(ptm)==0:
                    new_ptm = str(x+1)+'|Deamidated'
                elif len(ptm.split('|'))==2 or x+1<= int(ptm.split('|')[0]) or x+1>= int(ptm.split('|')[-2]):
                    if x+1 >= int(ptm.split('|')[-2]):
                        new_ptm = ptm+'|'+str(x+1)+'|Deamidated'
                    else:
                        new_ptm = str(x+1)+'|Deamidated|'+ptm
                else:
                    new_ptm = ptm.split('|')
                    for n in range(0,len(new_ptm),2):
                        if int(new_ptm[n]) <= x+1 and int(new_ptm[n+2]) >= x+1:
                            new_ptm.insert(n+2, str(x+1)+'|Deamidated')
                            new_ptm = '|'.join(new_ptm)
                            break
                temp_Ilist.append((new_element, new_ptm))
            I_list = I_list + temp_Ilist
            element = element[:x]+'N'+element[x+1:]
        cpl_DN = cpl_DN+ [(y, ptms, element2[2], element2[3], element2[4]) for y, ptms in I_list]
    
    cpl_fin=cpl_DN
    for element2 in cpl_DN:
        element = element2[0]
        ptm = element2[1]
        I_list = [(element, ptm)]
        while 'E' in element:
            x = element.find('E')
            temp_Ilist = []
            for i, ptm in I_list:
                new_element = i[:x]+'Q'+i[x+1:]
                if len(ptm)==0:
                    new_ptm = str(x+1)+'|Deamidated'
                elif len(ptm.split('|'))==2 or x+1<= int(ptm.split('|')[0]) or x+1>= int(ptm.split('|')[-2]):
                    if x+1 >= int(ptm.split('|')[0]):
                        new_ptm = ptm+'|'+str(x+1)+'|Deamidated'
                    else:
                        new_ptm = str(x+1)+'|Deamidated|'+ptm
                else:
                    new_ptm = ptm.split('|')
                    for n in range(0,len(new_ptm),2):
                        if int(new_ptm[n]) <= x+1 and int(new_ptm[n+2]) >= x+1:
                            new_ptm.insert(n+2, str(x+1)+'|Deamidated')
                            new_ptm = '|'.join(new_ptm)
                            break
                temp_Ilist.append((new_element, new_ptm))
            I_list = I_list + temp_Ilist
            element = element[:x]+'Q'+element[x+1:]
        cpl_fin = cpl_fin+ [(y, ptms, element2[2], element2[3], element2[4]) for y, ptms in I_list if (y, ptms, element2[2], element2[3], element2[4]) not in cpl_DN]
    return cpl_fin

def annotate(combined_peptide_list, end_result):
    found = False
    for a in combined_peptide_list:
        x = a
        a = [a[0].upper(), 'test']
        if a[0] in ['AETSELHTSLK', 'GAYVEVTAK', 'LGNEQGVSR', 'LVGTPAEER', 'LDSTSLPVAK', 'AGLLVAEGVTK',
                      'LGLDFDSFR', 'GFTAYYLPR', 'SGGLLWQLVR', 'AVGANPEQLTR', 'SAEGLDASASLR',
                      'VFTPELVDVAK','VGNELQYVALR', 'YLELAPGVDNSK', 'DGTFAVDGPGVLAK', 'YDSLNNTEVSGLR',
                      'SPYVLTGPGVVEYK', 'ALENDLGVPSDATVK', 'AVYFYAPQLPLYANK', 'TVESLFPEEAETPGSAVR',  'STQAALDQLNGK', 'ALLVASGHLK']:
                print('YESSSSSS!!!!!!!!!')
                end_result[title+'_found']=['YES', a[0], len(combined_peptide_list)]
                found = True
        for seq in crap.keys():
            if a[0] in seq.replace('I', 'L'):
                print('Found', crap[seq])
                found = True
                print(x)
                end_result[title+crap[seq]+a[0]]=(crap[seq], x)
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                break
        for seq in personal_db.keys():
            if a[0] in seq.replace('I', 'L'):
                print('Found', personal_db[seq])
                found = True
                print(x)
                end_result[title+personal_db[seq]+a[0]]=(personal_db[seq], x)
                print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                break
    return end_result, found

def tag_spectrum(y_peptide, b_peptide, spectrum, intensity, ptm, dev):
    ptm_dict = {}
    ptm = ptm.split('|')
    if len(ptm)>1:
        for i in range(0,len(ptm),2):
            if ptm[i] not in [0,-1]:
                ptm_dict[len(y_peptide)-int(ptm[i])] = ptm[i+1]
            elif ptm[i]==0:
                ptm_dict[len(y_peptide)-1]=ptm[i+1]
            elif ptm[i]==-1:
                ptm_dict[-1]=ptm[i+1]  
    if -1 in ptm_dict:
        yextra = float(unimod[ptm_dict[-1]])
    else:
        yextra = 0
    x_array_Y = []
    x_point_Y = ion_types['Y'] + yextra
    for loc, AA in enumerate(y_peptide):
        if loc in ptm_dict:
            adding = float(unimod[ptm_dict[loc]])
        else:
            adding = 0
        x_array_Y.append(AA_codes_tag[AA]+x_point_Y+adding)
        x_point_Y += AA_codes_tag[AA]+adding
    ptm_dict = {}
    if len(ptm)>1:
        for i in range(0,len(ptm),2):
            if ptm[i] not in [0,-1]:
                ptm_dict[int(ptm[i])-1] = ptm[i+1]
            elif ptm[i]==0:
                ptm_dict[-1]=ptm[i+1]
            elif ptm[i]==-1:
                ptm_dict[len(ptm)-1]=ptm[i+1]   
    if -1 in ptm_dict:
        bextra = float(unimod[ptm_dict[-1]])
    else:
        bextra = 0
    x_array_B = []
    x_point_B = ion_types['B'] + bextra
    for loc, AA in enumerate(b_peptide):
        if loc in ptm_dict:
            adding = float(unimod[ptm_dict[loc]])
        else:
            adding = 0
        x_array_B.append(AA_codes_tag[AA]+x_point_B+adding)
        x_point_B += AA_codes_tag[AA]+adding 
    x_array = set(x_array_B + x_array_Y) 
    x_array = sorted(x_array)
    pep_spectrum ={}
    for loc,num in enumerate(spectrum):
        former = 0
        for value in x_array:
            if value-dev <= num <= value+dev:
                if intensity[loc]>former:
                    pep_spectrum[value]=num
    final_y_array = []
    for y_loc, point in enumerate(spectrum):
        if point in pep_spectrum.values():
            final_y_array.append(intensity[y_loc])
        else:
            final_y_array.append(1e-15)
    return spectrum, np.array(final_y_array)

def deamidated_tag(tag_list):
    output_E = []
    for pep in tag_list:
        if 'E' not in pep[1]:#D->N
            output_E.append(pep)
        elif 'E' in pep[1]:
            new=[pep[1]]
            element = pep[1]
            while 'E' in element:
                x = element.find('E')
                temp_list = []
                for i in new:
                    new_element = i[:x]+'Q'+i[x+1:]
                    temp_list.append(new_element)
                new = new + temp_list
                element = element[:x]+'Q'+element[x+1:]
            output_E = output_E+ [(pep[0], y, pep[2]-0.984016, pep[3]) if pep[1] != y else (pep[0], y, pep[2], pep[3]) for y in new]
    output_D = []
    for pep in output_E:
        if 'D' not in pep[1]:
            output_D.append(pep)
        elif 'D' in pep[1]:
            new=[pep[1]]
            element = pep[1]
            while 'D' in element:
                x = element.find('D')
                temp_list = []
                for i in new:
                    new_element = i[:x]+'N'+i[x+1:]
                    temp_list.append(new_element)
                new = new + temp_list
                element = element[:x]+'N'+element[x+1:]
            output_D = output_D+ [(pep[0], y, pep[2]-0.984016, pep[3]) if pep[1] != y else (pep[0], y, pep[2], pep[3]) for y in new]
    return output_D

def do_tag(reader):
    column_values = ['peptide', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
    df_deeplc = pd.DataFrame(columns = column_values)
    nr_file =0
    for spectrum in reader:
        nr_file += 1
        if nr_file >180:
            continue
        rt_time = float(spectrum['params']['rtinseconds'])
        if 'charge' in spectrum['params']:
            charging = spectrum['params']['charge']
        else:
            charging  = [2]
        for test_charge in charging:
            spectrum['params']['charge']=[test_charge]
            sp_nr = str(nr_file)+'_'+str(test_charge)
            title = 'Spectrum_' + str(nr_file)+'_'+str(test_charge)
            print(title)
            xpoints_data, ypoints_data, pp_mass, dev = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'], spectrum, tag=True)
            print('Making TAG...')
            TAG_list = make_TAG(xpoints_data, dev, pp_mass, test_charge,ypoints_data)
            TAG_list = sorted(TAG_list, key=lambda x:x[3])[::-1]#len(x[1])
            if len(TAG_list)==0:
                continue
            # TAG_list = TAG_list[:2500]
            
            print('Combining TAG...')
            TAG_list = combos(TAG_list, pp_mass*test_charge)
            #TAG_list = deamidated_tag(TAG_list)
            TAG_list = reduce_taglist(TAG_list)
            # TAG_list = sorted(TAG_list, key=lambda x:x[4])[::-1]#x:len(x[1])
            # TAG_list = TAG_list[:500]
            print('completing sequence...')
            result = find_seq(crap, TAG_list)
            result = locate_seq(result, crap)
            ptm_check = [True for num in result if 0 in num[1] and ('K' in num[0][-1] or 'R' in num[0][-1]) and num[2]=='contamination']
            if len(ptm_check) > 0:
                print('here')
                result = [num for num in result if 0 in num[1] and ('K' in num[0][-1] or 'R' in num[0][-1]) and num[2]=='contamination']
            if len(result)>0:
                print(result)
                if 'charge' in spectrum['params'].keys():
                    charges = spectrum['params']['charge'][0]
                else:
                    charges = 2
                array = np.array([[pep[0], pep[1], 0.95, 0.95, 0.95, rt_time,xpoints_data, ypoints_data, sp_nr, charges, spectrum['params']['title']] for pep in result], dtype='object')
                column_values = ['peptide','PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
                df_add = pd.DataFrame(data = array, columns = column_values)
                df_deeplc = pd.concat([df_deeplc, df_add], ignore_index = True)
    ptm_list = []
    for PTM in df_deeplc['PTM'].values:
        if PTM == [0]:
            ptm_list.append('')
        else:
            ptm_list.append(PTM)
    df_deeplc['PTM']=ptm_list
    df_deeplc = ptm_combos(df_deeplc)
    
    PC_list =[]
    j_list = []
    cov_list = []
    for item, ptm, x, y in df_deeplc[['peptide', 'PTM', 'xpoints_data', 'ypoints_data']].values:
        a, c = tag_spectrum(item[::-1], item, x, y, ptm, max(20*find_mass(item)/1e6,0.01))
        j_list.append(jaccard(y, c, item))
        PC_list.append(pearson_correlation(y, c)[0][1])
        found_intensity = np.sum(c)
        all_intensity = np.sum(y)
        cov_list.append(found_intensity/all_intensity)
    df_deeplc['pc_score']=np.array(PC_list)
    df_deeplc['jaccard_score']=np.array(j_list)
    df_deeplc['coverage']=np.array(cov_list)
    keep = []
    for loc,i in enumerate(df_deeplc[['peptide', 'pc_score']].values):
        if i[0].endswith('R') or i[0].endswith('K') or float(i[1])>0.85:
            keep.append(loc)
    df_deeplc = df_deeplc.iloc[keep]
    df_deeplc = df_deeplc.reset_index()
    del df_deeplc['index']
    
    df_deeplc = df_deeplc.loc[df_deeplc['pc_score']>0.6]
    # df_deeplc = df_deeplc.loc[df_deeplc['coverage']>0.4]
    df_deeplc = df_deeplc.reset_index()
    del df_deeplc['index']
    keep = []
    for i in set(df_deeplc['spectrum_nr'].values):
        maximum = max(df_deeplc['jaccard_score'][df_deeplc['spectrum_nr']==i].values)
        keep.append(list(np.where(df_deeplc["jaccard_score"] == maximum)[0])[0])
    df_deeplc = df_deeplc.iloc[keep]
    df_deeplc = df_deeplc.reset_index()
    del df_deeplc['index']
    df_deeplc['pc_score']=np.array([1]*len(df_deeplc))
    return df_deeplc

def find_seq(crap, reform_peptide):
    final_result = {}
    for peptides in reform_peptide:
        peptide = peptides[1]
        for sequence, name in crap.items():
            if peptide in sequence:
                if str(name) in final_result.keys():
                    final_result[str(name)] = final_result[str(name)]+[peptides]
                else:
                    final_result[str(name)] = [peptides]
    return final_result

def reduce_taglist(TAG_list):#looking for smallest
    TAG_list = sorted(TAG_list, key=lambda x:len(x[1]))
    return_tag = []
    smallest = ['X']
    for element in TAG_list:
        add = True
        for seq in smallest:
            if seq in element[1]:
                add = False
        if add == True:
            smallest.append(element[1])
            return_tag.append(element) 
    print('reducing ',len(TAG_list), ' to ', len(return_tag))
    return return_tag

def make_TAG(x, dev, mass, charge,y):
    x = list(ion_types.values())+list(x)
    y = [0,0]+list(y)
    x = np.array(x)
    done = []
    y=np.array(y)
    min_tag = min(math.ceil(mass*charge/500)-1,4)
    min_tag = max(min_tag,3)#enkel voor histonen 2
    print('minimal tag =', min_tag)
    tag_list = []
    massa = mass*charge
    start_time = time.time()
    for begin in x[x<massa*75]:
        # for aa1 in AA_codes_tag.keys():
        iterate = True
        temp_list = [('',0)]
        if len([num for num in tag_list if len(num[1])>min_tag])>2000:
            break
        while iterate==True and abs(start_time - time.time())<60:
            new_list = []
            for peps in temp_list:
                pep = peps[0]
                cov = peps[1]
                if iterate == False:
                    break
                for aa2 in AA_codes_tag.keys():
                    new_pep = pep+aa2
                    former = True
                    if len(new_pep)>min(min_tag+3,6):
                        iterate = False
                        break
                    new_y = y[x>begin]
                    new_x = x[x>begin]
                    for z,xp in enumerate(new_x):
                        if xp-dev <= begin+find_mass(new_pep) <= xp+dev:
                            new_list.append((new_pep,cov+new_y[z]))
                            former = False
                            break
                    if former == True and len(pep)>min_tag and (begin,pep,(mass*charge)-(begin+find_mass(pep)),cov) not in tag_list and begin+find_mass(pep)<mass*charge and pep+str(begin) not in done:
                        tag_list.append((begin,pep,(mass*charge)-(begin+find_mass(pep)),find_mass(pep)))
                        done.append(pep+str(begin))
            if len(new_list) == 0:
                iterate = False
            else:
                temp_list = new_list
        temp_list = [(begin,element[0],((mass*charge)-begin-find_mass(element[0])),find_mass(element[0])) for element in temp_list if len(element[0])>min_tag]
        tag_list = tag_list+temp_list
    if abs(start_time - time.time())>60:
        print('out of time')
    return [num for num in set(tag_list)]

def combos(tag, mass):
    combo = []
    for tag1 in tag:
        tag_seq1 = tag1[1]
        if tag1[0]==ion_types['Y']:
            continue
        for tag2 in tag:
            tag_seq2 = tag2[1][::-1]
            for i in range(0,min(len(tag_seq1), len(tag_seq2))-1):
                
                if tag_seq1[-(i+1):]==tag_seq2[:i+1]:
                    if mass-0.018<=find_mass(tag_seq1+tag_seq2[:i])+tag1[0]+tag2[2]<=mass+0.018:
                        if (tag1[0],tag_seq1+tag_seq2[:i],tag2[2],0) not in combo and len(tag_seq1+tag_seq2[:i])>=min(len(tag_seq1),len(tag_seq2)):
                            combo.append((tag1[0],tag_seq1+tag_seq2[:i],tag2[2], 0))
                            break
            # m1 = find_mass(tag_seq1)
            # m2 = find_mass(tag_seq2)
            # if mass-1<=tag1[0]+tag2[0]+m2+m1+abs((tag1[2]-tag2[2]))<=mass+1 and tag2[0]>ion_types['Y'] and tag1[0]!=ion_types['Y']:
                
            #     for aa, el in AA_code.items():
            #         for ptm, num in unimod.items():
            #             if abs((tag1[2]-tag2[2]))-0.1<=float(el)+float(num)<=abs((tag1[2]-tag2[2]))+0.1 and tag2[0]>ion_types['Y'] and tag1[0]!=ion_types['Y']:
            #                 if (tag1[0],tag_seq1+aa+tag_seq2[::-1],tag2[0],ptm) not in combo and len(tag_seq1+aa+tag_seq2[::-1])>4:
            #                         combo.append((tag1[0],tag_seq1+aa+tag_seq2[::-1],tag2[0], 'between_tags_on_'+aa+'_with_'+str(num)))
            #                         break
    return combo

def locate_seq(tag, db):
    result = []
    crap_t = {str(u):i for i,u in db.items()}
    for name, peps in tag.items():
        seq = crap_t[name]
        for el in peps:
            match = re.finditer(str(el[1]), str(seq))
            for i in match:
                peptide = el[1]
                ptm=[el[3]]
                start=i.span()[0]
                end=i.span()[1]
                peptide_N_plus = ''
                ptm_1 = ptm
                if el[0] != ion_types['B']:
                    it = True
                    mass = el[0]
                    added = ''
                    while start >2 and it == True:
                        if seq[start-1]!='-' and seq[start-1]!='X' and mass - find_mass(seq[start-1])>0:
                            peptide = seq[start-1]+peptide
                            added = seq[start-1]+added
                            mass -= find_mass(seq[start-1])
                            start -= 1
                        elif seq[start-1]!='-' and seq[start-1]!='X' and 0.05>mass - find_mass(seq[start-1])>-0.05:
                            peptide = seq[start-1]+peptide
                            added = seq[start-1]+added
                            it=False
                        elif seq[start-1]!='-' and seq[start-1]!='X' and mass - find_mass(seq[start-1])==0:
                            it = False
                        else:
                            if ptm != [0] and abs(mass-ion_types['B'])>0.05:
                                ptm = ptm+[('N',control_ptm(mass-ion_types['B'], seq[start-1],added, 'N-term'))]
                                peptide_N_plus = seq[start-1]+peptide
                                ptm_1=[el[3]]+[('N',control_ptm(mass-ion_types['B']-find_mass(seq[start-1]), seq[start-2],seq[start-1]+added, 'N-term'))]
                            elif abs(mass-ion_types['B'])>0.05:
                                ptm = [('N',control_ptm(mass-ion_types['B'], seq[start-1],added, 'N-term'))]
                                peptide_N_plus = seq[start-1]+peptide
                                ptm_1=[('N',control_ptm(mass-ion_types['B']-find_mass(seq[start-1]), seq[start-2],seq[start-1]+added, 'N-term'))]
                            it = False
                        
                    if 0<start<=2:
                        if ptm != [0] and abs(mass-ion_types['B'])>0.05:
                            ptm = ptm+[('N',control_ptm(mass-ion_types['B'], seq[start],added, 'N-term'))]
                            peptide_N_plus = seq[start-1]+peptide
                            ptm_1=[el[3]]+[('N',control_ptm(mass-ion_types['B']-find_mass(seq[start-1]), seq[start-2],seq[start-1]+added, 'N-term'))]
                        elif abs(mass-ion_types['B'])>0.05:
                            ptm = [('N',control_ptm(mass-ion_types['B'], seq[start],added, 'N-term'))]
                            peptide_N_plus = seq[start-1]+peptide
                            ptm_1=[('N',control_ptm(mass-ion_types['B']-find_mass(seq[start-1]), seq[start-2],seq[start-1]+added, 'N-term'))]
                if len(peptide_N_plus)==0:
                    peptide_N_plus = peptide
                    ptm_1 = ptm
                peptide_N_plus_C = ''
                peptide_N_plus_C_plus = ''
                peptide_C_plus = ''
                if el[2] != ion_types['Y']:
                    it = True
                    mass = el[2]
                    added = ''
                    while end <len(seq)-2 and it == True:
                        if seq[end] != '-' and seq[end] != 'X' and mass - find_mass(seq[end])>0:
                            peptide = peptide+seq[end]
                            added = added+seq[end]
                            mass -= find_mass(seq[end])
                            end += 1
                        elif seq[end] != '-' and seq[end] != 'X' and 0.05>mass - find_mass(seq[end])>-0.05:
                            peptide = peptide+seq[end]
                            added = added+seq[end]
                            it=False
                        elif seq[end] != '-' and seq[end] != 'X' and mass - find_mass(seq[end])==0:
                            it = False
                        else:
                            if ptm != [0] and abs(mass-ion_types['Y'])>0.05:
                                ptm = ptm+[('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term',Clength=len(peptide)-len(added)))]
                                peptide_C_plus = peptide+seq[end+1]
                                ptm_2=ptm+[('C',control_ptm(mass-ion_types['Y']-find_mass(seq[end+1]), seq[end+2],added, 'C-term',Clength=len(peptide)-len(added+seq[end+1])))]
                            elif abs(mass-ion_types['Y'])>0.05:
                                ptm = [('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term',Clength=len(peptide)-len(added)))]
                                peptide_N_plus_C = peptide_N_plus + added
                            if ptm_1 != [0] and abs(mass-ion_types['Y'])>0.05:
                                peptide_N_plus_C = peptide_N_plus + added
                                ptm_1=ptm_1+[('C',control_ptm(mass-ion_types['Y']-find_mass(peptide[0]), seq[end+1],added, 'C-term',Clength=len(peptide_N_plus)-len(added)))]
                                peptide_N_plus_C_plus = peptide_N_plus+added+seq[end+1]
                                ptm_3=ptm_1+[('N',control_ptm(mass-ion_types['B']-find_mass(seq[end+1])-find_mass(peptide[0]), seq[end+2],added+seq[end+1],'C-term',Clength=len(peptide_N_plus)-len(added+seq[end+1])))]
                            elif abs(mass-ion_types['Y'])>0.05:
                                peptide_N_plus_C = peptide_N_plus + added
                                ptm_1=[('C',control_ptm(mass-ion_types['Y']-find_mass(peptide[0]), seq[end+1],added, 'C-term',Clength=len(peptide_N_plus)-len(added)))]
                                peptide_N_plus_C_plus = peptide_N_plus+added+seq[end+1]
                                ptm_3=[('N',control_ptm(mass-ion_types['B']-find_mass(seq[end+1])-find_mass(peptide[0]), seq[end+2],added+seq[end+1],'C-term',Clength=len(peptide_N_plus)-len(added+seq[end+1])))]
                            it = False
                    if end+2 == len(seq)-1:
                        if ptm != [0] and abs(mass-ion_types['Y'])>0.05:
                            ptm = ptm+[('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term',Clength=len(peptide)-len(added)))]
                            peptide_C_plus = peptide+seq[end+1]
                            ptm_2=ptm+[('C',control_ptm(mass-ion_types['Y']-find_mass(seq[end+1]), seq[end+2],added, 'C-term',Clength=len(peptide)-len(added+seq[end+1])))]
                        elif abs(mass-ion_types['Y'])>0.05:
                            ptm = [('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term',Clength=len(peptide)-len(added)))]
                            peptide_N_plus_C = peptide_N_plus + added
                        if ptm_1 != [0] and abs(mass-ion_types['Y'])>0.05:
                            peptide_N_plus_C = peptide_N_plus + added
                            ptm_1=ptm_1+[('C',control_ptm(mass-ion_types['Y']-find_mass(peptide[0]), seq[end+1],added, 'C-term',Clength=len(peptide_N_plus)-len(added)))]
                            peptide_N_plus_C_plus = peptide_N_plus+added+seq[end+1]
                            ptm_3=ptm_1+[('N',control_ptm(mass-ion_types['B']-find_mass(seq[end+1])-find_mass(peptide[0]), seq[end+2],added+seq[end+1],'C-term',Clength=len(peptide_N_plus)-len(added+seq[end+1])))]
                        elif abs(mass-ion_types['Y'])>0.05:
                            peptide_N_plus_C = peptide_N_plus + added
                            ptm_1=[('C',control_ptm(mass-ion_types['Y']-find_mass(peptide[0]), seq[end+1],added, 'C-term',Clength=len(peptide_N_plus)-len(added)))]
                            peptide_N_plus_C_plus = peptide_N_plus+added+seq[end+1]
                            ptm_3=[('N',control_ptm(mass-ion_types['B']-find_mass(seq[end+1])-find_mass(peptide[0]), seq[end+2],added+seq[end+1],'C-term',Clength=len(peptide_N_plus)-len(added+seq[end+1])))]
                if (peptide, ptm, name) not in result and ('N','MISS') not in ptm and ('C','MISS') not in ptm:
                    result.append((peptide, ptm, name))
                if len(peptide_N_plus_C)>0:
                    if (peptide_N_plus_C, ptm_1, name) not in result and ('N','MISS') not in ptm_1 and ('C','MISS') not in ptm_1:
                        result.append((peptide_N_plus_C, ptm_1, name))
                if len(peptide_N_plus_C_plus)>0:
                    if (peptide_N_plus_C_plus, ptm_3, name) not in result and ('N','MISS') not in ptm_3 and ('C','MISS') not in ptm_3:
                        result.append((peptide_N_plus_C_plus, ptm_3, name))
                if len(peptide_C_plus)>0:
                    if (peptide_C_plus, ptm_2, name) not in result and ('N','MISS') not in ptm_2 and ('C','MISS') not in ptm_2:
                        result.append((peptide_C_plus, ptm_2, name))
    return result

def control_ptm(miss_mass, next_AA, added_AA, term, stop=False,Clength=0): #C-term and N-term separately, different for middle!
    output = []
    reduced_unimod = unimod_db[unimod_db['AA'].isin(list(added_AA)+[term])]
    for ptm1, massa1, AA1 in reduced_unimod[['PTM', 'mass', 'AA']].values:#1ptm
        if miss_mass-0.015<=massa1<=miss_mass+0.015:
            if AA1 not in ['N-term', 'C-term']:
                found = re.finditer(AA1, added_AA)
                for i in found:
                    final_ptm = str(int(i.span()[0])+Clength+1)+'|'+ptm1
            elif AA1 == 'N-term':
                final_ptm = str(0)+'|'+ptm1
            elif AA1 == 'C-term':
                final_ptm = str(-1)+'|'+ptm1
            
            if final_ptm not in output:
                output.append(final_ptm)
    if len(output)>0:
        return output
    for ptm1, massa1, AA1 in reduced_unimod[['PTM', 'mass', 'AA']].values:#2ptm
        for ptm2, massa2, AA2 in reduced_unimod[['PTM', 'mass', 'AA']].values:
            if (ptm1, massa1, AA1)==(ptm2, massa2, AA2) or (added_AA.count(AA1)<2 and AA1 == AA2):
                continue
            if miss_mass-0.01<=massa1+massa2<=miss_mass+0.01:
                if AA1 not in ['N-term', 'C-term'] and AA2 not in ['N-term', 'C-term']:
                    final_ptm = str(int(re.search(AA1, added_AA).span()[0])+Clength+1)+'|'+ptm1+'|'+str(int(re.search(AA2, added_AA).span()[0])+Clength+1)+'|'+ptm2
                elif AA1 == 'N-term':
                    final_ptm = str(0)+'|'+ptm1+'|'+str(int(re.search(AA2, added_AA).span()[0])+Clength+1)+'|'+ptm2
                elif AA1 == 'C-term':
                    final_ptm = str(-1)+'|'+ptm1+'|'+str(int(re.search(AA2, added_AA).span()[0])+Clength+1)+'|'+ptm2
                elif AA2 == 'N-term':
                    final_ptm = str(int(re.search(AA1, added_AA).span()[0])+Clength+1)+'|'+ptm1+'|'+str(0)+'|'+ptm2
                elif AA2 == 'C-term':
                    final_ptm = str(int(re.search(AA1, added_AA).span()[0])+Clength+1)+'|'+ptm1+'|'+str(-1)+'|'+ptm2
                     
                if final_ptm not in output:
                    output.append(final_ptm)
        break
    if len(output)==0:
        output = 'MISS'
    return output

def make_unimod(include_label):
    raw_unimod = pd.DataFrame()
    with open('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/unimod.txt') as r:
        for line in r:
            line = line.replace('=', ',')
            line = line.strip()
            array = np.array(line.split(',')).reshape(1,-1)
            df_add = pd.DataFrame(array)
            raw_unimod = pd.concat([raw_unimod, df_add], ignore_index=True)
    AAs = []
    for aa,locaa in raw_unimod[[1,4]].values:
        
        if '[' in aa:
            i = aa.split(']')
        else:
            AAs.append(aa)
            continue
        if ''.join(c for c in i[1] if c.isdigit()==False) != i[1]:
            add = ''.join(c for c in i[1] if c.isdigit()==False)+'_!'#can change ! with locaa
            while add in AAs:
                add = add+'!'
                print(add)
            AAs.append(add)
        else:
            AAs.append(i[1])
    raw_unimod[1]=np.array(AAs)
    temp_unimod = pd.DataFrame()
    temp_unimod['PTM']=raw_unimod[1].values[2:]
    temp_unimod['delta_mass']=raw_unimod[2].values[2:]
    unimod = {}
    for i,u in temp_unimod[['PTM', 'delta_mass']].values:
        #i = i.split(']')
        # if abs(float(u)) > 120:
        #     continue
        # if ''.join(c for c in i[-1] if c.isdigit()==False) != i[-1]:
        #     continue
        unimod[i]=u

    unimod['AA']='0'
    for i, u in AA_codes_tag.items():
        unimod[i]=str(u)
    unimod_db = raw_unimod.iloc[2:]
    unimod_db = unimod_db.drop([0,3], axis=1)
    unimod_db.columns = ['PTM', 'mass', 'AA', 'type']
    unimod_db['mass']=np.array([float(m) for m in unimod_db['mass'].values])
    #unimod_db=unimod_db[abs(unimod_db['mass'])<max(AA_codes.values())]
    unimod_db['PTM']=np.array([num.split(']')[-1] for num in unimod_db['PTM'].values])
    unimod_db=unimod_db[unimod_db['PTM'] != '']
    if include_label == False:
        unimod_db=unimod_db[unimod_db['type'] != 'Isotopic label']
    unimod_db=unimod_db[unimod_db['type'] != 'Manual']
    unimod_db=unimod_db[unimod_db['AA'] != 'X']
    #unimod_db=unimod_db[unimod_db['mass'] <= 120]
    drop = [False if ''.join(c for c in element if c.isdigit()==False) != element else True for element in unimod_db['PTM'].values]
    unimod_db['digit']=np.array(drop)
    unimod_db=unimod_db[unimod_db['digit'] == True]
    return unimod_db, unimod

def ptm_combos(df):
    column_values = ['peptide', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
    df_new = pd.DataFrame(columns = column_values)
    for row in df.values:
        if len(row[1])==0:
            array = np.array(row).reshape(1,-1)
            df_add = pd.DataFrame(data = array, columns = column_values, dtype='object')
            df_new = pd.concat([df_new, df_add], ignore_index = True)
        elif len(row[1])==1:
            array = np.array([[row[0], ptm, row[2], row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10]] for ptm in row[1][0][1]], dtype='object')
            df_add = pd.DataFrame(data = array, columns = column_values)
            df_new = pd.concat([df_new, df_add], ignore_index = True)
        elif len(row[1])==2:
            for ptm_N in row[1][0][1]:
                array = np.array([[row[0], ptm_N+'|'+ptm, row[2], row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10]] for ptm in row[1][1][1]], dtype='object')
                df_add = pd.DataFrame(data = array, columns = column_values)
                df_new = pd.concat([df_new, df_add], ignore_index = True)
    return df_new

def discover_ptm(AA_codes_tag, PTM, N_terms,C_terms, df):
    AA = initial_AA.copy()
    done = []
    keepptm=[]
    for x in df['PTM'].values:
        if len(x)==0:
            continue
        x = x.split('|')
        for i in range(0,len(x),2):
            keepptm.append(x[i+1])
    include = []
    for i in set(keepptm):
        count = keepptm.count(i)
        include.append((i,count))
    include = sorted(include, key=lambda x:x[1])[::-1]
    include = [num for num,val in include[:15]]
    for x, pep in df[['PTM','peptide']].values:
        if len(x)==0:
            continue
        x = x.split('|')
        for i in range(0,len(x),2):
            if x[i+1] not in include:
                continue
            if int(x[i]) == 0:
                N_terms[float(unimod[x[i+1]])]=x[i+1]
            elif int(x[i]) == -1:
                C_terms[float(unimod[x[i+1]])]=x[i+1]
            for loc, n in enumerate(unimod_db['AA'][unimod_db['PTM']==x[i+1]].values):
                if n != pep[int(x[i])-1]:
                    continue
                if x[i+1] in ['Methyl', 'Acetyl'] and int(x[i]) not in [0,-1]:
                    if n == pep[int(x[i])-1]:
                        AA[x[i+1]+str(loc)]=AA_codes_tag[n]+float(unimod[x[i+1]])
                        PTM[x[i+1]+str(loc)]=n
                elif x[i+1] in ['Oxidation'] and int(x[i]) not in [0,-1]:
                    if n == pep[int(x[i])-1] or n in ['M', 'P']:
                        AA[x[i+1]+str(loc)]=AA_codes_tag[n]+float(unimod[x[i+1]])
                        PTM[x[i+1]+str(loc)]=n
                elif x[i+1] in ['Deamidated'] and int(x[i]) not in [0,-1]:
                    if n not in ['Q', 'N']:
                        AA[x[i+1]+str(loc)]=AA_codes_tag[n]+float(unimod[x[i+1]])
                        PTM[x[i+1]+str(loc)]=n
                elif x[i+1] in ['Dioxidation'] and int(x[i]) not in [0,-1]:
                    if n == pep[int(x[i])-1]:
                        AA[x[i+1]+str(loc)]=AA_codes_tag[n]+float(unimod[x[i+1]])
                        PTM[x[i+1]+str(loc)]=n
                elif x[i+1] in done or n == 'I' or n not in AA:
                    continue
                elif x[i+1] not in ['Methyl', 'Acetyl']:
                    AA[x[i+1]+str(loc)]=AA_codes_tag[n]+float(unimod[x[i+1]])
                    PTM[x[i+1]+str(loc)]=n
            done.append(x[i+1])
    return AA, PTM, N_terms, C_terms
    
AA_codes = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "K" : 128.09496, "L" : 113.08406, "M" : 131.04049,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "R" : 156.10111, "S" : 87.03203,
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333}#, 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199}

AA_codes_tag = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "K" : 128.09496, "L" : 113.08406, "M" : 131.04049,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "R" : 156.10111, "S" : 87.03203,
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,"I" : 113.08406,'X':10000}#, 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199}

           
extra_AA = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "L" : 113.08406, "M" : 131.04049, "K" : 128.09496,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "S" : 87.03203,"R" : 156.10111, 
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,'Glu->pyro-Glu':129.04259-18.010565,'Gln->pyro-Glu':128.05858-17.026549}#, 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199}#,
            # 'Phospho':79.966331+87.03203,'Oxidation':131.04049+15.994915,'Oxidation2': 97.05276+15.994915,'Methyl':14.015650+99.06841,#,'Carbonyl1':13.979265+128.09496,
            # 'Decarboxylation1':115.02694-30.010565,'Decarboxylation2':129.04259-30.010565,'Carbamidomethyl1':57.021464+103.00919,'Carbonyl2':13.979265+156.10111,#'dihydroxy2':31.989829+131.04049,'dihydroxy1':31.989829+186.07931,'kynurenin':3.994915+186.07931,
            # , 'Phospho2':79.966331+163.06333, 'Phospho3':79.966331+101.04768,
            # 'Carbamyl3': 163.06333+43.005814, 'Carbamyl4': 87.03203+43.005814, 'Carbamyl2': 43.005814+101.04768,'Carbamyl': 156.10111+43.005814, 'Carbamyl5': 43.005814+128.09496}

           # 'Carbamidomethyl10':128.09496 + 8.014199 + 57.021464,'r':156.10111 + 10.008269,'k':128.09496 + 8.014199,
           # 'Carbamidomethyl2':57.021464+128.09496, 'Carbamidomethyl3':57.021464+137.05891, 'Carbamidomethyl4':57.021464+115.02694, 'Carbamidomethyl5':57.021464+129.04259,
           # 'Carbamidomethyl6':57.021464+87.03203, 'Carbamidomethyl7':57.021464+101.04768, 'Carbamidomethyl8':57.021464+163.06333, 'Carbamidomethyl9':57.021464+131.04049}            
            #for the histons
            #'Gluratylation':	114.031694+128.09496,
            #'di-Methylation':28.031300+156.10111,'tri-Methylation':42.046950+128.09496, 'benzoyl':104.026215+128.09496, 'Acetyl':42.010565+128.09496, 'Carbamyl': 156.10111+43.005814,
            # 'Butyryl':	70.041865+128.09496,  'Formyl':27.994915+128.09496, 'Crotonyl':128.09496+68.026215,'Amide':156.10111-0.984016,,'Gln->pyro-Glu':128.05858-18.010565,
            #'hydroxyisobutyryl': 86.036779+128.09496, 'Propionyl_light':56.026215+128.09496, 'Suc_anh_light':100.016044+128.09496,'carboxyethyl': 128.09496+72.021129, 'Methylthio': 103.00919+45.987721
            
            #For the teeth
            # 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199,
            # 'Carbamidomethyl1':57.021464+103.00919,'Carbamidomethyl10':128.09496 + 8.014199 + 57.021464,
            # 'Carbamidomethyl2':57.021464+128.09496, 'Carbamidomethyl3':57.021464+137.05891, 'Carbamidomethyl4':57.021464+115.02694, 'Carbamidomethyl5':57.021464+129.04259,
            # 'Carbamidomethyl6':57.021464+87.03203, 'Carbamidomethyl7':57.021464+101.04768, 'Carbamidomethyl8':57.021464+163.06333, 'Carbamidomethyl9':57.021464+131.04049
             
N_terms = {}#43.005814: 'Carbamyl',42.010565:'Acetyl', }#,57.021464:'Carbamidomethyl', ,56.026215:'Propionyl',58.005479:'Carboxymethyl',

C_terms = {}

ion_types = {"Y" : 19.0226, "B" : 1.0073} #z/a reeks etc. Pas helemaal op het einde als optie laten aanduiden

mods = {"-H2O" : ['S', 'T', 'E', 'D','Glu->pyro-Glu','Carbamyl2','Carbamyl4','Phospho','Carbamyl2','Phospho3','Carbamidomethyl4', 'Carbamidomethyl5', 'Carbamidomethyl6', 'Carbamidomethyl7', 'Decarboxylation2', 'Decarboxylation1'], "-NH3":['Formyl','Gluratylation','Suc_anh_light','hydroxyisobutyryl','tri-Methylation','Propionyl_light', 'Crotonyl', 'benzoyl', 'Acetyl', 'Butyryl','Amide', 'N', 'Q', 'K', 'R', 'k', 'r','di-Methylation', 'Carbamidomethyl10', 'Carbamidomethyl2', 'Gln->pyro-Glu', 'Lys->CamCys', 'Lys->MetOx', 'Carbamyl', 'carboxyethyl','Carbonyl1', 'Carbonyl2', 'Carbamyl5']} #nog andere opties? Opnieuw implementeren in functie

PTM = {'Gln->pyro-Glu':'Q', 'Glu->pyro-Glu':'E'}#'Met->AspSA': 'M', 'Homocysteic_acid':'M', 'Lys->MetOx':'K', 'Leu->MetOx':'L', 'Phe->CamCys':'F',
#        'Lys->CamCys':'K', 'Trp->Kynurenin':'W', 'Trp->Hydroxykynurenin':'W', 'Met->Hsl': 'M', 'Met->Hse':'M',
#        'Gln->pyro-Glu':'Q', 'Carbamyl':'R', 'Asn->Asp':'N', 'Gln->Glu':'Q', 'Methylthio':'C', 'Oxidation':'M',
#        'Acetyl':'K', 'Butyryl':'K', 'Amide':'R', 'di-Methylation':'R', 'Formyl':'K', 'Phospho':'S', 'tri-Methylation':'K',
#        'benzoyl':'K', 'Crotonyl':'K', 'Gluratylation':'K','carboxyethyl':'K', 'Suc_anh_light':'K', 'hydroxyisobutyryl':'K',
#        'Ser->Ala':'S', 'Propionyl_light':'K', 'Oxidation2':'P', 'Carbamyl2':'T', 'Carbamidomethyl1':'C', 'Carbamidomethyl2':'K',
#        'Carbamidomethyl3':'H', 'Carbamidomethyl4':'D', 'Carbamidomethyl5':'E', 'Carbamidomethyl6':'S', 'Carbamidomethyl7':'T','Hydroxylation':'P',
#        'Carbamidomethyl8':'Y', 'Carbamidomethyl9':'M', 'Carbamidomethyl10':'k', 'Phospho2':'Y', 'Phospho3':'T', 'kynurenin':'W', 'dihydroxy1':'W',
#        'dihydroxy2':'M', 'Carbonyl1':'K', 'Carbonyl2':'R', 'Carbamyl3':'Y', 'Carbamyl4':'S', 'Glu->pyro-Glu':'E', 'Decarboxylation2':'E', 'Decarboxylation1':'D', 'Methyl':'V', 'Carbamyl5':'K'}
#option of heavy aminoacids, and option of acetylation, carbamylation at start AA. 

initial_AA = extra_AA.copy()
initial_N = {}
initial_C = {}
#####################################################################


#prosit data voor controle  PXD009565
################### load in data ####################################

if __name__ == '__main__':
    path = 'C:/Users/Gebruiker/Desktop/neolitic_protein_discovery'
    os.chdir(path)
    MGF_zip = ZipFile('./MGF_files.zip')
    make_mgf = []
    df = [mgf_file.filename
           for mgf_file in MGF_zip.infolist()
           if mgf_file.filename.endswith('.mgf')or mgf_file.filename.endswith('.txt')]
    personal_db = personal_input(path)
    crap, crap_individual = crap_f()
    unimod_db, unimod = make_unimod(include_label = False)
    heavy = True
    for files in df:
        print(files)
        if 'BASL' not in files:
            continue
        end_result = {}
        nr_file =0
        with mgf.read('.\\'+files.replace('/', '\\')) as reader:
            AA_code = AA_codes_tag
            df_deeplc = do_tag(reader)
            extra_AA, PTM, N_terms, C_terms = discover_ptm(AA_codes_tag,PTM, initial_N,initial_C, df_deeplc)
            #make plastic flag
        AA_code = AA_codes
        temp_Bnormal = make_list_pep(420, 'B', AA_codes)
        print('lengt B start', len(temp_Bnormal))
        temp_Ynormal = make_list_pep(420, 'Y', AA_codes)
        print('length Y start', len(temp_Ynormal))
        AA_code = extra_AA
        temp_Bmodifications = make_list_pep(420, 'B', extra_AA)
        print('lengt B start', len(temp_Bmodifications))
        temp_Ymodifications = make_list_pep(420, 'Y', extra_AA)
        print('length Y start', len(temp_Ymodifications))
        already = [int(num.split('_')[0]) for num in df_deeplc['spectrum_nr'].values]
        df_deeplc['spectrum_nr']=np.array(already)#renames the split tag search
        if len(df_deeplc)>2:
            tag_out = max(df_deeplc['rt_time'].values)+min(df_deeplc['rt_time'].values)*0.05
        else:
            tag_out=0
        with mgf.read('.\\'+files.replace('/', '\\')) as reader:
            for spectrum in reader:
                nr_file += 1
                if nr_file != 152:
                    continue
                # if nr_file in already or (float(spectrum['params']['rtinseconds']) >tag_out and tag_out != 0):
                #     continue
                spectral_match = True
                if 'charge' in spectrum['params']:
                    charging = spectrum['params']['charge']
                    spectral_match = False
                else:
                    charging  = [2,3]
                for test_charge in charging:
                    list_1 = []
                    list_2 = []
                    if test_charge == 3 and spectral_match == True and list(spectrum['params']['pepmass'])[0] > 380:
                        continue
                    spectrum['params']['charge']=[test_charge]

                    title = 'Spectrum_' + str(nr_file)+'_'+str(test_charge)
                    print(title)
                    xpoints_data, ypoints_data, pp_mass, dev, to_mgf, nx, ny = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'], spectrum)
                    if len(set([int(xs) for xs in xpoints_data]))<=pp_mass*test_charge/100 or len(xpoints_data)>234 or pp_mass*test_charge>2200:#or min(xpoints_data)>=400
                        print('wrong peaks')
                        #end_result[title+'_contamination'+spectrum['params']['title']] = 'Indeterminable'
                        spectral_match = False
                        continue
                    for extras in ['normal', 'modifications']:
                        if extras == 'normal':
                            AA_code = AA_codes
                        else:
                            if len(AA_code)==len(extra_AA):
                                break
                            AA_code = extra_AA
                        #locals()[str(nr_file)+'_to_mgf'] = to_mgf
                        #make_mgf.append(locals()[str(nr_file)+'_to_mgf'])
                        print('Start analysis of: ', spectrum['params']['title'])
                        rt_time = float(spectrum['params']['rtinseconds'])
                        if pp_mass>400 or spectrum['params']['charge'][0] > 2:
                            MZ_start = 420 #kan aangepast worden als meer PTMs vergeet niet bij iteration te kijken
                        else:
                            MZ_start = pp_mass
                        for ion in ['Y','B']:
                            if ion == 'B' and len(Y_list)==0:
                                print('no result possible')
                                B_list = []
                                break
                            print('begin ' +ion+'-peptide'+' '+extras)
                            if pp_mass*spectrum['params']['charge'][0] <= 1000:
                                max_window = (pp_mass*spectrum['params']['charge'][0])*0.7#((-pp_mass*spectrum['params']['charge'][0]*0.0001)+0.85-(spectrum['params']['charge'][0]-2)*0.025) #0.85kan verder gaan? Maar relevant??
                            else:
                                max_window = (pp_mass*spectrum['params']['charge'][0])*0.65#((0.0001*pp_mass*spectrum['params']['charge'][0])+0.6-(spectrum['params']['charge'][0]-2)*0.05) #0.35b-reek hierna niet meer betrouwbaar genoeg 
                            MZ_iteration = MZ_start
                            nr_iterations = math.ceil(abs(max_window-MZ_start)/75)
                            
                            
                            
                            locals()[ion+'_list'] = find_pep_ions_part1(xpoints_data, ypoints_data, MZ_start, ion, locals()['temp_'+ion+extras])
                            
                            while nr_iterations != 0 and max_window>MZ_start:
                                print('Still', nr_iterations, 'in-silico peptides iteration to go')
                                if nr_iterations != 1:
                                    MZ_iteration += abs(max_window - MZ_iteration)/((nr_iterations-1)*1.25)
                                else:
                                    MZ_iteration += abs(max_window - MZ_iteration)
                                locals()[ion + '_list']= find_pep_ions_iterations(locals()[ion+'_list'], MZ_iteration, ion)
                                nr_iterations -= 1
                                # for a in locals()[ion + '_list']:
                                #     if a[0].upper() in 'A|V|G|A|N|P|E|Q|L|T|R' or a[0].upper() in 'A|V|G|A|N|P|E|Q|L|T|R'[::-1]:
                                #         print(a)
                            if ion == 'Y':
                                cutoff = find_treshold(locals()[ion+'_list'], 750)
                                locals()[ion+'_list'] =[num for num in locals()[ion+'_list'] if num[-1]>=cutoff]
                            else: 
                                cutoff = find_treshold(locals()[ion+'_list'], 750)
                                locals()[ion+'_list'] =[num for num in locals()[ion+'_list'] if num[-1]>=cutoff]
                            print('end of ' + ion + '-peptide')
                        if extras == 'normal':
                            for nr, ion in enumerate(['Y', 'B']):
                                locals()['list_'+str(nr+1)] = locals()[ion+'_list']
                        else:
                            for nr, ion in enumerate(['Y', 'B']):
                                temp_list = locals()['list_'+str(nr+1)]
                                locals()['list_'+str(nr+1)] = locals()['list_'+str(nr+1)]+[num for num in locals()[ion+'_list'] if num not in temp_list]

                    overlap_output = find_overlap_sequence(list_1, list_2)
                    
                    combined_peptide_list = score_overlap(overlap_output, nx, ny)
                    
                    combined_peptide_list = deamidation(combined_peptide_list)
                    
                    end_result, found = annotate(combined_peptide_list, end_result)
                    if len(combined_peptide_list) == 0:
                        spectral_match = False
                        continue
                    # else:
                    #     if title not in end_result.keys() and found == False:
                    #         end_result[title+'_unkown'] = len(combined_peptide_list)
                    
                    #here we put all information needed to perform deeplc in one run
                    combined_peptide_list = IL_adjusted(combined_peptide_list)
                    mz, intens = nx, ny#xpoints_data, ypoints_data#pip_data(spectrum['m/z array'], spectrum['intensity array'], spectrum['charge array'], pp_mass)
                    if 'charge' in spectrum['params'].keys():
                        charges = spectrum['params']['charge'][0]
                    else:
                        charges = 2
                    array = np.array([[pep[0].upper(), pep[1], float(pep[2]), float(pep[3]), float(pep[4]), rt_time,mz, intens, nr_file, charges, spectrum['params']['title']] for pep in combined_peptide_list], dtype='object')
                    column_values = ['peptide', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
                    df_add = pd.DataFrame(data = array, columns = column_values)
                    df_deeplc = pd.concat([df_deeplc, df_add], ignore_index = True)
        #here do the deeplc on all spectra in file
        pip_score, s_calibration = filters('deeplc_file', df_deeplc) 
        #here do the final output ranking

        rank = rank_results(pip_score)
        print('save results of ', files)
        already = []
        with open('Results_thesis_'+files.split('/')[1]+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['Sequences used for calibration'])
            writer.writerow(['nr','seq', 'file_name', 'PTM', 'ID', 'covergae', 'ms2pip_score'])
            for spec, rest in s_calibration.items():
                for i in rest:
                    writer.writerow([i[6],i[0], i[4], i[1], i[2], i[5], i[3]])
            writer.writerow(['Sequences from personal database'])
            writer.writerow(['nr','seq', 'spectrum', 'PTM', 'ID', 'coverage', 'ms2pip_score'])
            rank['peptides']=rank.index
            for seq, spec, ptm, ID, cov, ms2, nr in rank[['peptides', 'file_name', 'PTM', 'ID', 'coverage', 'ms2pip_score', 'spectrum_nr']].values:
                if ID != 'unkown':
                    already.append(nr)
                    writer.writerow([nr, seq, spec, ptm, ID, cov, ms2])
            writer.writerow(['unkown sequences'])
            writer.writerow(['nr','seq', 'spectrum', 'PTM', 'ID', 'coverage', 'ms2pip_score'])
            for seq, spec, ptm, ID, cov, ms2, nr in rank[['peptides', 'file_name', 'PTM', 'ID', 'coverage', 'ms2pip_score', 'spectrum_nr']].values:
                if ID == 'unkown' and nr not in already:
                    writer.writerow([nr, seq, spec, ptm, ID, cov, ms2])
    #consensus = find_consensus(rank)
    #w_mgf(make_mgf)
    print('End of analysis')


