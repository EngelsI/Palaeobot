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
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings('ignore')
#####################################################################

################### Functions #######################################
def insilico_peptide(pep, new_pep, temp, MZ):
    new_pep_list = {}
    for AA in pep:
        for AA2 in new_pep:
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
    for num in spectrums:
        for value in x_array:

            if value-dev <= num <= value+dev:
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
        if '*' in i:
            element = i.split('*')
            if element[1][0]=='|':
                new.append((element[1][1:], float(element[0])))
            else: 
                new.append((element[1], float(element[0])))
        else:
            new.append((i,0))
    return new
        
def find_pep_ions_part1(xpoints_data, ypoints_data, MZ, ion):
    temp =set()
    pep = list(AA_code.keys())
    MZ_range = xpoints_data#spectrum['m/z array']
    ind = np.argsort(MZ_range)
    MZ_range = MZ_range[ind]
    
    if ion=='B':
        new_pep = {}
        for i,u in AA_code.items():
            new_pep[i]=u
        for item in N_terms.keys():
            new_pep[str(item)+'*']=item
    elif ion=='Y':
        new_pep = AA_code
    print('making the in-silico peptides')
    while len(new_pep)!=0:
        new_pep, temp = insilico_peptide(pep, new_pep, temp, MZ)
    temp = catch_asterix(temp)
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
    return overlap_output_final

def score_overlap(overlap_output, xpoints_data, ypoints_data):
    
    print('final pearson correlation')
    
    PC_list =[]
    for item in overlap_output:
        a, c = make_overlap_spectrum(item[2][0][::-1], item[2][0], xpoints_data, ypoints_data, item[2][1], item[2][2])
        j = jaccard(ypoints_data, c, item[2][0])
        pc = pearson_correlation(ypoints_data, c)[0][1]
        found_intensity = np.sum(c)
        all_intensity = np.sum(ypoints_data)
        add = tuple([item, j, found_intensity/all_intensity, pc]) #j*found_intensity/all_intensity
        PC_list.append(add)
        # if 'GAYVEVTAk' in item[0]:
        #     print(make_overlap_spectrum(item[2][::-1], item[2], xpoints_data, ypoints_data))
        #     print(jaccard(ypoints_data, c, item[0]))
        #     print(pearson_correlation(ypoints_data, c)[0][1])
        #     print(np.sum(c))
    
    PC_list = sorted(PC_list, key=lambda x: x[3])[::-1]

    if len(PC_list)==0:
    #     make_final_plots(xpoints_data, ypoints_data, PC_list[0][0][2][::-1], PC_list[0][0][2], PC_list)
    # else:
        return []
    #cutoff = find_treshold(PC_list)
    #PC_list = [num for num in PC_list if PC_list.count(num)>1]
    for a in PC_list:
        if 'AVGANPEQLTr' in a[0]:
            print(a)
    ##################DIT GAAT NOG AANGEPAST MOETEN WORDEN########################"""""
    # PC_list = [num for num in PC_list if num[3]>0]
    # thresh = [(num[0], num[3]) for num in PC_list]
    # threshold = find_treshold_overlap(thresh)
    # PC_list =[num for num in PC_list if num[3]>=threshold]
    begin = len(PC_list)
    PC_list = find_treshold_overlap(PC_list)
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
    crap = {}
    for record in SeqIO.parse("./crap.fasta.txt", "fasta"):
        crap[record.seq]=record.id
    for record in SeqIO.parse("./contaminants.fasta", "fasta"):
        crap[record.seq]=record.id
    lysyl = 'MHKRTYLNACLVLALAAGASQALAAPGASEMAGDVAVLQASPASTGHARFANPNAAISAAGIHFAAPPARRVARAAPLAPKPGTPLQVGVGLKTATPEIDLTTLEWIDTPDGRHTARFPISAAGAASLRAAIRLETHSGSLPDDVLLHFAGAGKEIFEASGKDLSVNRPYWSPVIEGDTLTVELVLPANLQPGDLRLSVPQVSYFADSLYKAGYRDGFGASGSCEVDAVCATQSGTRAYDNATAAVAKMVFTSSADGGSYICTGTLLNNGNSPKRQLFWSAAHCIEDQATAATLQTIWFYNTTQCYGDASTINQSVTVLTGGANILHRDAKRDTLLLELKRTPPAGVFYQGWSATPIANGSLGHDIHHPRGDAKKYSQGNVSAVGVTYDGHTALTRVDWPSAVVEGGSSGSGLLTVAGDGSYQLRGGLYGGPSYCGAPTSQRNDYFSDFSGVYSQISRYFAP'
    crap[lysyl]='lysyl'
    keratin1 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSEGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1] = 'keratin1'
    keratin1_2 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSGGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1_2] = 'keratin1_2'
    ##add below if no contaminants in files
    # for i, o in personal_db.items():
    #     crap[i]=o
    return crap

def personal_input(path):
    files_own = []
    for fastafile in os.walk(path+'/fasta_db'):
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta'):
                files_own.append(path+'/fasta_db/'+i)
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
        for record, ptm, rt, x, y, cs, nr, fn, cov, pc_sc in df[['peptide', 'PTM', 'rt_time', 'xpoints_data', 'ypoints_data', 'charge_state', 'spectrum_nr', 'file_name', 'coverage', 'pc_score']].values:
            for sequence in crap.keys():
                if record.upper() in sequence:
                    mz, intensity, annotation = ms2pip_sp.predict(record.upper(), '-', cs)
                    intensity = [num if num >0 else 1e-15 for num in intensity]
                    to_test, new_intensity = find_intensity(record,  x, y, intensity, ptm)
                    if len(to_test)==len(new_intensity): 
                        pc = pearson_correlation(to_test, new_intensity)[0][1]
                        if pc < 0.45 and (pc_sc<0.98 or len(record)<6 or pc <0):
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
    #Y = np.array([1 if i!='unkown' else 0 for i in df[['ID']].values])
    mini = min(df['rt_time'].values)*0.05
    Y2 = np.array([1 if i<=mini else 0 for i in df[['difference']].values])

    plt.title('plot of all predictions')
    #plt.scatter([num[0] for i, num in enumerate(X) if Y[i]==1], [num[1] for i, num in enumerate(X) if Y[i]==1], c='r')
    plt.scatter([num[0] for i, num in enumerate(X)], [num[1] for i, num in enumerate(X)], c='b', alpha=0.01)
    plt.scatter([num[0] for i, num in enumerate(X) if Y2[i]==1], [num[1] for i, num in enumerate(X) if Y2[i]==1], c='g', alpha=0.1)
    plt.show()
    df['error_interval'] = np.array([mini]*len(df))
    df = df[abs(df['difference'])<=mini] 
    df = df.set_index('peptide', drop=False)
    
    print('save results for ms2pip')
    nr = 0
    with open('peprec.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['spec_id','modifications','peptide','charge'])
        for seq, ch, ptm in df[['peptide', 'charge_state', 'PTM']].values:
            if ptm == '':
                ptm = '-'
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
    # "Carbamidomethyl,57.021464,opt,K",
    # "Carbamidomethyl,57.021464,opt,H",
    # "Carbamidomethyl,57.021464,opt,D",
    # "Carbamidomethyl,57.021464,opt,E",
    # "Carbamidomethyl,57.021464,opt,S",
    # "Carbamidomethyl,57.021464,opt,T",
    # "Carbamidomethyl,57.021464,opt,P",
    # "Carbamidomethyl,57.021464,opt,Y",
    # "Carbamidomethyl,57.021464,opt,M",
    # "Carbamidomethyl,57.021464,opt,C",
    # "Acetyl,42.010565,opt,N-term",
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
    return df, spectrum_calibration
         
def write_csv(title, combined_peptide_list):
    print('save in csv file')
    with open('records_'+title+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['seq', 'modifications'])
        for record, ptm in combined_peptide_list:
            writer.writerow([record, ptm])
    return 0
     
def set_data(xpoints_data, ypoints_data, pp_mass, charge):
          
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
    
    dev = 0.01+10*(spectrum['params']['charge'][0]*pp_mass)/1e6
    if dev <= 0.01:
          dev = dev*2
    if dev >0.0185:
        dev = 0.0185
    if len(xpoints_data)> 1:
        plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
        plt.title(title+ ' before mirror')
        plt.show()
    if len(xpoints_data)<=(spectrum['params']['charge'][0]*pp_mass)/15+25:
        print('only ghost performed')
        neutral_lossX, neutral_lossY, mdev = mirror_B(xpoints_data, ypoints_data, pp_mass, 1)
    else:
        print('mirroring')
        neutral_lossX, neutral_lossY, mdev = mirror_B(xpoints_data, ypoints_data, pp_mass)
    
    if (spectrum['params']['charge'][0]*pp_mass)/20+30 >= len(xpoints_data):
        print('only ghost performed 2')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass, 1)
    else:
        print('mirroring')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass)

    neutral_lossY = tic_normalize(neutral_lossY)
    if len(neutral_lossY)>0:
        ind = neutral_lossY >1e-3
        neutral_lossX = neutral_lossX[ind]
        neutral_lossY = neutral_lossY[ind]
    
    ypoints_data = tic_normalize(ypoints_data)
    if len(ypoints_data)>0:
        ind = ypoints_data >1e-3
        xpoints_data = xpoints_data[ind]
        ypoints_data = ypoints_data[ind]
    if len(xpoints_data)>1:
        plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
        
        plt.title(title+ ' after')
        plt.show()
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
        
def find_formula(file_name, a, predictions):
    title = file_name  
    preds = predictions#deeplc(title)
    preds = np.array(preds)
    #a = np.array([1190.295, 1287.507, 1291.205, 1452.976, 1531.312, 1605.268, 1580.191, 1623.252, 1654.071, 1701.912, 1783.793, 1786.38, 1992.298, 2275.02])
    
    i = np.argsort(a)
    new_a = a[i]
    new_preds = preds[i]
    log2_a = np.array([-math.log(x,2) for x in new_a])
    log2_preds = np.array([math.log(x,2) for x in new_preds])
    
    x_p = [0]
    y_p = [0]
    for i in range(0,len(new_a)-1):
        A = (log2_a[i], new_a[i])
        B = (log2_preds[i], new_preds[i])
        C = (log2_a[i+1], new_a[i+1])
        D = (log2_preds[i+1], new_preds[i+1])
        x, y = line_intersection((A, B), (C, D))
        x_p.append(x)
        y_p.append(y)
    
        
    # plt.scatter(log2_a, new_a, c='b')
    # plt.scatter(log2_preds, new_preds, c='r')
    # plt.show()
    
    # plt.scatter(log2_a, new_a, c='b')
    # plt.scatter(log2_preds, new_preds, c='r')
    # plt.scatter(x_p, y_p, c='g')
    # plt.show()
    
    column_values = ['real']
    df = pd.DataFrame(data = new_a, columns = column_values)
    df['log2_real'] = log2_a
    df['predicted'] = new_preds
    df['log2_predicted']=log2_preds
    df['x_interval'] = x_p
    df['y_interval'] = y_p
    
    return df

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
    
    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]
    
    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')
    
    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

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
        ranking = sorted(ranking)[::-1][:5] #top 5 so chimear also incuded?
        ranking = ranking[-1]
        peps = pip_results['peptide'][(pip_results['spectrum_nr']==spectra) & ((pip_results['rank_score'] >= ranking)&(pip_results['rank_score'] != 0))].values
        
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
        for sequence in crap.keys():
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

def mirror_A(x, y, peptide_mass, mdev = 0, takex=[], takey=[]):
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
        dev = 0.01+10*peptide_mass/1e6
        if dev <= 0.04:
              dev = dev*2
        if dev >0.06:
            dev = 0.08
    else:
        dev = mdev*0.95
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[])#1.0073,19.0226
        mirror_x = sorted([peptide_mass - u for u in x_withh2o])
        
        new_x =[]
        for element in x:
            added = False
            for value in mirror_x:
                if element-dev<=value <= element+dev:
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
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew) #peaks that were not accounted for can be double charged

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
            if xnew[i+1]-xnew[i] >0.01 and xnew[i] not in already:
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
    if len(finalx)>max(peptide_mass/15,200):
        finalx, finaly, mdev = mirror_A(finalx, finaly,peptide_mass/spectrum['params']['charge'][0], dev, takex, takey)
    elif 10 <= len(finalx)<=max(peptide_mass/15,200):
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
        dev = 0.01+10*peptide_mass/1e6
        if dev <= 0.03:
              dev = dev*2
        if dev >0.06:
            dev = 0.05
    else:
        dev = mdev*0.95
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[])#1.0073,19.0226
        mirror_x = sorted([peptide_mass - u for u in x_withh2o])
        
        new_x =[]
        for element in x:
            added = False
            for value in mirror_x:
                if element-dev<=value <= element+dev:
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
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew) #peaks that were not accounted for can be double charged

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
    if len(finalx)>peptide_mass/15+25:
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
            if xnew[i+1]-xnew[i] >0.01 and xnew[i] not in already:
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

def ghost_peaks(x, y, peptide_mass, xnew, ynew):
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

def adjust_intensity(x, y, d):
    x = list(x)
    y = list(y)
    newy = []
    for key, val in d.items():
        newy.append((y[x.index(val)]+y[x.index(key)]))
    return np.array(newy)

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
                    new_ptm = str(x+1)+'|Deamidation'
                elif len(ptm.split('|'))==2 or x+1<= int(ptm.split('|')[0]) or x+1>= int(ptm.split('|')[-2]):
                    if x+1 >= int(ptm.split('|')[-2]):
                        new_ptm = ptm+'|'+str(x+1)+'|Deamidation'
                    else:
                        new_ptm = str(x+1)+'|Deamidation|'+ptm
                else:
                    new_ptm = ptm.split('|')
                    for n in range(0,len(new_ptm),2):
                        if int(new_ptm[n]) <= x+1 and int(new_ptm[n+2]) >= x+1:
                            new_ptm.insert(n+2, str(x+1)+'|Deamidation')
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
                    new_ptm = str(x+1)+'|Deamidation'
                elif len(ptm.split('|'))==2 or x+1<= int(ptm.split('|')[0]) or x+1>= int(ptm.split('|')[-2]):
                    if x+1 >= int(ptm.split('|')[0]):
                        new_ptm = ptm+'|'+str(x+1)+'|Deamidation'
                    else:
                        new_ptm = str(x+1)+'|Deamidation|'+ptm
                else:
                    new_ptm = ptm.split('|')
                    for n in range(0,len(new_ptm),2):
                        if int(new_ptm[n]) <= x+1 and int(new_ptm[n+2]) >= x+1:
                            new_ptm.insert(n+2, str(x+1)+'|Deamidation')
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

AA_codes = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "K" : 128.09496, "L" : 113.08406, "M" : 131.04049,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "R" : 156.10111, "S" : 87.03203,
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333}#, 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199}
           
extra_AA = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "L" : 113.08406, "M" : 131.04049, "K" : 128.09496,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "S" : 87.03203,"R" : 156.10111, 
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,#, 'r':156.10111 + 10.008269,'k':128.09496 + 8.014199}#,
            'Phospho':79.966331+87.03203,'Oxidation':131.04049+15.994915,'Oxidation2': 97.05276+15.994915,'Methyl':14.015650+99.06841,#,'Carbonyl1':13.979265+128.09496,
            'Decarboxylation1':115.02694-30.010565,'Decarboxylation2':129.04259-30.010565,'Carbamidomethyl1':57.021464+103.00919,'Carbonyl2':13.979265+156.10111,#'dihydroxy2':31.989829+131.04049,'dihydroxy1':31.989829+186.07931,'kynurenin':3.994915+186.07931,
            'Glu->pyro-Glu':129.04259-18.010565,'Gln->pyro-Glu':128.05858-17.026549, 'Phospho2':79.966331+163.06333, 'Phospho3':79.966331+101.04768,
            'Carbamyl3': 163.06333+43.005814, 'Carbamyl4': 87.03203+43.005814, 'Carbamyl2': 43.005814+101.04768,'Carbamyl': 156.10111+43.005814, 'Carbamyl5': 43.005814+128.09496}

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
             
N_terms = {43.005814: 'Carbamyl'}#42.010565:'Acetyl', }#,57.021464:'Carbamidomethyl', ,56.026215:'Propionyl',58.005479:'Carboxymethyl',

ion_types = {"Y" : 19.0226, "B" : 1.0073} #z/a reeks etc. Pas helemaal op het einde als optie laten aanduiden

mods = {"-H2O" : ['S', 'T', 'E', 'D','Glu->pyro-Glu','Carbamyl2','Carbamyl4','Phospho','Carbamyl2','Phospho3','Carbamidomethyl4', 'Carbamidomethyl5', 'Carbamidomethyl6', 'Carbamidomethyl7', 'Decarboxylation2', 'Decarboxylation1'], "-NH3":['Formyl','Gluratylation','Suc_anh_light','hydroxyisobutyryl','tri-Methylation','Propionyl_light', 'Crotonyl', 'benzoyl', 'Acetyl', 'Butyryl','Amide', 'N', 'Q', 'K', 'R', 'k', 'r','di-Methylation', 'Carbamidomethyl10', 'Carbamidomethyl2', 'Gln->pyro-Glu', 'Lys->CamCys', 'Lys->MetOx', 'Carbamyl', 'carboxyethyl','Carbonyl1', 'Carbonyl2', 'Carbamyl5']} #nog andere opties? Opnieuw implementeren in functie

PTM = {'Met->AspSA': 'M', 'Homocysteic_acid':'M', 'Lys->MetOx':'K', 'Leu->MetOx':'L', 'Phe->CamCys':'F',
       'Lys->CamCys':'K', 'Trp->Kynurenin':'W', 'Trp->Hydroxykynurenin':'W', 'Met->Hsl': 'M', 'Met->Hse':'M',
       'Gln->pyro-Glu':'Q', 'Carbamyl':'R', 'Asn->Asp':'N', 'Gln->Glu':'Q', 'Methylthio':'C', 'Oxidation':'M',
       'Acetyl':'K', 'Butyryl':'K', 'Amide':'R', 'di-Methylation':'R', 'Formyl':'K', 'Phospho':'S', 'tri-Methylation':'K',
       'benzoyl':'K', 'Crotonyl':'K', 'Gluratylation':'K','carboxyethyl':'K', 'Suc_anh_light':'K', 'hydroxyisobutyryl':'K',
       'Ser->Ala':'S', 'Propionyl_light':'K', 'Oxidation2':'P', 'Carbamyl2':'T', 'Carbamidomethyl1':'C', 'Carbamidomethyl2':'K',
       'Carbamidomethyl3':'H', 'Carbamidomethyl4':'D', 'Carbamidomethyl5':'E', 'Carbamidomethyl6':'S', 'Carbamidomethyl7':'T','Hydroxylation':'P',
       'Carbamidomethyl8':'Y', 'Carbamidomethyl9':'M', 'Carbamidomethyl10':'k', 'Phospho2':'Y', 'Phospho3':'T', 'kynurenin':'W', 'dihydroxy1':'W',
       'dihydroxy2':'M', 'Carbonyl1':'K', 'Carbonyl2':'R', 'Carbamyl3':'Y', 'Carbamyl4':'S', 'Glu->pyro-Glu':'E', 'Decarboxylation2':'E', 'Decarboxylation1':'D', 'Methyl':'V', 'Carbamyl5':'K'}
#option of heavy aminoacids, and option of acetylation, carbamylation at start AA. 
#####################################################################


#prosit data voor controle
################### load in data ####################################

if __name__ == '__main__':
    MGF_zip = ZipFile('./MGF_files.zip') #teeth_thesis, MGF_files
    make_mgf = []
    df = [mgf_file.filename
           for mgf_file in MGF_zip.infolist()
           if mgf_file.filename.endswith('.mgf')]
    path = 'C:/Users/Gebruiker/Desktop/neolitic_protein_discovery'
    os.chdir(path)
    personal_db = personal_input(path)
    crap = crap_f()
    heavy = True
    other_retention_time = False
    for files in df:
        print(files)
        
        end_result = {}
        nr_file =0
        column_values = ['peptide', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
        df_deeplc = pd.DataFrame(columns = column_values)
        
        with mgf.read('.\\'+files.replace('/', '\\')) as reader:
            for spectrum in reader:
                nr_file += 1
                
                spectral_match = True
                if 'charge' in spectrum['params']:
                    charging = spectrum['params']['charge']
                    spectral_match = False
                else:
                    charging  = [2,3]
                for test_charge in charging:
                    list_1 = []
                    list_2 = []
                    if test_charge == 3 and spectral_match == True and list(spectrum['params']['pepmass'])[0] > 420:
                        continue
                    spectrum['params']['charge']=[test_charge]

                    title = 'Spectrum_' + str(nr_file)+'_'+str(test_charge)
                    print(title)
                    xpoints_data, ypoints_data, pp_mass, dev, to_mgf, nx, ny = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'])
                    if len(xpoints_data)<=pp_mass*test_charge/200 or len(xpoints_data)>234 or pp_mass*test_charge>2100:#or min(xpoints_data)>=400
                        print('wrong peaks')
                        #end_result[title+'_contamination'+spectrum['params']['title']] = 'Indeterminable'
                        spectral_match = False
                        continue
                    for extras in ['normal', 'modifications']:
                        if extras == 'normal':
                            AA_code = AA_codes
                        else:
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
                            locals()[ion+'_list'] = find_pep_ions_part1(xpoints_data, ypoints_data, MZ_start, ion)
                            
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
                    
                    #end_result, found = annotate(combined_peptide_list, end_result)
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


