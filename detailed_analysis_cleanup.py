# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 16:12:30 2023

@author: Ian
"""
from pyteomics import mgf
from pyteomics import mzid, auxiliary, pepxml, mzml, mgf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
ion_types = {"Y" : 19.0226, "B" : 1.0073}

AA_code = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "K" : 128.09496, "L" : 113.08406, "M" : 131.04049,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "R" : 156.10111, "S" : 87.03203,
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,"I" : 113.08406,'X':1000}

def set_data(xpoints_data, ypoints_data, pp_mass, charge,parameter1,parameter2):
    ind = np.argsort(xpoints_data)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    if 'charge' in spectrum['params'].keys():
        if spectrum['params']['charge'][0] > 2:
            pp_mass = pp_mass-((spectrum['params']['charge'][0]-2)**2)/spectrum['params']['charge'][0]#((pp_mass*3)/2)-0.5
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
    
    dev = parameter1*(spectrum['params']['charge'][0]*pp_mass)/1e6 #20ppm

    xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass,parameter1,parameter2)
    return xpoints_data, ypoints_data, pp_mass, dev

def tic_normalize(intensity):
        
        intensity = np.array(intensity)
        return intensity / np.sum(intensity)

def mirror_A(x,y,peptide_mass,parameter1,parameter2, mdev = 0, takex=[], takey=[]):
    peptide_mass = peptide_mass*spectrum['params']['charge'][0]
    if mdev == 0:
        takex=x
        takey=y
        dev = parameter1*peptide_mass/1e6 #20ppm
    else:
        dev = mdev*0.95
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
                if value-dev<=element <= value+dev:# and ((y[loc]+y[locx])>sorted(y)[int(len(y)*0.6)] or mdev not in [0,1]):
                    new_x.append(True)
                    added = True
                    for XXX in [18.010565, 17.026549]:
                        for element2 in x:
                            if element -XXX-0.05 <= element2 <= element-XXX+0.05:
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
    # if len(finalx)>max(len(x),len(takex))*parameter2:
    #     finalx, finaly, mdev = mirror_A(np.array(finalx), np.array(finaly),peptide_mass/spectrum['params']['charge'][0],parameter1,parameter2, dev, takex, takey)
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
            if xnew[i+1]-xnew[i] >0 and xnew[i] not in already:#0.01
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
    ghost = [num*2-1 for i, num in enumerate(x) if num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew]
    xghost = [num for i, num in enumerate(x) if num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew]
    ghost_intensity = [num for i, num in enumerate(y) if x[i]<=peptide_mass*(spectrum['params']['charge'][0]/2) and x[i] not in xnew]
    ghost = list(np.array(ghost))
    x = xghost+ghost
    
    y=ghost_intensity+ghost_intensity
    x= np.array(x)
    y=np.array(y)
    
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    return x, y


# df = mzid.DataFrame('.\F001266.mzid')
# df2 = mzid.DataFrame('.\F001266 (1).mzid')
if __name__ == '__main__':
    test = {}
    ionsID = {}
    for i in mzid.read(".\F010801.mzid"):
        mzs = i['SpectrumIdentificationItem'][0]['IonType']
        if i['SpectrumIdentificationItem'][0]['chargeState'] not in [2,3]:
            continue
        mz = []
        ions = []
        for n in mzs:
            # if 'frag: b ion' in n.keys() or 'frag: y ion' in n.keys():
            for io in list(n.keys()):
                ions.append(str(io))
            if n['charge'] ==2:
                for element in list(n['FragmentArray'][0]['values']):
                    mz.append(element*2-1) 
            else:
                mz = mz + list(n['FragmentArray'][0]['values'])
        test[i['spectrum title']]=mz
        ionsID[i['spectrum title']]=ions
        
    for parameter1 in [1,5,10,20,50]:
        nr_file =0
        result = {}
        x=[]
        y=[]
        intensity = []
        mascotid = []
        lijstje =[]
        y_ions = 0
        b_ions = 0
        both_ions = 0
        only_neutral = 0
        with mgf.read('.\mascot_daemon_merge.mgf') as reader:
            for spectrum in reader:
                nr_file += 1
                if len(spectrum['m/z array'])==0:
                    continue
                if spectrum['params']['title'] not in test.keys():
                    continue
                if 'charge' in spectrum['params']:
                    charging = spectrum['params']['charge']
                else:
                    charging  = [2,3]
                for test_charge in charging:
                    
                    parameter2 = 0.25
                    spectrum['params']['charge']=[test_charge]
                    title = 'Spectrum_' + str(nr_file)
                    print(title)
                    # if spectrum['params']['title'] not in test.keys():
                    #     continue
                    xpoints_data, ypoints_data, pp_mass, dev = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'],parameter1,parameter2)
                    counter = 0
                    # if len(set([int(xs) for xs in xpoints_data]))<=pp_mass*spectrum['params']['charge'][0]/100:
                    #     xpoints_data = []
                    #     ypoints_data = []
                    intens = sum(ypoints_data)*100
                    if spectrum['params']['title'] not in test.keys():
                        intensity.append(intens)
                        x.append(len(spectrum['m/z array']))
                        y.append(len(xpoints_data))
                        mascotid.append(0)
                        break
                    else:
                        for i in sorted(test[spectrum['params']['title']]):
                            for l,n in enumerate(xpoints_data):
                                if i-0.05<=n<=i+0.05:
                                    counter +=1
                                    break
                        if test_charge ==3:
                            if counter <= result[title][1]:
                                break
                            else:
                                x = x[:-1]
                                y = y[:-1]
                                intensity = intensity[:-1]
                                mascotid = mascotid[:-1]
                        result[title]=[len(test[spectrum['params']['title']]),counter]
                        intensity.append(intens)
                        x.append(len(spectrum['m/z array']))
                        y.append(len(xpoints_data))
                        mascotid.append(1)
                if spectrum['params']['title'] in test.keys():
                    if y[-1]==0:
                        if 'frag: b ion' in ionsID[spectrum['params']['title']] and 'frag: y ion' in ionsID[spectrum['params']['title']]:
                            both_ions += 1
                        elif 'frag: b ion' in ionsID[spectrum['params']['title']]:
                            b_ions += 1
                        elif 'frag: y ion' in ionsID[spectrum['params']['title']]:
                            y_ions += 1
                        else:
                            only_neutral += 1
        col = []
        for i in mascotid:
            if i == 0:
                col.append('orange')
            else:
                col.append('blue')
        ind = np.argsort(np.array(x))
        plt.scatter(np.array(x)[ind],np.array(y)[ind],color=np.array(col)[ind],alpha = 0.1)
        plt.title('Amount of datapoints remaining after filtering \n ppm error = '+str(parameter1))
        plt.xlabel('amount of datapoints before filtering')
        plt.ylabel('Amount of datapoints after filtering')
        plt.show()
        
        mid0 = []
        mid1 = []
        for i,l in enumerate(y):
            if l==0:
                mid0.append(mascotid[i])
            else:
                mid1.append(mascotid[i])
        plt.bar(['Not Mascot and 0 \n (' +str(mid0.count(0))+')', 'Mascot and 0 \n (' +str(mid0.count(1))+')', 'Not Mascot and >0 \n (' +str(mid1.count(0))+')', 'Mascot and >0 \n (' +str(mid1.count(1))+')'],[mid0.count(0),mid0.count(1),mid1.count(0),mid1.count(1)], color=['firebrick', 'skyblue', 'limegreen','salmon'])
        plt.title('Amount of spectra after filtering \n ppm error = '+str(parameter1))
        plt.xticks(rotation=45)
        plt.show()
        
        plt.bar(['Only B ion', 'Both ions', 'Only Y ions', 'Only neutral losses'],[b_ions, both_ions, y_ions, only_neutral], color=['firebrick', 'skyblue', 'limegreen', 'salmon'])
        plt.title('Found ions fragments MASCOT on \n lost spectra after filtering ppm error = '+str(parameter1))
        plt.show()
        
        
        t = []
        for i,l in enumerate(x):
            t.append(100*y[i]/x[i])
        plt.scatter(t, intensity,color = col, alpha=0.1)
        plt.title('Remainder of total Intensity \n ppm error = '+str(parameter1))
        plt.xlabel('percentage datapoints after filtering')
        plt.ylabel('percentage intensity remaining ')
        plt.show()
        
        xas = []
        yas=[]
        for i in result.values():
            xas.append(i[0])
            yas.append(100*i[1]/i[0])
        print(yas.count(0))
        col = []
        for i,n in enumerate(xas):
            if yas[i]==0:
                col.append('firebrick')
            else:
                col.append('blue')
        plt.scatter(xas,yas,c=col,alpha=0.1)
        #plt.plot([0]+xas,[0]+xas,c='r')
        plt.xlabel('Datapoints used for MASCOT scoring on unfiltered data')
        plt.ylabel('Percentage of remaining datapoints \n used by MASCOT after filtering')
        plt.title('Loss in identifiers used by mascot \n ppm error = '+str(parameter1))
        plt.show()
        
        locals()['yas'+str(parameter1)]=yas
        
df = pd.DataFrame()   
df['count']=[np.array(yas1),np.array(yas5),np.array(yas10),np.array(yas20),np.array(yas50)]
xje=['1ppm','5ppm','10ppm','20ppm','50ppm']
import seaborn as sns
ax = sns.violinplot(df['count'])
ax.set_xticklabels(xje)
ax.set_title('Distribution of retained datapoints of original MASCOT scoring')
plt.show()
    
    
    # yas = np.array(yas)
    # ind = np.argsort(np.array(x))
    # x = np.array(x)[ind]
    # y = np.array(y)[ind]
    # yas = yas[ind]
    # t = np.array(t)[ind]
    
    # plt.scatter(yas,t,alpha = 0.1)
    # plt.xlabel('percentage of MASCOT scoring datapoints after filtering')
    # plt.ylabel('Percentage of datapoints retained after filtering')
    # plt.title('to add')
    # plt.show()
    
    
    
                    