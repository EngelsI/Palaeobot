# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 20:51:36 2023

@author: Ian
"""

import PySimpleGUI as sg
import os.path
import os
from zipfile import ZipFile
from pyteomics import mgf
import numpy as np
import time

def ghost_peaks(x, y, peptide_mass, xnew, ynew, charge):
    x = list(x)
    y=list(y)
    xnew = list(xnew)
    ynew = list(ynew)
    ghost = [num*2-1 for i, num in enumerate(x) if num <= peptide_mass*(charge/2) and num not in xnew]
    xghost = [num for i, num in enumerate(x) if num <= peptide_mass*(charge/2) and num not in xnew]
    ghost_intensity = [num for i, num in enumerate(y) if x[i]<=peptide_mass*(charge/2) and x[i] not in xnew]
    ghost = list(np.array(ghost))
    x = xghost+ghost
    
    y=ghost_intensity+ghost_intensity
    x= np.array(x)
    y=np.array(y)
    
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    return x, y

def mirror(x,y,peptide_mass, charge):
    peptide_mass = peptide_mass*charge
    dev = 10*peptide_mass/1e6
    if dev <= 0.03:
          dev = dev*1.5
    if dev >0.06:
        dev = 0.05
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[1.0073,19.0226])
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
        if iteration == 1:
            x_it1 = xnew
            y_it1 = ynew
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew, charge)

    xnew = list(x_it1) + list(xnew)
    ynew = list(y_it1) + list(ynew)
    xnew= np.array(xnew)
    ynew=np.array(ynew)
    ind = np.argsort(xnew)
    xnew = xnew[ind]
    ynew = ynew[ind]
    if len(xnew)==0:
        return xnew, ynew
    
    # iterating = True
    finalx = xnew
    finaly = ynew
    # while iterating == True:
    #     keepx = []
    #     keepy = []
    #     already = []
    #     for i in range(0,len(xnew)-1):
    #         if xnew[i+1]-xnew[i] >dev*2 and xnew[i] not in already:
    #             keepx.append(xnew[i])
    #             keepy.append(ynew[i])
    #         elif xnew[i] not in already:
    #             keepx.append((xnew[i+1]+xnew[i])/2)
    #             keepy.append(ynew[i+1]+ynew[i])
    #             already.append(xnew[i+1])
    #     keepx.append(xnew[-1])
    #     keepy.append(ynew[-1])
    #     if len(finalx)==len(keepx):
    #         iterating = False
    #     finalx = keepx
    #     finaly = keepy
    #     xnew = keepx
    #     ynew = keepy
    return np.array(finalx), np.array(finaly)

def set_data(xpoints_data, ypoints_data, pp_mass, spectrum, charge):
    ind = np.argsort(xpoints_data)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]       
    
    if charge != 2:
        pp_mass = pp_mass-((charge-2)**2)/charge

    xpoints_data, ypoints_data = mirror(xpoints_data, ypoints_data, pp_mass, charge)
    if len(xpoints_data)<=5:
        return False
    spectrum_new = {'params':{}}
    spectrum_new['intensity array'] = ypoints_data
    spectrum_new['m/z array'] = xpoints_data
    spectrum_new['charge array'] = np.ma.MaskedArray([1]*len(xpoints_data))
    spectrum_new['params']['pepmass']=(pp_mass,0)
    return spectrum_new

def w_mgf(files, name):
    name = name.split('/')
    naming = name[-1]
    name = name[:-1]
    name = '/'.join(name)
    os.chdir = name
    mgf.write(files, output='filtered_'+naming)
 
def the_program(spectrum, test_charge, title):
    to_mgf = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum, test_charge)  
    return to_mgf

if __name__ == '__main__':
   #----- Full layout -----
    layout = [
        [sg.Text("Chose a MGF file:")],
        [sg.InputText(key="-FILE_PATH-"),sg.FileBrowse(file_types=[("MGF Files", "*.mgf")])],
        [sg.Button("Submit"), sg.Exit()],[sg.ProgressBar(100, orientation='h', expand_x=True, size=(20, 20),  key='-PBAR-')],
        [sg.Text('', key='-OUT-', enable_events=True, font=('Arial Bold', 16), justification='center', expand_x=True)]
    ]
    
    window = sg.Window("MGF mirror conversion", layout)
    # Run the Event Loop
    while True:
        event, values = window.read()
        if event in ("Exit",sg.WIN_CLOSED):
            break
        elif event == "Submit":
            mgf_adress = values["-FILE_PATH-"]
            with mgf.read(mgf_adress) as reader:
                start_time = time.time()
                make_mgf = []
                nr_file =0
                for spectrum in reader:
                    nr_file += 1
                    window['-PBAR-'].update(current_count=(nr_file/len(reader))*100)
                    window['-OUT-'].update(str(nr_file)+' of the '+str(len(reader))+' spectra')
                    if 'charge' in spectrum['params']:
                        charging = spectrum['params']['charge']   
                    else:
                        charging  = [2,3]
                    title = spectrum['params']['title']
                    for test_charge in charging:
                        to_mgf = the_program(spectrum, test_charge, title)
                        if to_mgf != False:
                            to_mgf['params']['title']=str(test_charge)+'_'+title
                            to_mgf['params']['charge'] = test_charge
                            to_mgf['params']['rtinseconds'] = spectrum['params']['rtinseconds']
                            locals()[title+str(test_charge)] = to_mgf
                            make_mgf.append(locals()[title+str(test_charge)])
                w_mgf(make_mgf, mgf_adress)
            sg.popup('The conversion has finished you can exit in %s seconds'% (time.time() - start_time))
            window['-PBAR-'].update(max=100)
    window.close()