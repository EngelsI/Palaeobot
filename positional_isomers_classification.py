# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 18:56:08 2022

@author: Ian
"""
import os
import pandas as pd
import csv
from deeplc import DeepLC
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import random
import re
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_curve, roc_auc_score

def write_csv(data):
    print('save in csv file')
    with open('records_pos_iso.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['seq', 'modifications'])
        for index in range(0,len(data)):
            writer.writerow([data.iloc[index][0], ' '])
    return 'done'

if __name__ == '__main__':
    path = 'C:/Users/Gebruiker/Desktop/neolitic_protein_discovery'
    os.chdir(path)
    data = pd.read_csv("./records_calibration.csv")
    calibration_file = "./test_train.csv"
    
    data['annotation'] = ['Original']*len(data)
    data['From']=data['seq']
    data['observed_rt']=[135, 145, 103, 128, 166, 191, 263, 222, 284,
                         155, 169, 236, 210, 195, 214, 179, 219, 198,
                         275, 233]
    
    data = data.fillna('No mod')
    data = data.loc[data["modifications"] == 'No mod']
    sequences = []
    mods = []
    annotation = []
    comes_from = []
    observed_rt = []
    random.seed(1234)
    for seq in data['seq']:
        i = 0
        while i !=20:
            for a in range(int(len(seq)/3),len(seq)+1, int(len(seq)/3)): #only ends needed in the model to be considered
                sub_seq = seq[a-int(len(seq)/3):a]
                sub_seq = random.sample(list(sub_seq), len(sub_seq))
                sub_seq = ''.join(sub_seq)
                new_seq = seq.replace(seq[a-int(len(seq)/3):a], sub_seq)
                if new_seq not in sequences and new_seq != seq:
                    sequences.append(new_seq)
                    mods.append('No mods')
                    annotation.append('positional_isomer')
                    comes_from.append(seq)
                    observed_rt.append(data['observed_rt'][data['seq']==seq].values[0])
            i += 1
        
    to_df = pd.DataFrame({'seq':sequences, 'modifications':mods, 'annotation':annotation,
                          'From':comes_from, 'observed_rt':observed_rt})
    
    data = pd.concat([data, to_df])

    to_pr = write_csv(data)  
    print(to_pr)
    
    print('Start DeepLC')
    
    peptide_file = "./records_pos_iso.csv"
    pep_df = pd.read_csv(peptide_file, sep=",")
    cal_df = pd.read_csv(calibration_file, sep=",")
    cal_df['modifications'] = cal_df['modifications'].fillna("")
    pep_df['modifications'] = pep_df['modifications'].fillna("")
    dlc = DeepLC(verbose=True)
    dlc.calibrate_preds(seq_df=cal_df)
    preds = dlc.make_preds(seq_df=pep_df)
    print('End  of deeplc')
    data['deeplc']=preds
    
    a = data['observed_rt'][data['annotation']=='Original'].values
    preds=data['deeplc'][data['annotation']=='Original'].values
    plt.figure();
    plt.suptitle('Scatterplot')
    plt.xlabel('real retention time')
    plt.ylabel('predicted retention time')
    plt.scatter(a, preds)
    z=np.polyfit(a.flatten(), preds.flatten(), 1)
    p = np.poly1d(z)
    plt.plot(a, p(a), "r--")
    plt.title("y=%.6fx+%.6f"%(z[0],z[1]))
    plt.show()
    formula1, formula2 = z[0], z[1]

    true_rt = []
    for seq in data['seq']:
        x1 = data['deeplc'][data['seq']==seq].values[0]
        true_rt.append((x1-formula2)/formula1)
    data['true_rt']=true_rt
    
    
    error = []
    annotation = []
    for sequence in data['seq']:
        dlc_orig = data['true_rt'][data['seq']==sequence].values[0]
        
        from_seq = data['From'][data['seq']==sequence].values[0]
        error_compare =(data['observed_rt'][data['seq']==from_seq].values[0]-data['true_rt'][data['seq']==from_seq].values[0])
        err = (data['observed_rt'][data['seq']==sequence].values[0]-dlc_orig)
        if data['annotation'][data['seq']==sequence].values[0]=='Original':
            annotation.append('Original')
        elif err >error_compare:
            annotation.append('positive PI')
        else:
            annotation.append('negative PI')
        error.append(err)
    
    
    data['annotation'] = annotation
    data['error']=error
    
    
    to_add = pd.DataFrame(data[data['annotation']=='Original'])
    for i in range(0,20, 5):
        i = i/100
        to_add1 = pd.DataFrame(data[data['annotation']=='Original'])
        to_add2 = pd.DataFrame(data[data['annotation']=='Original'])
        
        to_add1['error']= [x + i for x in to_add1['error']]
        to_add2['error']= [x - i for x in to_add2['error']]
        
        to_add = pd.concat([to_add, to_add1], ignore_index=True)
        to_add = pd.concat([to_add, to_add2], ignore_index=True)
    data = pd.concat([data, to_add], ignore_index=True)
    
    data_store = data
    
    data = data_store[['seq', 'annotation', 'true_rt', 'error', 'observed_rt']]
    for i in range(10, 16, 2):
        data = data_store[['annotation', 'observed_rt','error']][abs(data_store['error'])<=i]#true_rt   
        status = data['annotation']
        encoder = LabelEncoder().fit(status)
        y = encoder.transform(status)
        x = data.drop('annotation', axis=1).values
    
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.20)

        LRmodel = LogisticRegression(multi_class='ovr', C=1000)
        LRmodel.fit(X_train, y_train)
        
        predicted_classes = LRmodel.predict(X_test)
        predicted_probas = LRmodel.predict_proba(X_test)
        
        classification_accuracy = np.round(np.mean(y_test == predicted_classes)*100,2)
        
        fig, ax = plt.subplots(figsize=(15,10))
        
        colors=['#66c2a5', '#fc8d62', 'blue']

        for i, color in enumerate(colors):
            idx_train = np.where(y_train==i)
            idx_test = np.where(y_test==i)
            
            plt.scatter(X_train[idx_train,0], X_train[idx_train,1], c=color, edgecolor='black', s=30)
            plt.scatter(X_test[idx_test,0], X_test[idx_test, 1],c='white', edgecolor=color, s=70) 
        ax.legend(['Original train', 'Original test', 'Positional isomer - train', 'Positional isomer - test',
                   'Positional isomer + train', 'Positional isomer + test'])
        
        # add predictions
        for i, color in enumerate(colors):
            idx_predicted = np.where(predicted_classes==i)
            plt.scatter(X_test[idx_predicted,0], X_test[idx_predicted,1], c=color, marker='s', s=2)
        
        ax.set_xlabel('retention time')
        ax.set_ylabel('error')
        ax.set_title('classification accuracy: {}%'.format(classification_accuracy)).set_fontsize(20)
        plt.show()
        
    print('Start DeepLC test on spectrum')
    
    peptide_file = "./records_Spectrum_2.csv"
    pep_df = pd.read_csv(peptide_file, sep=",")
    cal_df = pd.read_csv(calibration_file, sep=",")
    cal_df['modifications'] = cal_df['modifications'].fillna("")
    pep_df['modifications'] = pep_df['modifications'].fillna("")
    dlc = DeepLC(verbose=True)
    dlc.calibrate_preds(seq_df=cal_df)
    preds = dlc.make_preds(seq_df=pep_df)
    print('End  of deeplc')
    Spectrum_data=pd.DataFrame({'seq':pep_df['seq'].values, 'deeplc':preds})
    true_rt = []
    for seq in Spectrum_data['seq']:
        x1 = Spectrum_data['deeplc'][Spectrum_data['seq']==seq].values[0]
        true_rt.append((x1-formula2)/formula1) 
    Spectrum_data['true_rt']=true_rt
    observed = [128]*len(Spectrum_data)
    Spectrum_data['observed_rt']=observed
    error = []
    for sequence in Spectrum_data['seq']:
        dlc_orig = Spectrum_data['true_rt'][Spectrum_data['seq']==sequence].values[0]
        err = (Spectrum_data['observed_rt'][Spectrum_data['seq']==sequence].values[0]-dlc_orig)
        error.append(err)
    Spectrum_data['error']=error
    Spectrum_data_test = Spectrum_data[['observed_rt', 'error']]
    predicted_probas_spectrum = LRmodel.predict_proba(Spectrum_data_test)
    predicted_probas_spectrum = pd.DataFrame(predicted_probas_spectrum)
    Spectrum_data['Class0_pred']=predicted_probas_spectrum[0]
        