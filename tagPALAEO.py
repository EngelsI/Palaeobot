# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 10:59:14 2023

@author: Ian
"""
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
from Bio.Seq import Seq 
from Bio import pairwise2
import re
import regex
import time
##################################################################

def make_matrix(codes):
    doubles = []
    for element1 in codes.keys():
        for element2 in codes.keys():
            doubles.append(element1+'|'+element2)
    names = list(codes.keys())+doubles#+triples
    reduced_matrix = []
    for element in names:
        for element2 in names:
            if find_mass(element2.split('|'))+0.05 >= find_mass(element.split('|')) >= find_mass(element2.split('|'))-0.05:
                if element != element2:
                    reduced_matrix.append((element,element2))
                    reduced_matrix.append((element2,element))
    reduced_matrix.append(('D','N'))
    reduced_matrix.append(('N','D'))
    reduced_matrix.append(('E','Q'))
    reduced_matrix.append(('Q','E'))
    
    return reduced_matrix   

def clusters(db):
    to_do = [str(num) for num in db.keys()]
    count = 0
    clusters = {}
    for element in to_do:
        print(len(to_do), db[element])
        add = [str(element)]
        for i,element2 in enumerate(to_do):
            if element==element2 or abs(len(element)-len(element2))>25:
                continue
            seq1 = Seq(element)
            seq2 = Seq(element2)
            alignments = pairwise2.align.localmx(seq1, seq2, 3,-1)
            for alignment in sorted(alignments,key=lambda x: x[-1]):
                if 1-(alignment[0].count('-')/2)/len(alignment[0])>0.8:
                    add.append(str(element2))
                break
        clusters[count]=add
        for el in add:
            to_do.remove(el)
        if len(to_do)==0:
            break
        count += 1
    return clusters  
  
def make_db(db, error, cluster):
    error = int(error/2)
    new_db = {}
    final_db = {}
    for number in cluster.keys():
        print(number, db[cluster[number][0]])
        for old_seq in cluster[number]:
            name = db[old_seq]
            print('aligning '+name+' ...')
            if len(new_db)==0:
                new_db[old_seq]=[name]
                continue
            found = False
            for seq,n in list(new_db.items()):
                if found == True:
                    break
                seq1 = Seq(seq)
                seq2 = Seq(old_seq)
                alignments = pairwise2.align.localmx(seq1, seq2, 3,-1)
                for alignment in sorted(alignments,key=lambda x: x[-1]):
                    if 1-(alignment[0].count('-')/2)/len(alignment[0])>0.9:
                        found = True
                        gaps_to_X = alignment[0].split('-')
                        gaps_to_specific = alignment[1].split('-')
                        #Make consensus and seperate seq1 specific peptides
                        consensus = []
                        alternative = []
                        count = 0
                        for parts in gaps_to_X:
                            if len(parts)<error and '-' not in parts:
                                consensus = consensus+['-']*len(parts)
                                #adjust below?
                                if count-error >=0 and count+error+len(parts)<len(seq):
                                    alternative.append(seq[count-error:count+error+len(parts)])
                                elif count-error >=0:
                                    alternative.append(seq[count-error:])
                                else:
                                    alternative.append(seq[:count+error+len(parts)])
                                count += len(parts)
                            else:
                                consensus = consensus+[parts]
                                count+=len(parts)
                        consensus = ''.join(consensus)
                        new_db[consensus]=[name]+[n]
                        for element in alternative:
                            new_db[element]=[n]
                        #keep only seq2 specific peptides
                        alternative = []
                        count = 0
                        for parts in gaps_to_specific:
                            if len(parts)<error:
                                if count-error >=0 and count+error+len(parts)<len(old_seq):
                                    alternative.append(old_seq[count-error:count+error+len(parts)])
                                elif count-error >=0:
                                    alternative.append(old_seq[count-error:])
                                else:
                                    alternative.append(old_seq[:count+error+len(parts)])
                                count += len(parts)
                        for element in alternative:
                            if element in new_db.keys():
                                new_db[element]=new_db[element]+[name]
                            else:
                                new_db[element]=[name]
                    else:
                        new_db[old_seq] = [name]
                        found = True
                    break
        final_db.update(new_db)
        new_db = {}
    return_db = {}
    for y, x in final_db.items():
        if set(y)!={'-'} and len(y.replace('-',''))>error:
            return_db[y]=x
    return_db = adj_db(return_db)
    return return_db

def adj_db(db):
    db_new = {}
    iterate = True
    while iterate == True:
        print('adjusting size database', len(db))
        iterate = False
        for key, val in db.items():
            found = False
            key_no = key.replace('-','')
            for keys in db.keys():
                if key_no in keys and key_no != keys.replace('-',''):
                    db_new[keys]=[el for el in db[keys]]+[num for num in val]
                    found = True
                    break
            if found == False:
                db_new[key]=val
        if len(db)>len(db_new):
            iterate = True
        db = db_new.copy()
        db_new.clear()
    return db

def crap_f(path):
    print('making database')
    crap = {}
    files_own = []
    crap_concat = {}
    
    for fastafile in os.walk('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/reference_db'):#fasta_db
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta'):
                files_own.append('C:/Users/Gebruiker/Desktop/neolitic_protein_discovery/reference_db/'+i) #fasta_db
    for i in files_own:
        concat = ''
        print(i.split('/')[-1])

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
    return crap_concat, crap

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

def set_data(xpoints_data, ypoints_data, pp_mass, charge, title):
    ind = np.argsort(xpoints_data)
    xpoints_data = xpoints_data[ind]
    ypoints_data = ypoints_data[ind]
    charge = charge[ind]
    # plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
    # plt.title(title+'before')
    # plt.show()
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
    
    
    dev = 20*(spectrum['params']['charge'][0]*pp_mass)/1e6 #20ppm
    if dev <= 0.01:
          dev = 0.01
    
    if 125 >= len(xpoints_data):
        print('only ghost performed 2')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass, 1)
    else:
        print('mirroring')
        xpoints_data, ypoints_data, mdev = mirror_A(xpoints_data, ypoints_data, pp_mass)
    
    ypoints_data = tic_normalize(ypoints_data)
        
    if len(ypoints_data)>0:
        # plt.stem(xpoints_data, ypoints_data, 'b', markerfmt=' ')
        # plt.title(title)
        # plt.show()
        ind = ypoints_data >1e-4
        xpoints_data = xpoints_data[ind]
        ypoints_data = ypoints_data[ind]
    if len(set([int(xs) for xs in xpoints_data]))<=pp_mass*spectrum['params']['charge'][0]/100:
        xpoints_data = []
        ypoints_data = []
    return xpoints_data, ypoints_data, pp_mass, dev

def tic_normalize(intensity):
        
        intensity = np.array(intensity)
        return intensity / np.sum(intensity)
    
def mirror_A(x, y, peptide_mass, mdev = 0, takex=[], takey=[]):
    peptide_mass = peptide_mass*spectrum['params']['charge'][0]
    if mdev == 1:
        dev = 0.07
        takex = x
        takey = y
    elif mdev == 0:
        takex = x
        takey = y
        dev = 20*peptide_mass/1e6 
        if dev <= 0.02:
              dev = 0.02
        if dev >0.05:
            dev = 0.05
    else:
        dev = mdev*0.97
    iteration = 0
    x_it1 = []
    y_it1 = []
    while iteration != 2:
        x_withh2o = np.array(list(x)+[])#[1.0073,19.0226]
        mirror_x = sorted([(peptide_mass - u,i) for i,u in enumerate(x_withh2o)])
        new_x ={}
        for locx,element in enumerate(x):
            for value,loc in mirror_x:
                if element-dev<=value <= element+dev and ((y[loc]+y[locx])>sorted(y)[int(len(y)*0.6)] or mdev not in [0,1]):
                    if value not in new_x:
                        new_x[element]=value
                    else:
                        if abs(value-element)< abs(element-new_x[element]):
                            new_x[element]=value
        iteration += 1
        ind = [True if i in new_x.keys() else False for i in x]
        xnew = x[ind]
        ynew = y[ind]
        if iteration == 1 and (mdev == 0 or mdev == 1):
            print('ghost')
            x_it1 = xnew
            y_it1 = ynew
            x, y = ghost_peaks(x,y, peptide_mass, xnew, ynew) 

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
            if xnew[i+1]-xnew[i] >0.01 and xnew[i] not in already:#dev/2
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
    if len(finalx)>min(peptide_mass/10+25,175):
        finalx, finaly, mdev = mirror_A(finalx, finaly,peptide_mass/spectrum['params']['charge'][0], dev, takex, takey)
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
    ghost = [num*2-1 for i, num in enumerate(x) if max(AA_code.values())<=num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew]
    xghost = [num for i, num in enumerate(x) if max(AA_code.values())<=num <= peptide_mass*(spectrum['params']['charge'][0]/2) and num not in xnew]
    ghost_intensity = [num for i, num in enumerate(y) if max(AA_code.values())<=x[i]<=peptide_mass*(spectrum['params']['charge'][0]/2) and x[i] not in xnew]
    ghost = list(np.array(ghost))
    x = xghost+ghost
    
    y=ghost_intensity+ghost_intensity
    x= np.array(x)
    y=np.array(y)
    
    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]
    return x, y

def find_mass(peptide):
    mass = 0
    for AA in peptide:
        mass += AA_code[AA]
    return mass

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
                for aa2 in AA_code.keys():
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
            # if done == False and max(len(tag_seq1),len(tag_seq2))<=4:
            #     m1 = find_mass(tag_seq1)
            #     m2 = find_mass(tag_seq2)
            #     if mass-1<=tag1[0]+tag2[0]+m2+m1+abs(tag1[2]-tag2[2])<=mass+1 and tag2[0]>ion_types['Y'] and tag1[0]!=ion_types['Y']:#and abs(m1-m2)<=max(AA_code.values())+max(unimod_db.values())
            #         for aa, el in AA_code.items():
            #             if tag2[0]>ion_types['Y'] and tag1[0]!=ion_types['Y']:#float(el)<abs(tag1[2]-tag2[2])
            #                 if (tag1[0],tag_seq1+aa+tag_seq2[::-1],tag2[0],'M|'+str(abs(tag1[2]-tag2[2])-el)) not in combo:#and len(tag_seq1+aa+tag_seq2[::-1])>4
            #                             combo.append((tag1[0],tag_seq1+aa+tag_seq2[::-1],tag2[0],'M|'+str(abs(tag1[2]-tag2[2])-el)))
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

def locate_center_ptm(sequence,m,Clength):
    miss_mass = m.split('|')[-1]
    output = []
    reduced_unimod = unimod_db[unimod_db['AA'].isin(list(sequence))]
    for ptm1, massa1, AA1 in reduced_unimod[['PTM', 'mass', 'AA']].values:
        if float(miss_mass)-0.02<=massa1<=float(miss_mass)+0.02:
            if AA1 not in ['N-term', 'C-term']:
                found = re.finditer(AA1, sequence)
                for i in found:
                    final_ptm = sequence+'|'+AA1+'|'+ptm1#str(int(i.span()[0])+Clength+1)+'|'+ptm1
            if final_ptm not in output:
                output.append(final_ptm)
    if len(output)>0:
        return output
    else:
        output = 'MISS'

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
        x_array_Y.append(AA_code[AA]+x_point_Y+adding)
        x_point_Y += AA_code[AA]+adding
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
        x_array_B.append(AA_code[AA]+x_point_B+adding)
        x_point_B += AA_code[AA]+adding 
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

def ptm_combos(df):
    column_values = ['peptide','protein', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
    df_new = pd.DataFrame(columns = column_values)
    for row in df.values:
        if len(row[2])==0:
            array = np.array(row).reshape(1,-1)
            df_add = pd.DataFrame(data = array, columns = column_values, dtype='object')
            df_new = pd.concat([df_new, df_add], ignore_index = True)
        elif len(row[2])==1:
            array = np.array([[row[0],row[1], ptm, row[3], row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11]] for ptm in row[2][0][1]], dtype='object')
            df_add = pd.DataFrame(data = array, columns = column_values)
            df_new = pd.concat([df_new, df_add], ignore_index = True)
        elif len(row[2])==2:
            for ptm_N in row[2][0][1]:
                array = np.array([[row[0],row[1], ptm_N+'|'+ptm, row[3], row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11]] for ptm in row[2][1][1]], dtype='object')
                df_add = pd.DataFrame(data = array, columns = column_values)
                df_new = pd.concat([df_new, df_add], ignore_index = True)
    return df_new

def jaccard(list1, list2, pep):
    mass = find_mass(pep)
    if list(list1).count(0)>1 or len(list1)==0:
        return 0
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    j = ((float(intersection) / union))
    return j*mass/len(pep)

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

ion_types = {"Y" : 19.0226, "B" : 1.0073}

AA_code = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
            "G" : 57.02146, "H" : 137.05891, "K" : 128.09496, "L" : 113.08406, "M" : 131.04049,
            "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "R" : 156.10111, "S" : 87.03203,
            "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,"I" : 113.08406,'X':1000}


if __name__ == '__main__':
    path = 'C:/Users/Gebruiker/Desktop/neolitic_protein_discovery'#neolitic_protein_discovery
    os.chdir(path)
    MGF_zip = ZipFile('./jhendy_msconvert.zip') #LFQ_bench
    include_label = False
    df = [mgf_file.filename
           for mgf_file in MGF_zip.infolist()
           if mgf_file.filename.endswith('.mgf')or mgf_file.filename.endswith('.txt')]
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
    for i, u in AA_code.items():
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
    crap, crap_individual = crap_f(path)
    # cluster = clusters(crap)
    # crap = make_db(crap, 25, cluster)
    
    for files in df:
        print(files)
        if 'B4_01.mgf' not in files:
            continue
        column_values = ['peptide','protein', 'PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
        df_deeplc = pd.DataFrame(columns = column_values)
        nr_file =0
        with mgf.read('.\\'+files.replace('/', '\\')) as reader:
            for spectrum in reader:
                nr_file += 1
                rt_time = float(spectrum['params']['rtinseconds'])
                if 'charge' in spectrum['params']:
                    charging = spectrum['params']['charge']
                else:
                    charging  = [2,3]
                for test_charge in charging:
                    spectrum['params']['charge']=[test_charge]
                    sp_nr = str(nr_file)+'_'+str(test_charge)
                    title = 'Spectrum_' + str(nr_file)+'_'+str(test_charge)
                    print(title)
                    xpoints_data, ypoints_data, pp_mass, dev = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'], spectrum)
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
                        array = np.array([[pep[0],pep[2], pep[1], 0.95, 0.95, 0.95, rt_time,xpoints_data, ypoints_data, sp_nr, charges, spectrum['params']['title']] for pep in result], dtype='object')
                        column_values = ['peptide','protein','PTM', 'jaccard_score', 'coverage', 'pc_score', 'rt_time', 'xpoints_data', 'ypoints_data', 'spectrum_nr', 'charge_state', 'file_name']
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
            df_deeplc_output = df_deeplc
    
        with open('Tag_results_'+files.split('/')[-1]+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['Sequences corresponding to tags'])
            writer.writerow(['peptide', 'protein','PTM', 'jaccard_score', 'coverage', 'pc_score', 'spectrum_nr', 'charge_state', 'file_name'])
            df_deeplc_output['peptides']=df_deeplc_output.index
            for a1, a2, a3,a4,a5,a6,a7,a8,a9 in df_deeplc_output[['peptide', 'protein','PTM', 'jaccard_score', 'coverage', 'pc_score', 'spectrum_nr', 'charge_state', 'file_name']].values:
                a2 = []
                for i,n in crap_individual.items():
                    if a1 in i:
                        a2.append(n)
                writer.writerow([a1, a2, a3,a4,a5,a6,a7,a8,a9])
    
    
    
    
    
    # end_result = {}
    # nr_file =0
    # for files in df:
    #     print(files)
    #     # if files != df[1]:
    #     #     continue
    #     column_values = [] #To add
    #     df_deeplc = pd.DataFrame(columns = column_values) #to do
    #     with mgf.read('.\\'+files.replace('/', '\\')) as reader:
    #         for spectrum in reader:
    #             nr_file += 1
                
    #             if 'charge' in spectrum['params']:
    #                 charging = spectrum['params']['charge']
    #                 spectral_match = False
    #             else:
    #                 charging  = [2,3]
    #             for test_charge in charging:
    #                 list_1 = []
    #                 list_2 = []
                    
    #                 spectrum['params']['charge']=[test_charge]

    #                 title = 'Spectrum_' + str(nr_file)+'_'+str(test_charge)
    #                 print(title)
    #                 xpoints_data, ypoints_data, pp_mass, dev = set_data(spectrum['m/z array'], spectrum['intensity array'], list(spectrum['params']['pepmass'])[0], spectrum['charge array'])
    #                 print('Making TAG...')
    #                 TAG_list = make_TAG(xpoints_data, dev, pp_mass, test_charge)
    #                 TAG_list = sorted(TAG_list, key=lambda x:len(x[1]))[::-1]
    #                 TAG_list = TAG_list[:2000]
    #                 print('Combining TAG...')
    #                 TAG_list = combos(TAG_list, pp_mass*test_charge)
    #                 print('checking database ...')
    #                 result = find_seq(crap, TAG_list)
    #                 print('completing sequence')
    #                 result = locate_seq(result, crap)
                    
    #                 result = sorted(result, key=lambda x:len(x[2]))
    #                 # result = [num for num in result if len(num[2])==len(result[0][2])]
    #                 ptm_check = [True for num in result if 0 in num[2]]
    #                 if len(ptm_check) > 0:
    #                     result = [num for num in result if 0 in num[2]]
    #                 if len(result)>0:
    #                     end_result[spectrum['params']['title']+str(test_charge)]=result
















# def locate_seq(tag, db):
#     result = []
#     crap_t = {str(u):i for i,u in db.items()}
#     for name, peps in tag.items():
#         seq = crap_t[name]
#         for el in peps:
#             match = re.finditer(str(el[1]), str(seq))
#             for i in match:
#                 peptide = el[1]
#                 ptm=[el[3]]
#                 start=i.span()[0]
#                 end=i.span()[1]
#                 if el[0] != ion_types['B']:
#                     it = True
#                     mass = el[0]
#                     added = ''
#                     while start >1 and it == True:
#                         if seq[start-1]!='-' and mass - find_mass(seq[start-1])>0:
#                             peptide = seq[start-1]+peptide
#                             added = added+seq[start-1]
#                             mass -= find_mass(seq[start-1])
#                             start -= 1
#                         elif seq[start-1]!='-' and 0.02>mass - find_mass(seq[start-1])>-0.02:
#                             peptide = seq[start-1]+peptide
#                             added = added+seq[start-1]
#                             it=False
#                         elif seq[start-1]!='-' and mass - find_mass(seq[start-1])==0:
#                             it = False
#                         else:
#                             if 0 not in ptm and abs(mass-ion_types['B'])>0.05:
#                                 ptm = ptm+[('N',control_ptm(mass-ion_types['B'], seq[start-1],added, 'N-term'))]
#                             elif abs(mass-ion_types['B'])>0.05:
#                                 ptm = [('N',control_ptm(mass-ion_types['B'], seq[start-1],added, 'N-term'))]
#                             it = False
#                     if start<=1:
#                         if 0 not in ptm and abs(mass-ion_types['B'])>0.05:
#                             ptm = ptm+[('N',control_ptm(mass-ion_types['B'], seq[start],added, 'N-term'))]
#                         elif abs(mass-ion_types['B'])>0.05:
#                             ptm = [('N',control_ptm(mass-ion_types['B'], seq[start],added, 'N-term'))]
#                 if el[2] != ion_types['Y']:
#                     it = True
#                     mass = el[2]
#                     added = ''
#                     while end <len(seq)-1 and it == True:
#                         if seq[end] != '-' and mass - find_mass(seq[end])>0:
#                             peptide = peptide+seq[end]
#                             added = added+seq[end]
#                             mass -= find_mass(seq[end])
#                             end += 1
#                         elif seq[end] != '-' and 0.02>mass - find_mass(seq[end])>-0.02:
#                             peptide = peptide+seq[end]
#                             added = added+seq[end]
#                             it=False
#                         elif seq[end] != '-' and mass - find_mass(seq[end])==0:
#                             it = False
#                         else:
#                             if 0 not in ptm and abs(mass-ion_types['Y'])>0.05:
#                                 ptm = ptm+[('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term'))]
#                             elif abs(mass-ion_types['Y'])>0.05:
#                                 ptm = [('C',control_ptm(mass-ion_types['Y'], seq[end+1],added, 'C-term'))]
#                             it = False
#                     if end+1 == len(seq):
#                         if 0 not in ptm and abs(mass-ion_types['Y'])>0.05:
#                             ptm = ptm+[('C',control_ptm(mass-ion_types['Y'], seq[end],added, 'C-term'))]
#                         elif abs(mass-ion_types['Y'])>0.05:
#                             ptm = [('C',control_ptm(mass-ion_types['Y'], seq[end],added, 'C-term'))]
#                 if (name,peptide, ptm) not in result and ('N','MISS') not in ptm and ('C','MISS') not in ptm:
#                     result.append((name,peptide, ptm))
#     return result

# def control_ptm(miss_mass, next_AA, added_AA, term, stop=False): #C-term and N-term seperately, different for middle!
#     output = []
#     reduced_unimod = unimod_db[unimod_db['AA'].isin(list(added_AA)+[term])]
#     for ptm1, massa1, AA1 in reduced_unimod[['PTM', 'mass', 'AA']].values:#1ptm
#         if miss_mass-0.02<=massa1<=miss_mass+0.02:
#             final_ptm = ptm1+'['+AA1+']'
#             if final_ptm not in output:
#                 output.append(final_ptm)
#     if len(output)>0:
#         return output
#     for ptm1, massa1, AA1 in reduced_unimod[['PTM', 'mass', 'AA']].values:#2ptm
#         for ptm2, massa2, AA2 in reduced_unimod[['PTM', 'mass', 'AA']].values:
#             if (ptm1, massa1, AA1)==(ptm2, massa2, AA2) or (added_AA.count(AA1)<2 and AA1 == AA2):
#                 continue
#             if miss_mass-0.02<=massa1+massa2<=miss_mass+0.02:
#                 final_ptm = ptm1+'['+AA1+']'+'|'+ptm2+'['+AA2+']'
#                 if final_ptm not in output:
#                     output.append(final_ptm)
#     #missing AA and ptm
#     if stop ==True:
#         if len(output)==0:
#             output = 'MISS'
#         return output
#     miss_mass = miss_mass-AA_codes[next_AA]
#     other_unimod = unimod_db[unimod_db['AA'].isin(list(added_AA+next_AA)+[term])]
#     for ptm1, massa1, AA1 in other_unimod[['PTM', 'mass', 'AA']].values:#1ptm
#         if miss_mass-0.02<=massa1<=miss_mass+0.02:
#             final_ptm = 'add_'+next_AA+'|'+ptm1+'['+AA1+']'
#             if final_ptm not in output:
#                 output.append(final_ptm)
#     if len(output)>0:
#         return output
#     for ptm1, massa1, AA1 in other_unimod[['PTM', 'mass', 'AA']].values:#2ptm
#         for ptm2, massa2, AA2 in other_unimod[['PTM', 'mass', 'AA']].values:
#             if (ptm1, massa1, AA1)==(ptm2, massa2, AA2) or (added_AA.count(AA1)<2 and AA1 == AA2):
#                 continue
#             if miss_mass-0.02<=massa1+massa2<=miss_mass+0.02:
#                 final_ptm = 'add_'+next_AA+'|'+ptm1+'['+AA1+']'+'|'+ptm2+'['+AA2+']'
#                 if final_ptm not in output:
#                     output.append(final_ptm)
#     if len(output)==0:
#         output = 'MISS'
#     return output

# def make_TAG(x, dev, mass, charge):
#     x = list(ion_types.values())+list(x)
#     x = np.array(x)
#     min_tag = min(math.ceil(mass*charge/500)-1,5)
#     min_tag = max(min_tag,3)
#     print('minimal tag =', min_tag)
#     tag_list = []
#     for begin in x:
#         for aa1 in AA_codes.keys():
#             iterate = True
#             temp_list = [aa1]
#             while iterate==True:
#                 new_list = []
#                 for pep in temp_list:
#                     for aa2 in AA_codes.keys():
#                         new_pep = pep+aa2
#                         former = True
#                         if len(new_pep)>7:
#                             iterate = False
#                             break
#                         for xp in x[x>begin]:
#                             if xp-dev <= begin+find_mass(new_pep) <= xp+dev:
#                                 new_list.append(new_pep)
#                                 former = False
#                                 break
#                         if former == True and len(pep)>min_tag and (begin,pep,(mass*charge)-(begin+find_mass(pep))) not in tag_list and begin+find_mass(pep)<mass*charge:
#                             tag_list.append((begin,pep,(mass*charge)-(begin+find_mass(pep))))                                    
#                 if len(new_list) == 0:
#                     iterate = False
#                 else:
#                     temp_list = new_list
#             temp_list = [(begin,element,((mass*charge)-begin-find_mass(element))) for element in temp_list if len(element)>min_tag]
#             tag_list = tag_list+temp_list
#     return [num for num in set(tag_list)]
        