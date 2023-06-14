# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 18:32:05 2023

@author: Ian
"""
import sys
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import warnings
from Bio import SeqIO
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings('ignore')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2
import re
from matplotlib import rc
import seaborn as sns
from pymzid.read_mzid import Mzid


def crap_f(nr_file):
    
    
    crap = {}
    for record in SeqIO.parse("./neolitic_protein_discovery/crap.fasta.txt", "fasta"):
        crap[record.seq]=record.id
    for record in SeqIO.parse("./neolitic_protein_discovery/contaminants.fasta", "fasta"):
        crap[record.seq]=record.id
    lysyl = 'MHKRTYLNACLVLALAAGASQALAAPGASEMAGDVAVLQASPASTGHARFANPNAAISAAGIHFAAPPARRVARAAPLAPKPGTPLQVGVGLKTATPEIDLTTLEWIDTPDGRHTARFNISAAGAASLRAAIRLETHSGSLPDDVLLHFAGAGKEIFEASGKDLSVNRPYWSPVIEGDTLTVELVLPANLQPGDLRLSVPQVSYFADSLYKAGYRDGFGASGSCEVDAVCATQSGTRAYDNATAAVAKMVFTSSADGGSYICTGTLLNNGNSPKRQLFWSAAHCIEDQATAATLQTIWFYNTTQCYGDASTINQSVTVLTGGANILHRDAKRDTLLLELKRTPPAGVFYQGWSATPIANGSLGHDIHHPRGDAKKYSQGNVSAVGVTYDGHTALTRVDWPSAVVEGGSSGSGLLTVAGDGSYQLRGGLYGGPSYCGAPTSQRNDYFSDFSGVYSQISRYFAP'
    crap[lysyl]='lysyl'
    keratin1 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSEGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1] = 'keratin1'
    keratin1_2 = 'MSRQFSSRSGYRSGGGFSSGSAGIINYQRRTTSSSTRRSGGGGGRFSSCGGGGGSFGAGGGFGSRSLVNLGGSKSISISVARGGGRGSGFGGGYGGGGFGGGGFGGGGFGGGGIGGGGFGGFGSSGGGGFGGGGFGGGGYGGGYGPVCPPGGIQEVTINQSLLQPLNVEIDPEIQKVKSREREQIKSLNNQFASFIDKVRFLEQQNQVLQTKWELLQQVDTSTRTHNLEPYFESFINNLRRRVDQLKSDQSRLDSELKNMQDMVEDYRNKYEDEINKRTNAENEFVTIKKDVDGAYMTKVDLQAKLDNLQQEIDFLTALYQAELSQMQTQISETNVILSMDNNRSLDLDSIIAEVKAQYEDIAQKSKAEAESLYQSKYEELQITAGRHGDSVRNSKIEISELNRVIQRLRSEIDNVKKQISNLQQSISDAEQRGENALKDAKNKLNDLEDALQQAKEDLARLLRDYQELMNTKLALDLEIATYRTLLEGEESRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSGGGGGGHGSYGSGSSSGGYRGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGGVKSSGGSSSVKFVSTTYSGVTR'
    crap[keratin1_2] = 'keratin1_2'
    crap_benchmark = crap #{}
    mascot = pd.read_csv('.\mascot_output\mascot_'+nr_file+'_output.csv', sep=';')
    # for record, pr in mascot[['pep_seq', 'prot_acc']].values:
        
    #     for seq in crap.keys():
    #         if record in seq and 'unkown' not in pr:
    #             crap_benchmark[record]=crap[seq]
    #             break
    return crap_benchmark, crap, mascot

def make_venn(file_nr, unique_Palaeo, unique_casa, casa_Palaeo, unique_casaun, Palaeo_casaun, casa_casaun, all_seq, what, tester):
    v = venn3(subsets = (unique_Palaeo, unique_casa, casa_Palaeo, unique_casaun, Palaeo_casaun, casa_casaun, all_seq))#, set_labels=('Palaeobot', 'casanovo_filtered', 'casanovo_unfiltered'))
    # v.get_patch_by_id('100').set_alpha(0)
    # v.get_patch_by_id('100').set_color('white')
    #v.get_label_by_id('100').set_text('Unknown')
    if [unique_Palaeo,casa_Palaeo,Palaeo_casaun,all_seq].count(0)!=4:
        v.get_label_by_id('A').set_text('Palaeo')
    else:
        v.get_label_by_id('A').set_text('')
    if [unique_casa,casa_Palaeo,casa_casaun, all_seq].count(0)!=4:
        v.get_label_by_id('B').set_text('Mascot')
    else:
        v.get_label_by_id('B').set_text('')
    if [unique_casaun,Palaeo_casaun,casa_casaun,all_seq].count(0)!=4:
        v.get_label_by_id('C').set_text(tester)
    else:
        v.get_label_by_id('C').set_text('')
    # c = venn3(subsets = (unique_Palaeo, unique_casa, casa_Palaeo, unique_casaun, Palaeo_casaun, casa_casaun, all_seq), set_labels=('Palaeobot', 'casanovo_filtered', 'casanovo_unfiltered'))
    # c[0].set_lw(1.0)
    # c[0].set_ls('dotted')
    plt.title("Contamination found in "+file_nr)
    plt.annotate('', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
    ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),)
    #arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.show()
    return 0

def find_unique(psm_casa, psm_Palaeo, psm_casaun):
    unique_casa = 0
    for seq in set(psm_casa):
        if seq not in psm_Palaeo and seq not in psm_casaun:
            unique_casa += psm_casa.count(seq)
        else:
            i = psm_casaun.count(seq)
            u = psm_Palaeo.count(seq)
            y = max(i,u)
            if psm_casa.count(seq) > y:
                unique_casa += psm_casa.count(seq) - y
    unique_Palaeo = 0
    for seq in set(psm_Palaeo):
        
        if seq not in psm_casa and seq not in psm_casaun:
            unique_Palaeo += psm_Palaeo.count(seq)
        else:
            i = psm_casaun.count(seq)
            u = psm_casa.count(seq)
            y = max(i,u)
            if psm_Palaeo.count(seq) > y:
                unique_Palaeo += psm_Palaeo.count(seq) - y
    unique_casaun = 0
    for seq in set(psm_casaun):

        if seq not in psm_Palaeo and seq not in psm_casa:
            unique_casaun += psm_casaun.count(seq)
        else:
            i = psm_casa.count(seq)
            u = psm_Palaeo.count(seq)
            y = max(i,u)
            if psm_casaun.count(seq) > y:
                unique_casaun += psm_casaun.count(seq) - y
    return unique_casa, unique_casaun, unique_Palaeo

def find_shared(psm_casa, psm_Palaeo, psm_casaun):
    casa_Palaeo = 0
    for seq in set(list(set(psm_Palaeo).intersection(set(psm_casa)))):
        if seq not in psm_casaun:
            in_Palaeo = psm_Palaeo.count(seq)
            in_casa = psm_casa.count(seq)
            casa_Palaeo += min(in_Palaeo, in_casa)
        else:
            in_casaun = psm_casaun.count(seq)
            in_Palaeo = psm_Palaeo.count(seq)
            in_casa = psm_casa.count(seq)
            if in_casaun < min(in_Palaeo, in_casa):
                casa_Palaeo += min(in_Palaeo, in_casa)-in_casaun
    casa_casaun = 0
    for seq in set(list(set(psm_casa).intersection(set(psm_casaun)))):
        if seq not in psm_Palaeo:
            in_casa = psm_casa.count(seq)
            in_casaun = psm_casaun.count(seq)
            casa_casaun += min(in_casa, in_casaun)
        else:
            in_casaun = psm_casaun.count(seq)
            in_Palaeo = psm_Palaeo.count(seq)
            in_casa = psm_casa.count(seq)
            if in_Palaeo < min(in_casaun, in_casa):
                casa_casaun += min(in_casaun, in_casa)-in_Palaeo
    Palaeo_casaun = 0
    for seq in set(list(set(psm_Palaeo).intersection(set(psm_casaun)))):
        if seq not in psm_casa:
            in_Palaeo = psm_Palaeo.count(seq)
            in_casaun = psm_casaun.count(seq)
            Palaeo_casaun += min(in_Palaeo, in_casaun)
        else:
            in_casaun = psm_casaun.count(seq)
            in_Palaeo = psm_Palaeo.count(seq)
            in_casa = psm_casa.count(seq)
            if in_casa < min(in_Palaeo, in_casaun):
                Palaeo_casaun += min(in_Palaeo, in_casaun)-in_casa
    return casa_Palaeo, casa_casaun, Palaeo_casaun

def collect_files():
    csv_files_casa = []
    for casa_file in os.walk(path+'/casanovo_output_filtered'):
        for i in casa_file[-1]:
            if i.endswith('.csv'):
                csv_files_casa.append(path+'/casanovo_output_filtered/'+i)
    csv_files_Palaeo = []
    for Palaeo_file in os.walk(path+'/Paleobot_crap'):
        for i in Palaeo_file[-1]:
            if i.endswith('.csv'):
                csv_files_Palaeo.append(path+'/Paleobot_crap/'+i)
    csv_files_casaun = []
    for Palaeo_file in os.walk(path+'/casanovo_output_unfiltered'):
        for i in Palaeo_file[-1]:
            if i.endswith('.csv'):
                csv_files_casaun.append(path+'/casanovo_output_unfiltered/'+i)
    pepnovo = []
    for pepnovo_file in os.walk(path+'/pepnovo/unfiltered'):

        for i in pepnovo_file[-1]:
            if i.endswith('.txt'):
                pepnovo.append(path+'/pepnovo/unfiltered/'+i)
    pepnovo_filtered = []
    for pepnovo_filtered_file in os.walk(path+'/pepnovo/filtered'):
        for i in pepnovo_filtered_file[-1]:
            if i.endswith('.txt'):
                pepnovo_filtered.append(path+'/pepnovo/filtered/'+i)
    return csv_files_casa, csv_files_Palaeo, csv_files_casaun, pepnovo_filtered, pepnovo

def find_truepos(df, df2, df3, file_nr, types):
    count_mascot = 0
    FN_mascot = len(df)
    psm_mascot = []
    mascot_titles = []
    for seq, scan, pr, z in df[['pep_seq', 'pep_scan_title', 'prot_acc', 'pep_exp_z']].values:#'sequence', 'search_engine_score[1]']].values:
        if scan in mascot_titles:
            continue
        if (2<= z<=3) == False:
            continue
        mascot_titles.append(scan)
        if isinstance(seq, str) == False or isinstance(pr, str) == False:
            continue
        for i in crap_db.keys():
            if ''.join(c for c in seq if c.isalpha()==True) in i:
                count_mascot += 1
                psm_mascot.append(''.join(c for c in seq if c.isalpha()==True))
                break
    FN_mascot = FN_mascot-count_mascot
    print('PSMs found in db of mascot:'+file_nr,count_mascot, '/', len(df))#len([num for num in df['search_engine_score[1]'].values if '-' in str(num)]))
    if types == 'casanovo':
        count_casaun = 0
        FN_casaun = len(df3)
        psm_casaun = []
        for seq, sc in df3[['sequence', 'search_engine_score[1]']].values:
            if '-' in str(sc):
                FN_casaun -= 1
                continue
            for i in crap_db.keys():
                if ''.join(c for c in seq if c.isalpha()==True) in i:
                    count_casaun += 1
                    psm_casaun.append(''.join(c for c in seq if c.isalpha()==True))
                    break
        FN_casaun = FN_casaun-count_casaun
        print('PSMs found in db of casa unfiltered:'+file_nr,count_casaun, '/', len([num for num in df3['search_engine_score[1]'].values if '-' in str(num)]))
    if types == 'pepnovo':
        count_casaun = 0
        FN_casaun = len(df3)
        psm_casaun = []
        pep_titles = []
        for seq, scan, N, C in df3[['Sequence', 'Title', 'N-Gap', 'C-Gap']].values:
            if float(N)>0 or float(C)>0:
                continue
            # if '-' in str(sc):
            #     FN_casaun -= 1
            #     continue
            for i in crap_db.keys():
                if ''.join(c for c in seq if c.isalpha()==True) in i:
                    if scan in pep_titles:
                        break
                    pep_titles.append(scan)
                    count_casaun += 1
                    psm_casaun.append(''.join(c for c in seq if c.isalpha()==True))
                    break
        FN_casaun = FN_casaun-count_casaun
        print('PSMs found in db of pepnovo unfiltered:'+file_nr,count_casaun, '/', len(df3))
    if types == 'MSGF':
        count_casaun = 0
        FN_casaun = len(df3)
        psm_casaun = []

        for seq in df3['seq'].values:
            
            for i in crap_db.keys():
                if ''.join(c for c in seq if c.isalpha()==True) in i:
                    count_casaun += 1
                    psm_casaun.append(''.join(c for c in seq if c.isalpha()==True))
                    break
        FN_casaun = FN_casaun-count_casaun
        print('PSMs found in db of pepnovo unfiltered:'+file_nr,count_casaun, '/', len(df3))
        
    count_Palaeo = 0
    FN_Palaeo = len(df2)
    psm_Palaeo = []
    already = []
    for seq, nr in df2[['seq', 'nr']].values:
        for i in crap_db.keys():
            if ''.join(c for c in seq if c.isalpha()==True) in i and nr not in already:
                count_Palaeo += 1
                already.append(nr)
                psm_Palaeo.append(''.join(c for c in seq if c.isalpha()==True))
                break
    FN_Palaeo = FN_Palaeo-count_Palaeo
    print('PSMs found in db of Palaeo:'+file_nr, count_Palaeo, '/', len(df2)+no_id[file_nr])
    
    psm_mascot, psm_Palaeo, psm_casaun, crap = keep_crap(psm_mascot, psm_Palaeo, psm_casaun, crap_db)
    
    globals()['crap_db'] = crap
    
    unique_mascot, unique_casaun, unique_Palaeo = find_unique(psm_mascot, psm_Palaeo, psm_casaun)
    mascot_Palaeo, mascot_casaun, Palaeo_casaun = find_shared(psm_mascot, psm_Palaeo, psm_casaun)
    all_seq = 0
    for seq in set(psm_Palaeo):
        if seq in psm_mascot and seq in psm_casaun:
            in_Palaeo = psm_Palaeo.count(seq)
            in_mascot = psm_mascot.count(seq)
            in_casaun = psm_casaun.count(seq)
            all_seq += min(in_Palaeo, in_mascot, in_casaun)
    make_venn(file_nr, unique_Palaeo, unique_mascot, mascot_Palaeo, unique_casaun, Palaeo_casaun, mascot_casaun, all_seq, 'TP', types)
    return psm_Palaeo, psm_mascot, psm_casaun

def find_falsepos(df, df2, df3, file_nr):
    psm_casa_FN = []
    for seq, sc in df[['sequence', 'search_engine_score[1]']].values:
        done = True
        if '-' in str(sc):
            done = True
            continue
        for i in crap.keys():
            if ''.join(c for c in seq if c.isalpha()==True) in i and ''.join(c for c in seq if c.isalpha()==True) in crap_db.keys():
                done = True
                break
            elif ''.join(c for c in seq if c.isalpha()==True) in i:
                done = False
                break
        if done == False:
            psm_casa_FN.append(''.join(c for c in seq if c.isalpha()==True))
    print('False positives found of casa'+file_nr, len(psm_casa_FN))
    print(psm_casa_FN)
    psm_casaun_FN = []
    for seq, sc in df3[['sequence', 'search_engine_score[1]']].values:
        done = True
        if '-' in str(sc):
            done = True
            continue
        for i in crap.keys():
           if ''.join(c for c in seq if c.isalpha()==True) in i and ''.join(c for c in seq if c.isalpha()==True) in crap_db.keys():
               done = True
               break
           elif ''.join(c for c in seq if c.isalpha()==True) in i:
               done = False
        if done == False:
            psm_casaun_FN.append(''.join(c for c in seq if c.isalpha()==True))
    print('False positives found of casaun'+file_nr, len(psm_casaun_FN))
    print(psm_casaun_FN)
    psm_Palaeo_FN = []
    for seq in df2['seq'].values:
        done = True
        for i in crap.keys():
            if ''.join(c for c in seq if c.isalpha()==True) in i and ''.join(c for c in seq if c.isalpha()==True) in crap_db.keys():
                done = True
                break
            elif ''.join(c for c in seq if c.isalpha()==True) in i:
                done = False
        if done == False:
            psm_Palaeo_FN.append(''.join(c for c in seq if c.isalpha()==True))
    print('False positives found of Palaeo'+file_nr, len(psm_Palaeo_FN))
    print(psm_Palaeo_FN)
    unique_casa_FN, unique_casaun_FN, unique_Palaeo_FN = find_unique(psm_casa_FN, psm_Palaeo_FN, psm_casaun_FN)
    casa_Palaeo_FN, casa_casaun_FN, Palaeo_casaun_FN = find_shared(psm_casa_FN, psm_Palaeo_FN, psm_casaun_FN)
    all_seq = 0
    for seq in set(psm_Palaeo_FN):
        if seq in psm_casa_FN and seq in psm_casaun_FN:
            in_Palaeo = psm_Palaeo_FN.count(seq)
            in_casa = psm_casa_FN.count(seq)
            in_casaun = psm_casaun_FN.count(seq)
            all_seq += min(in_Palaeo, in_casa, in_casaun)
    make_venn(file_nr, unique_Palaeo_FN, unique_casa_FN, casa_Palaeo_FN, unique_casaun_FN, Palaeo_casaun_FN, casa_casaun_FN, all_seq, 'FN')
    return psm_Palaeo_FN, psm_casa_FN, psm_casaun_FN
 
def parse_pepnovo(txt_pepnovo, txt_pepnovo_filtered,file_nr):
    column_values = ['Index','RnkScr','PnvScr','N-Gap','C-Gap','[M+H]','Charge','Sequence','Title']
    pepnovo = pd.DataFrame(columns = column_values)
    title = ''
    with open([i for i in txt_pepnovo if file_nr in i][0]) as p_file:
        for line in p_file:
            if line.startswith('0'):
                line1 = line.strip()
                line1 = line1.split('\t')

                if str(title.split('_191114_FoodProteomics_')[0][-1])!=line1[6]: #consider only the good charges
                    continue
                line1 = line1 +[title]
                line1 = np.array(line1).reshape(1,-1)

                column_values = ['Index','RnkScr','PnvScr','N-Gap','C-Gap','[M+H]','Charge','Sequence', 'Title']
                df_add = pd.DataFrame(data = line1, columns = column_values)
                pepnovo = pd.concat([pepnovo, df_add], ignore_index = True)
            if line.startswith('>'):
                line2 = line.strip()
                title = line2
    column_values = ['Index','RnkScr','PnvScr','N-Gap','C-Gap','[M+H]','Charge','Sequence','Title']
    pepnovo_filtered = pd.DataFrame(columns = column_values)
    title = ''
    with open([i for i in txt_pepnovo_filtered if file_nr in i][0]) as p_file:
        for line in p_file:
            if line.startswith('0'):
                line1 = line.strip()
                line1 = line1.split('\t')

                if str(title.split('_191114_FoodProteomics_')[0][-1])!=line1[6]: #consider only the good charges
                    continue
                line1 = line1 +[title]
                line1 = np.array(line1).reshape(1,-1)

                column_values = ['Index','RnkScr','PnvScr','N-Gap','C-Gap','[M+H]','Charge','Sequence', 'Title']
                df_add = pd.DataFrame(data = line1, columns = column_values)
                pepnovo = pd.concat([pepnovo, df_add], ignore_index = True)
            if line.startswith('>'):
                line2 = line.strip()
                title = line2
    return pepnovo, pepnovo_filtered
  
def protein_coverage(mascot, casa, Palaeo, pep,msgf, crap):
    # Names of group and bar width
    dfs = {}
    dfs['Mascot']=mascot
    dfs['MSGF+']=msgf
    dfs['casanovo']=casa
    dfs['Palaeo']=Palaeo
    dfs['pep']=pep
    m_dict = {1:'Mascot',2:'casanovo',3:'Palaeo',4:'pep',5:'MSGF+'}
    names = []
    for u in set(mascot+casa+Palaeo+pep+msgf):
        for n, i in crap.items():
            if u in n:
                names.append(i)
                
    names = set(names)
    names = list(names)
    sequence_dict={}
    for seq,i in crap.items():
        if i in names:
            sequence_dict[i]=seq
    for element in names:
        column_values = ['begin','start', 'stop','end', 'denovo_type']
        temp = pd.DataFrame(columns = column_values)
        for name, lijst in dfs.items():
            for peptide in lijst:
                for i in re.finditer(peptide, str(sequence_dict[element])):
                    start = i.span()[0]
                    stop = i.span()[1]
                    array = np.array([0,start,stop, len(sequence_dict[element]),name], dtype='object').reshape(1,-1)
                    df_add = pd.DataFrame(data = array, columns = column_values)
                    temp = pd.concat([temp, df_add], ignore_index = True)
        
        ordered_df = temp#.sort_values(by='start')
        # my_range=range(1,len(temp.index)+1)
        # print(my_range) 
        # The horizontal plot is made using the hline function
        plt.hlines(y=ordered_df['denovo_type'], xmin=ordered_df['start'], xmax=ordered_df['stop'], color='grey', alpha=0.4)
        plt.hlines(y=ordered_df['denovo_type'], xmin=ordered_df['begin'], xmax=ordered_df['end'], color='grey', alpha=0.4)
        plt.scatter(ordered_df['start'], ordered_df['denovo_type'], color='skyblue', alpha=1)
        plt.scatter(ordered_df['stop'], ordered_df['denovo_type'], color='green', alpha=0.4)
        plt.scatter(ordered_df['begin'], ordered_df['denovo_type'], color='white', alpha=1)
        plt.scatter(ordered_df['end'], ordered_df['denovo_type'], color='white', alpha=0.4)

        # Add title and axis names
        #plt.yticks(ordered_df['denovo_type'], ordered_df['denovo_type'])
        plt.title("Sequence coverage of "+element, loc='center')
        plt.xlabel('protein sequence')
        plt.ylabel('algorithm')
        
        # Show the graph
        plt.show()
    return 0

def keep_crap(mascot, Palaeo, casa, crap):
    names = {}
    for u in set(mascot+casa+Palaeo):
        for n, i in crap.items():
            if u in n:
                count = 0
                if u in mascot:
                    count += 1
                if u in Palaeo:
                    count += 1
                if u in casa:
                    count += 1
                if count > 1:
                    names[n]=i
    mascot_new = []
    for i in mascot:
        for n, u in names.items():
            if i in n:
                mascot_new.append(i)
                break
    Palaeo_new = []
    for i in Palaeo:
        for n, u in names.items():
            if i in n:
                Palaeo_new.append(i)
                break
    casa_new = []
    for i in casa:
        for n, u in names.items():
            if i in n:
                casa_new.append(i)
                break
    return mascot_new, Palaeo_new, casa_new, names

def keep_casanovo(df):
    for i in range(1,len(df)+1,2):
        if df['search_engine_score[1]'][df['PSM_ID']==i].values>df['search_engine_score[1]'][df['PSM_ID']==i].values:
            df = df.drop(i)
        elif df['search_engine_score[1]'][df['PSM_ID']==i].values<df['search_engine_score[1]'][df['PSM_ID']==i].values:
            df = df.drop(i-1)
        else:
            df = df.drop(i)
    return(df)

def different_histos(psm_Palaeo, psm_mascot, psm_casanovo,psm_pepnovo,psm_msgf, crap_db):
    df = pd.DataFrame()
    Palaeo = []
    mascot= []
    casa = []
    pep = []
    msgf=[]
    peptide = []
    for i in set(psm_Palaeo+psm_mascot+psm_casanovo+psm_pepnovo+psm_msgf):
        Palaeo.append(psm_Palaeo.count(i))
        casa.append(psm_casanovo.count(i))
        pep.append(psm_pepnovo.count(i))
        mascot.append(psm_mascot.count(i))
        msgf.append(psm_msgf.count(i))
        peptide.append(i)
    df['Palaeo']=Palaeo
    df['casanovo']=casa
    df['pepnovo']=pep
    df['mascot']=mascot
    df['MSGF+']=msgf
    
    for i in range(0, len(df),19):
        df_temp = df.loc[i:i+18]
        peptides = peptide[i:i+19]
        df_temp['peptide']=peptides
        df_temp = df_temp.set_index('peptide')
        sns.heatmap(df_temp, vmin=0,vmax=max(casa+pep+Palaeo+mascot+msgf), annot=True)
        plt.title('peptide PSM abundance')
        plt.show()
     
    # done = []
    # while len(done) != len(crap_db.keys()):
    #     count = 0
    #     names=[]
    #     Palaeo =[]
    #     mascot =[]
    #     casa = []
    #     pep = []
    #     msgf = []
    #     for n, u in crap_db.items():
    #         if count==10:
    #             break
    #         if u in done:
    #             continue
    #         count += 1
    #         done.append(u)
    #         for i in set(psm_Palaeo+psm_mascot+psm_casanovo+psm_pepnovo+psm_msgf):
    #             if i in n:
    #                 names.append(u)
    #             else:
    #                 continue
    #             if i in psm_Palaeo:
    #                 Palaeo.append(1)
    #             else:
    #                 Palaeo.append(0)
    #             if i in psm_casanovo:
    #                 casa.append(1)
    #             else:
    #                 casa.append(0)
    #             if i in psm_mascot:
    #                 mascot.append(1)
    #             else:
    #                 mascot.append(0)
    #             if i in psm_msgf:
    #                 msgf.append(1)
    #             else:
    #                 msgf.append(0)
    #             if i in psm_pepnovo:
    #                 pep.append(1)
    #             else:
    #                 pep.append(0)
        
    #     r = [i for i in range(0,len(names))]
    #     barWidth = 0.75
    
    #     plt.bar(r, mascot, color='#7f6d5f', edgecolor='white', width=barWidth, label='Mascot')
    #     plt.bar(r, Palaeo, bottom=mascot, color='#557f2d', edgecolor='white', width=barWidth, label='Palaeo')
    #     temp = np.add(mascot, Palaeo).tolist()
    #     plt.bar(r, casa, bottom=temp, color='#2d7f5e', edgecolor='white', width=barWidth, label='Casanovo')
    #     temp = np.add(temp, casa).tolist()
    #     plt.bar(r, msgf, bottom=temp, color='skyblue', edgecolor='white', width=barWidth, label='MSGF+')
    #     temp = np.add(temp, msgf).tolist()
    #     plt.bar(r, pep, bottom=temp, color='y', edgecolor='white', width=barWidth, label='Pepnovo')
    
    #     plt.xticks(r, names, fontweight='bold', rotation=90)
    #     plt.xlabel("protein")
    #     plt.ylabel("unique peptides")
    #     plt.title('unique peptides found per protein per tool')
    #     plt.legend()
    #     plt.show()
    return 0

def make_length_violins(psm_Palaeo, psm_mascot, psm_casanovo,psm_pepnovo,psm_msgf, crap_db):
    df = pd.DataFrame()
    add = []
    tool = []
    for i in psm_Palaeo:
        add.append(len(i))
        tool.append('Palaeo')
    for i in psm_casanovo:
        add.append(len(i))
        tool.append('Casanovo')
    for i in psm_pepnovo:
        add.append(len(i))
        tool.append('Pepnovo')
    for i in psm_mascot:
        add.append(len(i))
        tool.append('Mascot')
    for i in psm_msgf:
        add.append(len(i))
        tool.append('MSGF+')
    df['Length of peptides']=add
    df['Annotation Tool'] = tool
    sns.violinplot(y=df["Annotation Tool"], x=df["Length of peptides"], palette = "mako")
    plt.title('Length of peptides found per tool')
    plt.show()
    return 0

def open_msgf(file):
    mzid = Mzid(".\output_msgfplus\MSGFplus_"+file+".mzid")
    mzid.read_psm()
    mzid.read_peptide()
    df = mzid.peptide_df
    return df

no_id = {'BASL':58, 'BS09':51 , 'BS16':56, 'BS18':64, 'BS23':53}

if __name__ == '__main__':
    path = 'C:/Users/Gebruiker/Desktop'
    os.chdir(path)
    
    csv_files_casa, csv_files_Palaeo, csv_files_casaun, txt_pepnovo, txt_pepnovo_filtered = collect_files()
    
    for x in csv_files_casa:
        df = pd.read_csv(x, sep = ';')
        file_nr = x.split('_')
        for i in file_nr:
            if '.csv' in i:
                file_nr = i[:-4]
        crap_db, crap, mascot = crap_f(file_nr)
        df_Palaeo = pd.read_csv([i for i in csv_files_Palaeo if file_nr in i][0], sep=';')
        df_casanovo = pd.read_csv([i for i in csv_files_casaun if file_nr in i][0], sep = ';')
        df_casanovo = keep_casanovo(df_casanovo)
        df_pepnovo, df_pepnovo_f = parse_pepnovo(txt_pepnovo, txt_pepnovo_filtered,file_nr)
        df_mascot = mascot #vervangt filtered door mascot
        df_msgf = open_msgf(file_nr)
        psm_Palaeo, psm_mascot, psm_casanovo = find_truepos(df_mascot, df_Palaeo, df_casanovo, file_nr, 'casanovo')
        psm_Palaeo, psm_mascot, psm_pepnovo = find_truepos(df_mascot, df_Palaeo, df_pepnovo, file_nr, 'pepnovo')
        psm_Palaeo, psm_mascot, psm_msgf = find_truepos(df_mascot, df_Palaeo, df_msgf, file_nr, 'MSGF')
        #psm_Palaeo_FN, psm_casa_FN, psm_casaun_FN = find_falsepos(df, df2, df3, file_nr)
        
        
        #TO DO
        
        #novor  werkt niet zonder charge state, en ook niet op mijn nieuwe files. Ookal werken al de rest er wel mee. https://www.rapidnovor.com/sequencing-protein-identification-software/

        #protein coverage 
        #protein_coverage(psm_Palaeo, psm_mascot, psm_casanovo,psm_pepnovo,psm_msgf, crap_db)
        #histogram of psms, unique, psm/protein
        different_histos(psm_Palaeo, psm_mascot, psm_casanovo,psm_pepnovo,psm_msgf, crap_db)
        #plot of peptide lengths violin
        make_length_violins(psm_Palaeo, psm_mascot, psm_casanovo,psm_pepnovo,psm_msgf, crap_db)
        
        ##ANOVA
        from scipy.stats import f_oneway
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        a_Palaeo = []
        a_mascot = []
        a_casa = []
        a_pep = []
        a_msgf = []
        for element in set(psm_Palaeo+psm_mascot+psm_casanovo+psm_pepnovo+psm_msgf):
            a_Palaeo.append(psm_Palaeo.count(element))
            a_mascot.append(psm_mascot.count(element))
            a_casa.append(psm_casanovo.count(element))
            a_pep.append(psm_pepnovo.count(element))
            a_msgf.append(psm_msgf.count(element))
        df = pd.DataFrame({'score': a_Palaeo+a_mascot+a_casa+a_pep+a_msgf,
                           'group': ['Palaeobot']*len(a_Palaeo)+['Mascot']*len(a_mascot)+['Casanovo']*len(a_casa)+['Pepnovo']*len(a_pep)+['MSGF+']*len(a_msgf)}) 
        
        from statsmodels.formula.api import ols   
        import statsmodels.api as sm
        model = ols('score ~ group', data=df).fit()
        test = sm.stats.anova_lm(model, typ=2)
        print(test)
        #f_oneway(a_Palaeo, a_mascot, a_casa, a_pep, a_msgf)
        #create DataFrame to hold data
        
        # perform Tukey's test
        tukey = pairwise_tukeyhsd(endog=df['score'],
                                  groups=df['group'],
                                  alpha=0.05)
        
        #display results
        print(tukey)
        
        df.boxplot('score', by='group')
        
        #plot of unique peptides
        a_Palaeo = []
        a_mascot = []
        a_casa = []
        a_pep = []
        a_msgf = []
        for element in set(psm_Palaeo+psm_mascot+psm_casanovo+psm_pepnovo+psm_msgf):
            if element in psm_Palaeo:
                a_Palaeo.append(1)
            else:
                a_Palaeo.append(0)
            if element in psm_mascot:
                a_mascot.append(1)
            else:
                a_mascot.append(0)
            if element in psm_casanovo:
                a_casa.append(1)
            else:
                a_casa.append(0)
            if element in psm_pepnovo:
                a_pep.append(1)
            else:
                a_pep.append(0)
            if element in psm_msgf:
                a_msgf.append(1)
            else:
                a_msgf.append(0)
        from venn import pseudovenn
        df = {'Palaeobot': set(psm_Palaeo),
              'Mascot':set(psm_mascot),
              'Casanovo':set(psm_casanovo),
              'Pepnovo':set(psm_pepnovo),
              'MSGF+':set(psm_msgf)}
        from venn import venn
        venn(df)
        plt.title('Unique peptides per algorithm')
        plt.show()
        def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
            """ 
            Annotate barplot with p-values.
        
            :param num1: number of left bar to put bracket over
            :param num2: number of right bar to put bracket over
            :param data: string to write or number for generating asterixes
            :param center: centers of all bars (like plt.bar() input)
            :param height: heights of all bars (like plt.bar() input)
            :param yerr: yerrs of all bars (like plt.bar() input)
            :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
            :param barh: bar height in axes coordinates (0 to 1)
            :param fs: font size
            :param maxasterix: maximum number of asterixes to write (for very small p-values)
            """
        
            if type(data) is str:
                text = data
            else:
                # * is p < 0.05
                # ** is p < 0.005
                # *** is p < 0.0005
                # etc.
                text = ''
                p = .05
        
                while data < p:
                    text += '*'
                    p /= 10.
        
                    if maxasterix and len(text) == maxasterix:
                        break
        
                if len(text) == 0:
                    text = 'n. s.'
        
            lx, ly = center[num1], height[num1]
            rx, ry = center[num2], height[num2]
        
            if yerr:
                ly += yerr[num1]
                ry += yerr[num2]
        
            ax_y0, ax_y1 = plt.gca().get_ylim()
            dh *= (ax_y1 - ax_y0)
            barh *= (ax_y1 - ax_y0)
        
            y = max(ly, ry) + dh
        
            barx = [lx, lx, rx, rx]
            bary = [y, y+barh, y+barh, y]
            mid = ((lx+rx)/2, y+barh)
        
            plt.plot(barx, bary, c='black')
        
            kwargs = dict(ha='center', va='bottom')
            if fs is not None:
                kwargs['fontsize'] = fs
        
            plt.text(*mid, text, **kwargs)
            
        names=['MSGF+', 'Mascot','Palaeobot', 'Casanovo', 'Pepnovo']
        color=['darkorchid', 'limegreen','firebrick', 'skyblue', 'salmon']
        colors = {names[i]:col for i,col in enumerate(color)}         
        plt.bar(names, [len(psm_msgf), len(psm_mascot),len(psm_Palaeo), len(psm_casanovo), len(psm_pepnovo)], color=color, label=names)
        plt.title('Comparison of identified PSMs')
        #plt.legend(names, loc="upper right")
        plt.ylabel('PSMs identified')
        plt.xlabel('Algorithm')
        plt.ylim(0,180)
        heights = [ len(psm_msgf), len(psm_mascot),len(psm_Palaeo), len(psm_casanovo), len(psm_pepnovo)]
        names = [0,1,2,3,4]
        barplot_annotate_brackets(0, 1, .0368, names, heights,dh=0.005)
        barplot_annotate_brackets(0, 2, .0044, names, heights, dh=0.005)
        #Wbarplot_annotate_brackets(2, 4, .0403, names, heights,dh=0.02)
        barplot_annotate_brackets(2, 4, .0087, names, heights, dh=0.08)
        plt.show()
        
        names=['Palaeobot', 'Mascot', 'Casanovo', 'Pepnovo','MSGF+']
        color=['firebrick', 'limegreen', 'skyblue', 'salmon','darkorchid']
        colors = {names[i]:col for i,col in enumerate(color)} 
        plt.bar(names, [len(set(psm_Palaeo)), len(set(psm_mascot)), len(set(psm_casanovo)), len(set(psm_pepnovo)), len(set(psm_msgf))], color=color, label=names)
        plt.title('Comparison of identified unique peptides')
        #plt.legend(names, loc="lower left")
        plt.ylabel('Unique peptides identified')
        plt.xlabel('Algorithm')
        plt.show()
        
        
        
        
        
        #plots for peakpickers
        #need: results (heatmap), peaks/spectrum (barplot), 
        names = ['IGNEQGVSR','LVGTPAEER','AETSELHTSLK','GAYVEVTAK','AVGANPEQLTR','LDSTSIPVAK',
                  'SAEGLDASASLR','YDSINNTEVSGLR','AGLIVAEGVTK','YIELAPGVDNSK','ALENDIGVPSDATVK',
                  'VGNEIQYVALR','DGTFAVDGPGVIAK','SPYVITGPGVVEYK','GFTAYYIPR','TVESLFPEEAETPGSAVR',
                  'VFTPELVDVAK','LGLDFDSFR','AVYFYAPQIPLYANK','SGGLLWQLVR']
        Sciex = [3,3,3,3,3,3,3,3,3,1,3,3,3,3,3,1,1,3,2,3]
        bars1 = [151,132,198,103,155,202,139,163,113,79,362,239,207,213,93,265,99,159,131,224,
                 151,162,205,140,158,203,155,179,242,135,65,337,219,199,204,86,229,145,133,201,
                 157,148,204,132,165,211,143,166,243,129,77,356,182,162,205,88,234,146,130,225]
        
        Msconvert = [3,3,0,3,2,3,3,3,2,3,3,2,2,3,0,3,2,3,2,2]
        bars2 = [855,903,2125,662,825,1628,524,692,527,253,3600,1839,1543,1174,374,1820,740,1322,378,2456,
                 804,928,2155,1041,984,1672,584,894,2325,1551,352,3158,1597,1285,993,277,1549,966,758,2087,
                 827,894,2080,933,839,1675,556,786,2179,1529,322,3377,1078,898,1069,326,1529,1060,373,2357]
        
        Distiller = [2,0,0,2,3,3,3,3,3,0,0,3,0,3,3,0,0,3,0,3]
        bars3 = [41,45,70,40,72,111,58,66,44,43,54,83,108,83,72,93,15,39,85,64,
                 42,46,70,40,70,78,64,76,52,28,39,55,93,109,76,64,84,39,96,58,
                 43,44,70,40,66,84,56,72,55,23,35,57,71,113,66,59,82,39,89,61]
        df = pd.DataFrame()
        df['peptide']=names
        df = df.set_index('peptide')
        df['Sciex']=Sciex
        df['MsConvert']=Msconvert
        df['Distiller']=Distiller
        sns.heatmap(df, vmin=0,vmax=3, annot=True)
        plt.title('Peptides found per Peak Picker ')
        plt.show()
        
        #amount of average datapoints
        df = pd.DataFrame()
        
        df['average amount of datapoints']=bars1+bars2+bars3
        df['peptides']=names*9
        df['Peak Picker']=['Sciex']*len(bars1)+['MsConvert']*len(bars1)+['Distiller']*len(bars1)
        a = sns.barplot(y="peptides", x="average amount of datapoints", hue="Peak Picker", data=df, ci=None)
        #a.tick_params(axis='x', rotation=90)
        a.set(title='Average amount of datapoints from different peak pickers')
        plt.show()
        
        #bar for compare dataloss
        names = ['IGNEQGVSR','LVGTPAEER','GAYVEVTAK','AVGANPEQLTR']
        bars1 = [151,132,103,155,
                 151,162,140,158,
                 157,148,132,165]
        bars4 = [41,23,30,36,
                 42, 23, 29, 36,
                 45,22,32, 36]
        bars2 = [855,903,662,825,
                 804,928,1041,984,
                 827,894,933,839]
        bars5 = [47,21,34,41,
                 47,22,33,41,
                 50,26,32,43]
        bars3 = [41,45,40,72,
                 42,46,40,70,
                 43,44,40,66]
        bars6 = [13,9,14,21,
                 18,10,17,23,
                 15,9,17,23]
                 
        df = pd.DataFrame()
        
        df['average amount of datapoints']=bars1+bars2+bars3+bars4+bars5+bars6

        df['peptides']=names*9*2
        df['Peak Picker']=['Sciex']*len(bars1)+['MsConvert']*len(bars1)+['Distiller']*len(bars1)+['Sciex_filtered']*len(bars1)+['MsConvert_filtered']*len(bars1)+['Distiller_filtered']*len(bars1)
        
        
        fig, ax = plt.subplots()
        sns.barplot(x="peptides", y="average amount of datapoints", hue="Peak Picker", data=df, ci=None)
        for i in ax.containers:
            ax.bar_label(i,rotation=90, fmt='%.1f', padding = 5)
        # sns.barplot(x="peptides", y="average amount of datapoints", hue="Peak Picker", data=df, ci=None)
        
        
        #a.tick_params(axis='x', rotation=90)
        ax.set(title='Average amount of datapoints before \n and after spectral clean-up', ylim=(0,1100))
        plt.show()
                

        # #Palaeoresults
        # # Values of each group
        # bars1 = [2, 0, 0, 0, 0] #Milk
        # bars2 = [0,0,2,0,0]#blood mammal
        # bars3 = [2,0,0,0,0]#blood fish
        # bars4 = [0,2,0,0,0]#seed protein
        
         
        # # Heights of bars1 + bars2
        # bars = np.add(bars1, bars2).tolist()
        # bars_2 = np.add(bars, bars3).tolist()
        
        # # The position of the bars on the x-axis
        # r = [0,1,2,3,4]
         
        # # Names of group and bar width
        # names = ['BASL','BS09','BS16','BS18','BS23']
        # barWidth = 0.75
         
        # # Create brown bars
        # plt.barh(names, bars1, color='skyblue', edgecolor='white')
        # # Create green bars (middle), on top of the first ones
        # plt.barh(names, bars2, left=bars1, color='firebrick', edgecolor='white')
        # # Create green bars (top)
        # plt.barh(names, bars3, left=bars, color='salmon', edgecolor='white')
        # plt.barh(names, bars4, left=bars_2, color='limegreen', edgecolor='white')
         
        # # Custom X axis
        # plt.xlabel('Number of unique identified peptides')
        # #plt.xticks(rotation=90)
        # plt.ylabel("Name file pottery shard")
        # plt.xticks([0,1,2,3,4])
        # #plt.legend([bars1, bars2, bars3, bars4], ["Milk", "mammal blood", 'fish blood', 'seed protein'], loc="upper right")
        # # Show graphic
        # plt.title('Amount of unique peptides found in food crusts')
        # plt.show()