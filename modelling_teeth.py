# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 14:00:19 2023

@author: Ian
"""

import os
from zipfile import ZipFile
import pandas as pd
import numpy as np
import csv
from sklearn.ensemble import HistGradientBoostingRegressor
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import BaggingClassifier, AdaBoostClassifier, RandomForestClassifier, VotingClassifier,  HistGradientBoostingClassifier
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings('ignore')
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from Bio import SeqIO
import re
from sklearn.svm import SVC
import scikitplot as skplt
from sklearn.decomposition import PCA

def transform_data(df):
    names = []
    for i, u in df[['Protein', 'Variable modifications ([position] description)']].values:
        if 'Label' in u:
            names.append(i+'C13label')
        else:
            names.append(i)
    df['Protein'] = np.array(names)
    df.replace('nothing', 0, inplace=True)
    samples = [col for col in df if 'Teeth' in col]
    s_id = np.array(['Protein']+samples)
    df_new = pd.DataFrame(columns = s_id)
    for protein in set(df['Protein'].values):
        count = [0]*len(samples)
        for numbers in df[samples][df['Protein']==protein].values:
            count = [sum(x) for x in zip(count, numbers)]
        count = np.array([protein]+count).reshape(1,-1)
        df_add = pd.DataFrame(data = count, columns =s_id)
        df_new = pd.concat([df_new, df_add])
    df_new=df_new.set_index('Protein')
    return df_new.transpose()

def filter_data(df):
    labeled = [col for col in df if 'C13label' in col]
    to_drop = []
    for name in labeled:
        for loc,value in enumerate(df[name].values):
            if float(value) == 0:
                to_drop.append(loc)
    names = df.index
    drop_it = []

    for i in set(to_drop):
        drop_it.append(names[i])
    df = df.drop(drop_it)
    return df

def make_fake_data(df):
    s_id = df.columns.values
    df_new = pd.DataFrame(columns = s_id)
    for i in range(0,len(df)):
        row = df.iloc[i].values.astype(float)
        for i in range(5,200,2):
            data = row*(i/100)
            df_add = pd.DataFrame(data = data.reshape(1,-1), columns =s_id)
            df_new = pd.concat([df_new, df_add], ignore_index=True)
    return df_new

def make_db(path, NO_Y):
    files_own = []
    for fastafile in os.walk(path):
        for i in fastafile[-1]:
            if i.endswith('teeth_db.txt'):
                files_own.append(path+'/'+i)
    db = {}
    for i in files_own:
        for record in SeqIO.parse(i, "fasta"):
            if 'Y' not in str(record.description) or NO_Y==False:
                db[record.seq]=record.description
    return db

def make_df(df, path, concat, db, NO_Y=False):
    df.replace(np.nan, '0', inplace=True)
    to_drop=[]
    for i in df.columns:
        if 'Raw' not in i:
            other = False
            for x in df[[i]].values:
                if 'Sequence' in x or 'Variable modifications ([position] description)' in x or "Protein" in x:
                    other = True
            if other == False:
                to_drop.append(i)
        else:
            break
    df = df.drop(to_drop, axis=1)
    to_drop=[]
    y = False
    for i in df.columns:
        for x in df[[i]].values:
            if 'Sequence' in x or 'Variable modifications ([position] description)' in x or "Protein" in x:
                y2 = False
                break
            else:
                y2 = True
        if 'Intensity' in i or 'Sample retention time (min)' in i:
            y = True
        if (y == True and y2== True):
            to_drop.append(i)
        
    df = df.drop(to_drop, axis=1)
    df.columns = df.loc[1]
    df=df.drop([0,1])
    df.set_axis(range(len(df)), inplace=True)

    to_drop = []
    for loc, x in enumerate(df['Sequence'].values):
        if len(x)==1:
            to_drop.append(loc)
    df = df.drop(to_drop)
    df.set_axis(range(len(df)), inplace=True)
    
    to_drop = []
    for loc, x in enumerate(df['Sequence'].values):
        drop = True
        for seq in db.keys():
            if x in seq:
                drop = False
                break
        if drop == True:
            to_drop.append(loc)
    df = df.drop(to_drop)
    df.set_axis(range(len(df)), inplace=True)
    df_labels = df[df['Variable modifications ([position] description)'].str.contains('Label')]
    df_labels.set_axis(range(len(df_labels)), inplace=True)
    df = df[df['Variable modifications ([position] description)'].str.contains('Label')==False]
    df.set_axis(range(len(df)), inplace=True)
    #return df
    df_t = df
    df_t=df_t.drop(columns=['Variable modifications ([position] description)', 'Protein'])
    df_t=df_t.set_axis(df['Sequence'].values)
    df_t = df_t.drop(columns=['Sequence'])
    df_t = df_t.transpose()
    
    # ######
    # names = df_t.index.values
    # labels = df_t.index.values
    # labels = ['Sample_'+num.split('_')[-1] for num in labels]
    # labels = labels[:-2]
    # names = names[:-2]
    # height = []
    # for i, name in enumerate(names):
    #     datas = df_t.loc[str(name)].values
    #     datas=datas.astype(np.float)
    #     av = np.average(datas)/np.sum(np.array([float(num) for num in df_labels.iloc[:,i]]))/len(set(df_labels['Sequence']))
    #     height.append(av) 
    
    # #names = ['Sample_'+str(nr+1) for nr, num in enumerate(names)]
    # plt.bar(labels,height,color = 'skyblue')
    # plt.title('Average normamlized peptide abundances')
    # plt.xlabel('Sample')
    # plt.ylabel('Avergae raw ion abundance')
    # plt.xticks(rotation=90)
    # plt.show()
    # ######
    
    column_values = [el+'|'+str(lo) for lo,el in enumerate(concat)]
    df_features = pd.DataFrame(columns = column_values)
    location = {}
    for peptide in df_t.columns:
        itera = re.finditer(peptide, concat)
        itera = [(tupl.span()[0], tupl.span()[1]) for tupl in itera]
        location[peptide] = itera
    for i in range(0,len(df_t)):
        row = df_t.iloc[i].values
        new_row = np.array([0]*len(concat))
        for i3,pep in enumerate(df_t.columns):
            if (pep in 'IALVLTPLK' or pep in 'WYQSMIRPPYS') and NO_Y==True:
                continue
            add_row = [0]*len(concat)
            for element in location[pep]:
                for i2 in range(element[0], element[1]+1):
                    
                    if (np.sum(np.array([float(num) for num in df_labels.iloc[:,i]]))/len(set(df_labels['Sequence'])))>0:
                        count = (float(row[i3])/len(location[pep]))/(np.sum(np.array([float(num) for num in df_labels.iloc[:,i]]))/len(set(df_labels['Sequence'])))
                    else:
                        count = (float(row[i3])/len(location[pep]))
                    if count < 0.01:
                        count = 0
                    add_row[i2]=count
                new_row = np.sum([np.array(add_row), new_row], axis=0)

        df_add = pd.DataFrame(data = new_row.reshape(1,-1), columns = column_values)
        df_features = pd.concat([df_features, df_add], ignore_index = True)
    df_features=df_features.set_axis(df_t.index)
    not_drop=[]
    for i in range(0,len(df_features.columns)):
        if np.sum(df_features.iloc[:,i])!=0:
            not_drop.append(i)
    df_features = df_features.iloc[:,not_drop]
    return df_features

def protein_coverage(seq, crap):
    # Names of group and bar width
    dfs = {}
    dfs['teeth']=seq
    
    m_dict = {1:'teeth'}
    names = []
    for u in set(seq):
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
                    array = np.array([0,start,stop, len(sequence_dict[element]),peptide], dtype='object').reshape(1,-1)
                    df_add = pd.DataFrame(data = array, columns = column_values)
                    temp = pd.concat([temp, df_add], ignore_index = True)

        ordered_df = temp#

        plt.hlines(y=ordered_df['denovo_type'], xmin=ordered_df['start'], xmax=max(ordered_df['end'].values)+10, color='grey', alpha=0.4)
        plt.hlines(y=ordered_df['denovo_type'], xmin=ordered_df['begin'], xmax=ordered_df['end'], color='grey', alpha=0.4)
        plt.scatter(ordered_df['start'], ordered_df['denovo_type'], color='skyblue', alpha=1)
        plt.scatter(ordered_df['stop'], ordered_df['denovo_type'], color='green', alpha=0.4)
        plt.scatter(ordered_df['begin'], ordered_df['denovo_type'], color='white', alpha=1)
        plt.scatter(ordered_df['end'], ordered_df['denovo_type'], color='white', alpha=0.4)
        plt.xlim(min(ordered_df['start'].values)-5,max(ordered_df['stop'].values)+5)
        plt.title("Sequence coverage of "+element, loc='center')
        plt.xlabel('Protein sequence')
        plt.ylabel('Peptide sequence')
        
        # Show the graph
        plt.show()
    return 0


if __name__ == '__main__':
    teeth_zip = ZipFile('.\\trainingset1.zip')
    
    file_list = [teeth_file.filename
           for teeth_file in teeth_zip.infolist()
           if teeth_file.filename.endswith('.csv')]
    path = 'C:/Users/Gebruiker/Desktop/neolitic_teeth'
    os.chdir(path)
    
    for file in file_list:
        NO_Y = False
        #I need males where no amely is found!!!!!!!!!!!!!!!!!!!!!!! otherwise the top feature is amely, no good accuracy when no amely considered
        db = make_db(path, NO_Y)
        #for plots
        #db={i:n for i,n in db.items() if 'Y' in n}
        concat = ''.join([str(num) for num in db.keys()])
        df = pd.read_csv(file, header=0, sep=',')
        df_features = make_df(df, path, concat, db)
        #temp = protein_coverage(df_features['Sequence'].values, db)
        
        
        
        
        
        df_features['Y'] = np.array([1,0,0,0,0,0,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,1,0])#M=1, F=0#np.array([1 if num >0 else 0 for num in df_features['I|492']])#
        NO_Y = False
        df2=pd.read_excel(path+'/Export_all_Peptide_Data_MD.xlsx', header=0)
        df_features2 = make_df(df2, path, concat, db, NO_Y)
        df_features2 = df_features2.iloc[0:len(df_features2)-2]
        df_features2['Y']=np.array([0,0,0,0,1,1,0,0,1,0,0])
        if NO_Y == True:
            df_features2 = df_features2[df_features2['Y']==1]
            print('checking only'+str(len(df_features2))+' males')
        new = pd.DataFrame()
        
        for i in df_features.columns:
            if i in df_features2.columns:
                new[i]=df_features2[i]
            else:
                new[i]=np.array([0]*len(df_features2))
        new = new.set_index(df_features2.index)
        df_features2 = new
        
        #df_features = pd.concat([df_features, df_features2], ignore_index = True)
        
        ################################################################################
        sns.countplot(data = df_features, x='Y', orient='v')
        plt.title('distribution female-male 0-1')
        plt.ylabel('# Samples')
        plt.xlabel('Class')
        plt.show()
        
        y_train = df_features['Y'].values
        df_features = df_features.drop(columns=['Y'])
        X_train = StandardScaler().fit_transform(df_features)
        y_test = df_features2['Y'].values
        df_features2 = df_features2.drop(columns=['Y'])
        X_test = StandardScaler().fit_transform(df_features2)
        
        # y=df_features['Y'].values
        # df_features = df_features.drop(columns=['Y'])
        # X_train, X_test, y_train, y_test = train_test_split(df_features, y, test_size=0.3,stratify=y, random_state=1)
        # X_test = StandardScaler().fit_transform(X_test)
        # X_train = StandardScaler().fit_transform(X_train)

        print('make PCA')
        pca = PCA(random_state=1)
        pca.fit(X_train)
        skplt.decomposition.plot_pca_component_variance(pca, figsize=(8,6), title="PCA plot of classes to expect before")
        skplt.decomposition.plot_pca_2d_projection(pca, X_train, y_train,
                                                    figsize=(10,10),
                                                    cmap="tab10")
        plt.show()
        min_features_to_select = 5  # Minimum number of features to consider
        
        from sklearn.feature_selection import RFECV
        from sklearn.model_selection import StratifiedKFold

        clf = LogisticRegression(C=1000, solver='saga', class_weight='balanced', random_state=0)
        
        cv = StratifiedKFold(3)
        
        rfecv = RFECV(
            estimator=clf,
            step=1,
            cv=cv,
            scoring="accuracy",
            min_features_to_select=min_features_to_select,
            n_jobs=2,
        )
        selector = rfecv.fit(X_train, y_train)
        ind = selector.support_
        save = []
        for loc, i in enumerate(df_features.columns):
            if ind[loc]==True:
                print(i)
                save.append(i)
        print("Optimal number of features:", rfecv.n_features_)
        X_train = selector.transform(X_train)
        X_test = selector.transform(X_test)
        
        n_scores = len(rfecv.cv_results_["mean_test_score"])
        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Mean test accuracy")
        plt.errorbar(
            range(min_features_to_select, n_scores + min_features_to_select),
            rfecv.cv_results_["mean_test_score"],
            #yerr=rfecv.cv_results_["std_test_score"],alpha=0.05
        )
        plt.title("Recursive Feature Elimination \nwith correlated features")
        plt.axvline(rfecv.n_features_, c='r', alpha=0.5)
        plt.show()
        
        # min_features_to_select = 3
        # clf = RandomForestClassifier(class_weight='balanced', random_state=0)
        # cv = StratifiedKFold(3)
        
        # rfecv = RFECV(
        #     estimator=clf,
        #     step=1,
        #     cv=cv,
        #     scoring="accuracy",
        #     min_features_to_select=min_features_to_select,
        #     n_jobs=2,
        # )
        # selector = rfecv.fit(X_train, y_train)
        # ind = selector.support_
        # new_save = []
        # for loc, i in enumerate(save):
        #     if ind[loc]==True:
        #         print(i)
        #         new_save.append(i)
        # print("Optimal number of features:", rfecv.n_features_)
        # n_scores = len(rfecv.cv_results_["mean_test_score"])
        # plt.figure()
        # plt.xlabel("Number of features selected")
        # plt.ylabel("Mean test accuracy")
        # plt.errorbar(
        #     range(min_features_to_select, n_scores + min_features_to_select),
        #     rfecv.cv_results_["mean_test_score"],
        #     yerr=rfecv.cv_results_["std_test_score"],
        # )
        # plt.title("Recursive Feature Elimination \nwith correlated features")
        # plt.axvline(rfecv.n_features_, c='r', alpha=0.5)
        # plt.show()
        # X_train = selector.transform(X_train)
        # X_test = selector.transform(X_test)
        
        print('make PCA')
        pca = PCA(random_state=1)
        pca.fit(X_train)
        skplt.decomposition.plot_pca_component_variance(pca, figsize=(8,6), title="PCA plot of classes to expect after")
        skplt.decomposition.plot_pca_2d_projection(pca, X_train, y_train,
                                                    figsize=(10,10),
                                                    cmap="tab10")
        plt.show()
        ####IALVLTPLK####WYQSMIRPPYS######
        #######################################################################################
        lr = LogisticRegression(C=1000, solver='saga', class_weight='balanced')
        lr.fit(X_train, y_train)
        preds = lr.predict(X_test)
        
        F1=lr.score(X_test, y_test)
        print(F1)
        confusion_matrix_show = confusion_matrix(y_test, preds)
        print("Confusion Matrix:\n",confusion_matrix_show)
        print("Classification Report LR:\n",classification_report(y_test, preds,zero_division=1))
        probas = lr.predict_proba(X_test)
        
        a = [num[0] if preds[i]==0 else num[1] for i, num in enumerate(probas)]

        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,preds,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) LR')
        plt.xlabel('probability score')
        plt.ylabel('class')
        for line in range(0,len(a)):
            plt.text(a[line], preds[line], df_features2.index[line], rotation='vertical')
        plt.show()
        a = [num[0] for i, num in enumerate(probas)]
        b = [num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,b,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) LR')
        plt.xlabel('probability of being female')
        plt.ylabel('probability of being male')
        labeling = list(df_features2.index.values)
        labeling = ['Sample_'+num.split('_')[-1] for num in labeling]
        
        for line in range(0,len(a)):
            plt.text(a[line], b[line], labeling[line])
        plt.show()
        if NO_Y==False:
            Y_test_probs = probas
            skplt.metrics.plot_roc_curve(y_test, Y_test_probs,
                                   title="Digits ROC Curve", figsize=(12,6))
            plt.show()
        ###################################################################################
        rf = RandomForestClassifier(class_weight='balanced', random_state=0)
        rf.fit(X_train, y_train)
        preds = rf.predict(X_test)
        
        F1=rf.score(X_test, y_test)
        print(F1)
        confusion_matrix_show = confusion_matrix(y_test, preds)
        print("Confusion Matrix:\n",confusion_matrix_show)
        print("Classification Report RF:\n",classification_report(y_test, preds,zero_division=1))
        probas = rf.predict_proba(X_test)
        
        a = [num[0] if preds[i]==0 else num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,preds,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) RF')
        plt.xlabel('probability score')
        plt.ylabel('class') 
        for line in range(0,len(a)):
            plt.text(a[line], preds[line], df_features2.index[line], rotation='vertical')
        plt.show()
        a = [num[0] for i, num in enumerate(probas)]
        b = [num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,b,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) RF')
        plt.xlabel('probability of being female')
        plt.ylabel('probability of being male')
        labeling = list(df_features2.index.values)
        labeling = ['Sample_'+num.split('_')[-1] for num in labeling]
        
        for line in range(0,len(a)):
            plt.text(a[line], b[line], labeling[line])
        plt.show()
        
        if NO_Y==False:
            Y_test_probs = probas
            skplt.metrics.plot_roc_curve(y_test, Y_test_probs,
                                   title="Digits ROC Curve", figsize=(12,6))
            plt.show()
        ######################################################################################
        from sklearn.linear_model import SGDClassifier
        lr = SGDClassifier(loss='log_loss')
        lr.fit(X_train, y_train)
        preds = lr.predict(X_test)
        
        F1=lr.score(X_test, y_test)
        print(F1)
        confusion_matrix_show = confusion_matrix(y_test, preds)
        print("Confusion Matrix:\n",confusion_matrix_show)
        print("Classification Report SGDC:\n",classification_report(y_test, preds,zero_division=1))
        probas = lr.predict_proba(X_test)
        
        a = [num[0] if preds[i]==0 else num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,preds,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) SGDC')
        plt.xlabel('probability score')
        plt.ylabel('class') 
        for line in range(0,len(a)):
            plt.text(a[line], preds[line], df_features2.index[line], rotation='vertical')
        plt.show()
        a = [num[0] for i, num in enumerate(probas)]
        b = [num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,b,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) SGDC')
        plt.xlabel('probability of being female')
        plt.ylabel('probability of being male')
        labeling = list(df_features2.index.values)
        labeling = ['Sample_'+num.split('_')[-1] for num in labeling]
        
        for line in range(0,len(a)):
            plt.text(a[line], b[line], labeling[line])
        plt.show()
        
        if NO_Y==False:
            Y_test_probs = probas
            skplt.metrics.plot_roc_curve(y_test, Y_test_probs,
                                   title="Digits ROC Curve", figsize=(12,6))
            plt.show()
        #####################################################################################
        lr = AdaBoostClassifier(svm.LinearSVC(class_weight='balanced'), algorithm="SAMME")
        lr.fit(X_train, y_train)
        preds = lr.predict(X_test)
        
        F1=lr.score(X_test, y_test)
        print(F1)
        confusion_matrix_show = confusion_matrix(y_test, preds)
        print("Confusion Matrix:\n",confusion_matrix_show)
        print("Classification Report AdaSVM:\n",classification_report(y_test, preds,zero_division=1))
        probas = lr.predict_proba(X_test)
        
        a = [num[0] if preds[i]==0 else num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,preds,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) AdaSVM')
        plt.xlabel('probability score')
        plt.ylabel('class') 
        for line in range(0,len(a)):
            plt.text(a[line], preds[line], df_features2.index[line], rotation='vertical')
        plt.show()
        a = [num[0] for i, num in enumerate(probas)]
        b = [num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,b,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) AdaSVM')
        plt.xlabel('probability of being female')
        plt.ylabel('probability of being male')
        labeling = list(df_features2.index.values)
        labeling = ['Sample_'+num.split('_')[-1] for num in labeling]
        
        for line in range(0,len(a)):
            plt.text(a[line], b[line], labeling[line])
        plt.show()
        
        if NO_Y==False:
            Y_test_probs = probas
            skplt.metrics.plot_roc_curve(y_test, Y_test_probs,
                                   title="Digits ROC Curve", figsize=(12,6))
            plt.show()
        #####################################################################################
        clf1 = RandomForestClassifier(class_weight='balanced')
        clf2 = LogisticRegression(C=1000, solver='saga', class_weight='balanced')
        clf3 = AdaBoostClassifier(svm.LinearSVC(class_weight='balanced'), algorithm="SAMME")
        #clf4 = SGDClassifier(loss='log_loss')
        lr = VotingClassifier(
            estimators=[("rf", clf1), ("lr", clf2), ("svm", clf3)],#, ('sdg',clf4)
            voting="soft",
            #weights=[2,1,2],
        )
        lr.fit(X_train, y_train)
        preds = lr.predict(X_test)
        
        F1=lr.score(X_test, y_test)
        print(F1)
        confusion_matrix_show = confusion_matrix(y_test, preds)
        print("Confusion Matrix:\n",confusion_matrix_show)
        print("Classification Report Voting:\n",classification_report(y_test, preds,zero_division=1))
        probas = lr.predict_proba(X_test)
        
        a = [num[0] if preds[i]==0 else num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,preds,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) Voting')
        plt.xlabel('probability score')
        plt.ylabel('class')
        for line in range(0,len(a)):
            plt.text(a[line], preds[line], df_features2.index[line], rotation='vertical')
        plt.show()
        
        a = [num[0] for i, num in enumerate(probas)]
        b = [num[1] for i, num in enumerate(probas)]
        colors=['b' if i==1 else 'r' for i in y_test]
        plt.scatter(a,b,color=colors, alpha=0.5)
        plt.title('Classification of male (b) female (r) Voting')
        plt.xlabel('probability of being female')
        plt.ylabel('probability of being male')
        labeling = list(df_features2.index.values)
        labeling = ['Sample_'+num.split('_')[-1] for num in labeling]
        
        for line in range(0,len(a)):
            plt.text(a[line], b[line], labeling[line])
        plt.show()

        if NO_Y==False:
            Y_test_probs = probas
            skplt.metrics.plot_roc_curve(y_test, Y_test_probs,
                                   title="Digits ROC Curve", figsize=(12,6))
            plt.show()
        
        
        


