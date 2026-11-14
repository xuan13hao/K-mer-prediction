import os
import sys
import pandas as pd
import numpy as np
import pickle

from sklearn.cluster import MiniBatchKMeans, KMeans, AgglomerativeClustering
from sklearn.preprocessing import MinMaxScaler, Normalizer

from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import SGDClassifier

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import SelectKBest, mutual_info_classif

from sklearn.pipeline import Pipeline

from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import pairwise_distances
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef

import warnings
warnings.filterwarnings('ignore')

helper =  "call by\n    ~ <rs_model> <rs_sp> <sp_code> <profile> <feature1> <feature2> ..."
helper += "\t\t<rs_model> - random_state for best model checkpoint"
helper += "\t\t<rs_sp>    - random_state for training-testing split"
helper += "\t\t<sp_code>  - code for training-testing split, must from {mix, dds3, ds3, ds12}"
helper += "\t\t               mix,  [1,2,3|1,2,3], mix all samples from all three datasets"
helper += "\t\t               dds3, [1,2,3|3], using all DS1 and DS2 in addition to part of DS3, same number of training sample as [1,2,3|1,2,3]"
helper += "\t\t               ds3,  [3|3]"
helper += "\t\t               ds12, [1,2|1,2]"
helper += "\t\t<profile>  - code for profile, must from {pfam, tax, k30}"
helper += "\t\t<feature1> <feature2> ...  "
helper += "               - selected feature from given profile"

if len(sys.argv) > 5:
    model_code = sys.argv[1]
    split_code = sys.argv[2]
    ds = sys.argv[3]
    pf = sys.argv[4]
    selected = sys.argv[5:]
else:
    print(helper)
    exit()

if not model_code.isdigit():
    print("<rs_model> must be <int>")
    exit()
else:
    rs_model = int(model_code)

if not split_code.isdigit():
    print("<rs_sp> must be <int>")
    exit()
else:
    rs_sp = int(split_code)

if ds not in ['mix', 'dds3', 'ds3', 'ds12']:
    print("you must choose <sp_code> from {mix, dds3, ds3, ds12}!")
    print(helper)
    exit()

if pf not in ['pfam', 'k30', 'tax']:
    print("you must choose <profile> from {pfam, tax, k30}!")
    exit()

if pf == 'k30':
    for i,v in enumerate(selected):
        selected[i] = int(v)
"""== Parameters ============================================================"""
ff_skipNormalize = False
map_class = { "PR":1, "CR":1, "PD":0, "DD":0, "SD":1}
script_dir = os.path.dirname(os.path.abspath(__file__))
# data/ lives at the repository root, one level above analysis/
base_dir   = os.path.normpath(os.path.join(script_dir, "..", "data"))
work_space = f"{base_dir}/reproduce"
os.system("mkdir -p " + work_space) # create work_space if not exist

"""== Utility ==============================================================="""
def getDS(labels):
    gt = np.zeros(len(labels))
    for (i, v) in enumerate(labels):
        if v[0] == 'J':
            gt[i] = 0
        elif v[0] == 'E':
            gt[i] = 1
        elif v[0] == 'S':
            gt[i] = 2
    return gt

def getGT(labels):
    gt = np.zeros(len(labels))
    for (i, v) in enumerate(labels):
        gt[i] = map_class[v[-2:]]
    return gt

def normalize_feature(pd):
    ans = pd.copy()
    for col in ans.columns:
        ans[col] = (ans[col] - ans[col].mean()) / (ans[col].std()+1e-17)
    return ans

def normalize_test(pd, pd_train):
    ans = pd.copy()
    for col in ans.columns:
        ans[col] = (ans[col] - pd_train[col].mean()) / (pd_train[col].std()+1e-17)
    return ans

renameSample_test = {}
renameSample_pfs = {}
stable_list = []
with open(f"{base_dir}/metadata.csv") as fin:
    for line in fin:
        if not line[0] == '#':
            cont = line.strip().split(',')
            renameSample_test[cont[0]] = cont[0] + '_' +cont[5]
            if len(cont[12]) > 0:
                renameSample_pfs[cont[0]] = cont[0] + '_' +cont[12]
            if cont[5] == 'SD':
                stable_list.append(cont[0] + '_' +cont[5])


acc_length_pfam = {}
acc_name_pfam = {}
with open(f'{base_dir}/pfam.acclengname.csv', 'r') as fin:
    for line in fin:
        cont = line.strip().split(',')
        acc_length_pfam[cont[0]] = float(cont[1])/1000
        acc_name_pfam[cont[0]] = cont[2]

def loadPFAM(fname, ff_fast=False):
    df_pfam = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    if not ff_fast:
        for col in df_pfam.columns:
            df_pfam[col] = df_pfam[col]*1e6 / df_pfam[col].sum()
        for row in df_pfam.index:
            df_pfam.loc[row] = df_pfam.loc[row] / acc_length_pfam[row]
    df_pfam.sort_index(axis=1, inplace=True)
    return df_pfam


def loadTAX(fname):
    df_tax  = pd.DataFrame(pd.read_csv(fname, sep=',', index_col=0))
    df_tax.sort_index(axis=1, inplace=True)
    return df_tax

def loadKMC(fname):
    df_kmc  = pd.DataFrame(pd.read_csv(fname, sep='\t', index_col=0))
    df_kmc.sort_index(axis=1, inplace=True)
    return df_kmc

map_data = {}
if pf == 'pfam':
    df_pfam = loadPFAM(f"{base_dir}/pfam.count.tsv", ff_skipNormalize)
    df_pfam.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_pfam
elif pf == 'tax':
    df_tax = loadTAX(f"{base_dir}/otu.count.csv")
    df_tax.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_tax
elif pf == 'k30':
    df_kmc = loadKMC(f"{base_dir}/k30.count.tsv")
    df_kmc.rename(columns=renameSample_test,inplace=True)
    map_data[pf] = df_kmc

"""== Test bestK ============================================================"""
def pickFeature(df):
    if len(set(df.columns) & set(selected)) == 0:
        print("None of the selected feature was found in data")
    return df.T.reindex(selected).fillna(0).T.copy()

def my_training_testing_split(df_join, Y, ds='mix'):
    ds_label = getDS(df_join.index)
    if ds == 'dds3':
        df_train = df_join[ds_label!=2]
        y_train = Y[ds_label!=2]
        df_join = df_join[ds_label==2]
        Y = Y[ds_label==2]
        N = len(ds_label)
    elif ds == 'ds3':
        df_join = df_join[ds_label==2]
        Y = Y[ds_label==(int(ds[-1])-1)]
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    elif ds == 'ds12':
        df_join = df_join[(ds_label==0) | (ds_label==1)]
        Y = Y[(ds_label==0) | (ds_label==1)]
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    else:
        df_train = pd.DataFrame()
        y_train = np.array([])
        N = len(Y)
    kf = StratifiedShuffleSplit(n_splits=1, test_size=int(np.ceil(0.3*N)), random_state=rs_sp)
    for T_index, t_index in kf.split(df_join, Y):
        df_tmp , y_tmp   = df_join.iloc[T_index], Y[T_index]
        df_test,  y_test = df_join.iloc[t_index], Y[t_index]
    df_train = pd.concat([df_train, df_tmp])
    y_train = np.concatenate((y_train, y_tmp))
    return (df_train, y_train, df_test, y_test)

def crossValidateSingle(pf, ds):
    df_join = pickFeature(map_data[pf].T.copy())
    Y = getGT(df_join.index, sp)

    (df_train, y_train, df_test, y_test) = my_training_testing_split(df_join, Y, ds)
    print(df_train.shape, df_test.shape, set(y_train), set(y_test))
    X_train = normalize_feature(df_train)
    X_test  = normalize_test(df_test, df_train)

    clf = MLPClassifier(hidden_layer_sizes=(100,100), max_iter=400, random_state=rs_model, early_stopping=True)
    clf.fit(X_train, y_train)

    p_test = clf.predict(X_test)
    tn, fp, fn, tp = confusion_matrix( y_test, p_test).ravel()
    recall         = recall_score(     y_test, p_test)
    precision      = precision_score(  y_test, p_test, zero_division=1)
    f1_score       = f1_score(         y_test, p_test)
    macro_f1_score = f1_score(         y_test, p_test, average='macro')
    MCC            = matthews_corrcoef(y_test, p_test)
    print(f"TN\tFP\tFN\tTP\n{tn}\t{fp}\t{fn}\t{tp}")
    print(f"recall: {recall}, precision: {precision}, f1_score: {f1_score}")
    print(f"macro_f1_score: {macro_f1_score}, MCC: {MCC}")


crossValidateSingle(pf, ds)
