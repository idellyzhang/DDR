import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
import numpy
from sklearn import datasets
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC, NuSVC
from sklearn.naive_bayes import GaussianNB

from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC

from sklearn.multiclass import OneVsOneClassifier

import pprint
###########MEAN###############
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

########The SELECTION of BIOMARKERS#############
def biomarker_selection(data, size):
    top_biomarker = list(data.iloc[:,0:size]) + ['Group']
    data_marker = data[top_biomarker]
    return data_marker

###classifer Methods########
def classifer (train_X, train_Y, test_X,test_Y, method):
    cl = method.fit(train_X, train_Y)
    score = cl.score(test_X,test_Y)
    return score

#####cross validation #####################
def cross_validation(feature, annotation, method,test_size, cv):
    accuracy = []
    for x in range(cv):
        train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
        score = classifer (train_X,train_Y, test_X,test_Y, method = method)
        accuracy.append(score)
    return accuracy
##############Sample Size Selection ################
def sample_selection(data,len1,size):
    idx = np.random.randint(len(data)-len1, size = size)
    idx1 = np.random.randint(len1, size = size)
    s1 = data.ix[idx]
    s2 = data.ix[idx1+len(data)-len1]
    s = pd.concat([s1, s2])
    return s

############Evaluation Metrics##############
###classifer Methods########
def classify (train_X, train_Y, test_X,test_Y, method):
    cl = method.fit(train_X, train_Y)
    return cl


def cross_recall(feature, annotation, method,test_size, cv):
    recall = []
    for x in range(cv):
        train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
        clf = classify (train_X,train_Y, test_X,test_Y, method = method)
        y_pred = clf.predict(test_X)
        rec = recall_score(test_Y, y_pred, average='macro')
        recall.append(rec)
    return recall

def cross_accuracy(feature, annotation, method,test_size, cv):
    accuracy = []
    for x in range(cv):
        train_X, test_X, train_Y, test_Y = train_test_split(feature, annotation, test_size = test_size)
        clf = classify (train_X,train_Y, test_X,test_Y, method = method)
        y_pred = clf.predict(test_X)
        rec = accuracy_score(test_Y, y_pred)
        accuracy.append(rec)
    return accuracy

pattern = 'TEP_pan_feature.csv'
data_r = pd.read_csv(pattern)

y = data_r.Group
X = data_r.drop('Group', axis=1)

data = biomarker_selection(data_r, 596)

HC = data['Group'] == "HC"
Br = data['Group'] == "Br"
CRC = data['Group'] == "CRC"
GBM = data['Group'] == "GBM"
HBC = data['Group'] == "HBC"
lung = data['Group'] == "lung"
PAAD = data['Group'] == "PAAD"


y_HC = data[HC] .Group
X_HC = data[HC].drop('Group', axis=1)


y_Br = data[Br] .Group
X_Br = data[Br].drop('Group', axis=1)

y_CRC = data[CRC] .Group
X_CRC = data[CRC].drop('Group', axis=1)


y_GBM = data[GBM] .Group
X_GBM = data[GBM].drop('Group', axis=1)

y_HBC = data[HBC] .Group
X_HBC = data[HBC].drop('Group', axis=1)


y_lung = data[lung] .Group
X_lung = data[lung].drop('Group', axis=1)

y_PAAD = data[PAAD] .Group
X_PAAD = data[PAAD].drop('Group', axis=1)



accuracy = []
recall = []

HC_acc = []
HC_Br_p = []
HC_CRC_p = []
HC_GBM_p =[]
HC_HBC_p = []
HC_lung_p = []
HC_PAAD_p = []

Br_acc = []
Br_HC_p = []
Br_CRC_p = []
Br_GBM_p = []
Br_HBC_p = []
Br_lung_p = []
Br_PAAD_p = []

CRC_acc = []
CRC_HC_p = []
CRC_Br_p = []
CRC_GBM_p = []
CRC_HBC_p =[]
CRC_lung_p = []
CRC_PAAD_p = []

GBM_acc = []
GBM_HC_p = []
GBM_Br_p = []
GBM_CRC_p = []
GBM_HBC_p = []
GBM_lung_p = []
GBM_PAAD_p = []

HBC_acc = []
HBC_HC_p = []
HBC_Br_p = []
HBC_CRC_p = []
HBC_GBM_p = []
HBC_lung_p =[]
HBC_PAAD_p = []

lung_acc = []
lung_HC_p = []
lung_Br_p = []
lung_CRC_p = []
lung_GBM_p = []
lung_HBC_p = []
lung_PAAD_p = []

PAAD_acc = []
PAAD_HC_p = []
PAAD_Br_p = []
PAAD_CRC_p = []
PAAD_GBM_p = []
PAAD_HBC_p = []
PAAD_lung_p = []


for x in range(1000):
	    X_HC_train, X_HC_test, y_HC_train, y_HC_test = train_test_split(X_HC, y_HC,test_size=0.1)
	    X_Br_train, X_Br_test, y_Br_train, y_Br_test = train_test_split(X_Br, y_Br,test_size=0.1)
	    X_CRC_train, X_CRC_test, y_CRC_train, y_CRC_test = train_test_split(X_CRC, y_CRC,test_size=0.1)
	    X_GBM_train, X_GBM_test, y_GBM_train, y_GBM_test = train_test_split(X_GBM, y_GBM,test_size=0.1)
	    X_HBC_train, X_HBC_test, y_HBC_train, y_HBC_test = train_test_split(X_HBC, y_HBC,test_size=0.1)
	    X_lung_train, X_lung_test, y_lung_train, y_lung_test = train_test_split(X_lung, y_lung,test_size=0.1)
	    X_PAAD_train, X_PAAD_test, y_PAAD_train, y_PAAD_test = train_test_split(X_PAAD, y_PAAD,test_size=0.1)
	    
	    
	    
	    y_train = pd.concat([y_HC_train, y_Br_train, y_CRC_train, y_GBM_train,y_HBC_train,y_lung_train, y_PAAD_train])
	    X_train = pd.concat([X_HC_train, X_Br_train, X_CRC_train, X_GBM_train, X_HBC_train, X_lung_train, X_PAAD_train])
	    
	    y_test = pd.concat([y_HC_test, y_Br_test, y_CRC_test, y_GBM_test, y_HBC_test, y_lung_test, y_PAAD_test])
	    X_test = pd.concat([X_HC_test,X_Br_test, X_CRC_test, X_GBM_test, X_HBC_test, X_lung_test, X_PAAD_test])

	    y_pred = OneVsOneClassifier(LinearSVC(C=100.)).fit(X_train, y_train).predict(X_test)
	   
	    acc = accuracy_score(y_test, y_pred)

	   
	    rec = recall_score(y_test, y_pred, average='macro')
	    
	    accuracy.append(acc)
	    recall.append(rec)
	    
	    HC_pred = y_pred[0:len(y_HC_test)]
	    Br_pred = y_pred[len(y_HC_test):len(y_HC_test)+len(y_Br_test)]
	    CRC_pred = y_pred[len(y_HC_test)+len(y_Br_test):len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)]
	    GBM_pred = y_pred[len(y_HC_test)+len(y_Br_test)+len(y_CRC_test):len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test)]
	    HBC_pred = y_pred[len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test):len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test)+len(y_HBC_test)]
	    lung_pred = y_pred[len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test)+len(y_HBC_test):len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test)+len(y_HBC_test)+len(y_lung_test)]
	    PAAD_pred = y_pred[len(y_HC_test)+len(y_Br_test)+len(y_CRC_test)+len(y_GBM_test)+len(y_HBC_test)+len(y_lung_test):len(y_test)]
	    
	    #print(MU_pred)
	    HC_corr = HC_pred == "HC"
	    HC_Br = HC_pred == "Br"
	    HC_CRC = HC_pred == "CRC"
	    HC_GBM = HC_pred == "GBM"
	    HC_HBC = HC_pred == "HBC"
	    HC_lung = HC_pred == "lung"
	    HC_PAAD = HC_pred == "PAAD"
	    
	    Br_corr = Br_pred == "Br"
	    Br_HC = Br_pred == "HC"
	    Br_CRC = Br_pred == "CRC"
	    Br_GBM = Br_pred == "GBM"
	    Br_HBC = Br_pred == "HBC"
	    Br_lung = Br_pred == "lung"
	    Br_PAAD = Br_pred == "PAAD"
	    
	    CRC_corr = CRC_pred == "CRC"
	    CRC_HC = CRC_pred == "HC"
	    CRC_Br = CRC_pred == "Br"
	    CRC_GBM = CRC_pred == "GBM"
	    CRC_HBC = CRC_pred == "HBC"
	    CRC_lung = CRC_pred == "lung"
	    CRC_PAAD = CRC_pred == "PAAD"
	    
	    GBM_corr = GBM_pred == "GBM"
	    GBM_HC = GBM_pred == "HC"
	    GBM_Br = GBM_pred == "Br"
	    GBM_CRC = GBM_pred == "CRC"
	    GBM_HBC = GBM_pred == "HBC"
	    GBM_lung = GBM_pred == "lung"
	    GBM_PAAD = GBM_pred == "PAAD"
	    
	    HBC_corr = HBC_pred == "HBC"
	    HBC_HC = HBC_pred == "HC"
	    HBC_Br = HBC_pred == "Br"
	    HBC_CRC = HBC_pred == "CRC"
	    HBC_GBM = HBC_pred == "GBM"
	    HBC_lung = HBC_pred == "lung"
	    HBC_PAAD = HBC_pred == "PAAD"
	    
	    lung_corr = lung_pred == "lung"
	    lung_HC = lung_pred == "HC"
	    lung_Br = lung_pred == "Br"
	    lung_CRC = lung_pred == "CRC"
	    lung_GBM = lung_pred == "GBM"
	    lung_HBC = lung_pred == "HBC"
	    lung_PAAD = lung_pred == "PAAD"
	    
	    PAAD_corr = PAAD_pred == "PAAD"
	    PAAD_HC = PAAD_pred == "HC"
	    PAAD_Br = PAAD_pred == "Br"
	    PAAD_CRC = PAAD_pred == "CRC"
	    PAAD_GBM = PAAD_pred == "GBM"
	    PAAD_HBC = PAAD_pred == "HBC"
	    PAAD_lung = PAAD_pred == "lung"

	    HC_accuracy = len(HC_pred[HC_corr])/len(y_HC_test)
	    HC_Br_pr = len(HC_pred[HC_Br])/len(y_HC_test)
	    HC_CRC_pr = len(HC_pred[HC_CRC])/len(y_HC_test)
	    HC_GBM_pr = len(HC_pred[HC_GBM])/len(y_HC_test)
	    HC_HBC_pr = len(HC_pred[HC_HBC])/len(y_HC_test)
	    HC_lung_pr = len(HC_pred[HC_lung])/len(y_HC_test)
	    HC_PAAD_pr = len(HC_pred[HC_PAAD])/len(y_HC_test)
	    
	    Br_accuracy = len(Br_pred[Br_corr])/len(y_Br_test)
	    Br_HC_pr = len(Br_pred[Br_HC])/len(y_Br_test)
	    Br_CRC_pr = len(Br_pred[Br_CRC])/len(y_Br_test)
	    Br_GBM_pr = len(Br_pred[Br_GBM])/len(y_Br_test)
	    Br_HBC_pr = len(Br_pred[Br_HBC])/len(y_Br_test)
	    Br_lung_pr = len(Br_pred[Br_lung])/len(y_Br_test)
	    Br_PAAD_pr = len(Br_pred[Br_PAAD])/len(y_Br_test)
	    
	    CRC_accuracy = len(CRC_pred[CRC_corr])/len(y_CRC_test)
	    CRC_HC_pr = len(CRC_pred[CRC_HC])/len(y_CRC_test)
	    CRC_Br_pr = len(CRC_pred[CRC_Br])/len(y_CRC_test)
	    CRC_GBM_pr = len(CRC_pred[CRC_GBM])/len(y_CRC_test)
	    CRC_HBC_pr = len(CRC_pred[CRC_HBC])/len(y_CRC_test)
	    CRC_lung_pr = len(CRC_pred[CRC_lung])/len(y_CRC_test)
	    CRC_PAAD_pr = len(CRC_pred[CRC_PAAD])/len(y_CRC_test)
	    
	    GBM_accuracy = len(GBM_pred[GBM_corr])/len(y_GBM_test)
	    GBM_HC_pr = len(GBM_pred[GBM_HC])/len(y_GBM_test)
	    GBM_Br_pr = len(GBM_pred[GBM_Br])/len(y_GBM_test)
	    GBM_CRC_pr = len(GBM_pred[GBM_CRC])/len(y_GBM_test)
	    GBM_HBC_pr = len(GBM_pred[GBM_HBC])/len(y_GBM_test)
	    GBM_lung_pr = len(GBM_pred[GBM_lung])/len(y_GBM_test)
	    GBM_PAAD_pr = len(GBM_pred[GBM_PAAD])/len(y_GBM_test)
	    
	    HBC_accuracy = len(HBC_pred[HBC_corr])/len(y_HBC_test)
	    HBC_HC_pr = len(HBC_pred[HBC_HC])/len(y_HBC_test)
	    HBC_Br_pr = len(HBC_pred[HBC_Br])/len(y_HBC_test)
	    HBC_CRC_pr = len(HBC_pred[HBC_CRC])/len(y_HBC_test)
	    HBC_GBM_pr = len(HBC_pred[HBC_GBM])/len(y_HBC_test)
	    HBC_lung_pr = len(HBC_pred[HBC_lung])/len(y_HBC_test)
	    HBC_PAAD_pr = len(HBC_pred[HBC_PAAD])/len(y_HBC_test)
	    
	    lung_accuracy = len(lung_pred[lung_corr])/len(y_lung_test)
	    lung_HC_pr = len(lung_pred[lung_HC])/len(y_lung_test)
	    lung_Br_pr = len(lung_pred[lung_Br])/len(y_lung_test)
	    lung_CRC_pr = len(lung_pred[lung_CRC])/len(y_lung_test)
	    lung_GBM_pr = len(lung_pred[lung_GBM])/len(y_lung_test)
	    lung_HBC_pr = len(lung_pred[lung_HBC])/len(y_lung_test)
	    lung_PAAD_pr = len(lung_pred[lung_PAAD])/len(y_lung_test)
	    
	    PAAD_accuracy = len(PAAD_pred[PAAD_corr])/len(y_PAAD_test)
	    PAAD_HC_pr = len(PAAD_pred[PAAD_HC])/len(y_PAAD_test)
	    PAAD_Br_pr = len(PAAD_pred[PAAD_Br])/len(y_PAAD_test)
	    PAAD_CRC_pr = len(PAAD_pred[PAAD_CRC])/len(y_PAAD_test)
	    PAAD_GBM_pr = len(PAAD_pred[PAAD_GBM])/len(y_PAAD_test)
	    PAAD_HBC_pr = len(PAAD_pred[PAAD_HBC])/len(y_PAAD_test)
	    PAAD_lung_pr = len(PAAD_pred[PAAD_lung])/len(y_PAAD_test)
	    
	    
	    HC_acc.append(HC_accuracy)
	    HC_Br_p.append(HC_Br_pr)
	    HC_CRC_p.append(HC_CRC_pr)
	    HC_GBM_p.append(HC_GBM_pr)
	    HC_HBC_p.append(HC_HBC_pr)
	    HC_lung_p.append(HC_lung_pr)
	    HC_PAAD_p.append(HC_PAAD_pr)
	    
	    Br_acc.append(Br_accuracy)
	    Br_HC_p.append(Br_HC_pr)
	    Br_CRC_p.append(Br_CRC_pr)
	    Br_GBM_p.append(Br_GBM_pr)
	    Br_HBC_p.append(Br_HBC_pr)
	    Br_lung_p.append(Br_lung_pr)
	    Br_PAAD_p.append(Br_PAAD_pr)

	    CRC_acc.append(CRC_accuracy)
	    CRC_HC_p.append(CRC_HC_pr)
	    CRC_Br_p.append(CRC_Br_pr)
	    CRC_GBM_p.append(CRC_GBM_pr)
	    CRC_HBC_p.append(CRC_HBC_pr)
	    CRC_lung_p.append(CRC_lung_pr)
	    CRC_PAAD_p.append(CRC_PAAD_pr)
	    
	    GBM_acc.append(GBM_accuracy)
	    GBM_HC_p.append(GBM_HC_pr)
	    GBM_Br_p.append(GBM_Br_pr)
	    GBM_CRC_p.append(GBM_CRC_pr)
	    GBM_HBC_p.append(GBM_HBC_pr)
	    GBM_lung_p.append(GBM_lung_pr)
	    GBM_PAAD_p.append(GBM_PAAD_pr)
	    	    
	    HBC_acc.append(HBC_accuracy)
	    HBC_HC_p.append(HBC_HC_pr)
	    HBC_Br_p.append(HBC_Br_pr)
	    HBC_CRC_p.append(HBC_CRC_pr)
	    HBC_GBM_p.append(HBC_GBM_pr)
	    HBC_lung_p.append(HBC_lung_pr)
	    HBC_PAAD_p.append(HBC_PAAD_pr)
	    
	    lung_acc.append(lung_accuracy)
	    lung_HC_p.append(lung_HC_pr)
	    lung_Br_p.append(lung_Br_pr)
	    lung_CRC_p.append(lung_CRC_pr)
	    lung_GBM_p.append(lung_GBM_pr)
	    lung_HBC_p.append(lung_HBC_pr)
	    lung_PAAD_p.append(lung_PAAD_pr)
	    
	    PAAD_acc.append(PAAD_accuracy)
	    PAAD_HC_p.append(PAAD_HC_pr)
	    PAAD_Br_p.append(PAAD_Br_pr)
	    PAAD_CRC_p.append(PAAD_CRC_pr)
	    PAAD_GBM_p.append(PAAD_GBM_pr)
	    PAAD_HBC_p.append(PAAD_HBC_pr)
	    PAAD_lung_p.append(PAAD_lung_pr)
	    
	    #print(y_test)
print("Accuracy is : ", mean(accuracy))
print("Recall is : ", mean(recall))

print("HC Accuracy: ", mean(HC_acc))
print("HC to Br: ", mean(HC_Br_p))
print("HC to CRC: ", mean(HC_CRC_p))
print("HC to GBM: ", mean(HC_GBM_p))
print("HC to HBC: ", mean(HC_HBC_p))
print("HC to lung: ", mean(HC_lung_p))
print("HC to PAAD: ", mean(HC_PAAD_p))


print("Br Accuracy: ", mean(Br_acc))
print("Br to HC: ", mean(Br_HC_p))
print("Br to CRC: ", mean(Br_CRC_p))
print("Br to GBM: ", mean(Br_GBM_p))
print("Br to HBC: ", mean(Br_HBC_p))
print("Br to lung: ", mean(Br_lung_p))
print("Br to PAAD: ", mean(Br_PAAD_p))

print("CRC Accuracy: ", mean(CRC_acc))
print("CRC to HC: ", mean(CRC_HC_p))
print("CRC to Br: ", mean(CRC_Br_p))
print("CRC to GBM: ", mean(CRC_GBM_p))
print("CRC to HBC: ", mean(CRC_HBC_p))
print("CRC to lung: ", mean(CRC_lung_p))
print("CRC to PAAD: ", mean(CRC_PAAD_p))

print("GBM Accuracy: ", mean(GBM_acc))
print("GBM to HC: ", mean(GBM_HC_p))
print("GBM to Br: ", mean(GBM_Br_p))
print("GBM to CRC: ", mean(GBM_CRC_p))
print("GBM to HBC: ", mean(GBM_HBC_p))
print("GBM to lung: ", mean(GBM_lung_p))
print("GBM to PAAD: ", mean(GBM_PAAD_p))

print("HBC Accuracy: ", mean(HBC_acc))
print("HBC to HC: ", mean(HBC_HC_p))
print("HBC to Br: ", mean(HBC_Br_p))
print("HBC to CRC: ", mean(HBC_CRC_p))
print("HBC to GBM: ", mean(HBC_GBM_p))
print("HBC to lung: ", mean(HBC_lung_p))
print("HBC to PAAD: ", mean(HBC_PAAD_p))


print("lung Accuracy: ", mean(lung_acc))
print("lung to HC: ", mean(lung_HC_p))
print("lung to Br: ", mean(lung_Br_p))
print("lung to CRC: ", mean(lung_CRC_p))
print("lung to GBM: ", mean(lung_GBM_p))
print("lung to HBC: ", mean(lung_HBC_p))
print("lung to PAAD: ", mean(lung_PAAD_p))

print("PAAD Accuracy: ", mean(PAAD_acc))
print("PAAD to HC: ", mean(PAAD_HC_p))
print("PAAD to Br: ", mean(PAAD_Br_p))
print("PAAD to CRC: ", mean(PAAD_CRC_p))
print("PAAD to GBM: ", mean(PAAD_GBM_p))
print("PAAD to HBC: ", mean(PAAD_HBC_p))
print("PAAD to lung: ", mean(PAAD_lung_p))
