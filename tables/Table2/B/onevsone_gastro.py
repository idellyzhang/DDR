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

pattern = 'feature_gastro.csv'
data_r = pd.read_csv(pattern)
#sample = sample_selection(data,39,50)

y = data_r.Group
X = data_r.drop('Group', axis=1)

data = biomarker_selection(data_r, 144)

CRC = data['Group'] == "CRC"
HBC = data['Group'] == "HBC"
PAAD = data['Group'] == "PAAD"

y_CRC = data[CRC] .Group
X_CRC = data[CRC].drop('Group', axis=1)

y_HBC = data[HBC] .Group
X_HBC = data[HBC].drop('Group', axis=1)

y_PAAD = data[PAAD] .Group
X_PAAD = data[PAAD].drop('Group', axis=1)


accuracy = []
recall = []


CRC_acc = []
CRC_HBC_p =[]
CRC_PAAD_p = []

HBC_acc = []
HBC_CRC_p = []
HBC_PAAD_p = []

PAAD_acc = []
PAAD_CRC_p = []
PAAD_HBC_p = []


for x in range(1000):
	    X_CRC_train, X_CRC_test, y_CRC_train, y_CRC_test = train_test_split(X_CRC, y_CRC,test_size=0.1)
	    X_HBC_train, X_HBC_test, y_HBC_train, y_HBC_test = train_test_split(X_HBC, y_HBC,test_size=0.1)
	    X_PAAD_train, X_PAAD_test, y_PAAD_train, y_PAAD_test = train_test_split(X_PAAD, y_PAAD,test_size=0.1)
	    
	    
	    
	    y_train = pd.concat([y_CRC_train, y_HBC_train, y_PAAD_train])
	    X_train = pd.concat([X_CRC_train, X_HBC_train, X_PAAD_train])
	    
	    y_test = pd.concat([y_CRC_test, y_HBC_test, y_PAAD_test])
	    X_test = pd.concat([X_CRC_test, X_HBC_test, X_PAAD_test])
	    

	    y_pred = OneVsOneClassifier(LinearSVC(C=100.)).fit(X_train, y_train).predict(X_test)

	    acc = accuracy_score(y_test, y_pred)

	   
	    rec = recall_score(y_test, y_pred, average='macro')
	    
	    accuracy.append(acc)
	    recall.append(rec)
	    
	    CRC_pred = y_pred[0:len(y_CRC_test)]
	    HBC_pred = y_pred[len(y_CRC_test):len(y_CRC_test)+len(y_HBC_test)]
	    PAAD_pred = y_pred[len(y_CRC_test)+len(y_HBC_test):len(y_CRC_test)+len(y_HBC_test)+len(y_PAAD_test)]
	    
	    #print(MU_pred)
	    CRC_corr = CRC_pred == "CRC"
	    CRC_HBC = CRC_pred == "HBC"
	    CRC_PAAD = CRC_pred == "PAAD"
	    
	    HBC_corr = HBC_pred == "HBC"
	    HBC_CRC = HBC_pred == "CRC"
	    HBC_PAAD = HBC_pred == "PAAD"
	    
	    PAAD_corr = PAAD_pred == "PAAD"
	    PAAD_CRC = PAAD_pred == "CRC"
	    PAAD_HBC = PAAD_pred == "HBC"

	    
	    CRC_accuracy = len(CRC_pred[CRC_corr])/len(y_CRC_test)
	    CRC_HBC_pr = len(CRC_pred[CRC_HBC])/len(y_CRC_test)
	    CRC_PAAD_pr = len(CRC_pred[CRC_PAAD])/len(y_CRC_test)
	    
	    HBC_accuracy = len(HBC_pred[HBC_corr])/len(y_HBC_test)
	    HBC_CRC_pr = len(HBC_pred[HBC_CRC])/len(y_HBC_test)
	    HBC_PAAD_pr = len(HBC_pred[HBC_PAAD])/len(y_HBC_test)

	    PAAD_accuracy = len(PAAD_pred[PAAD_corr])/len(y_PAAD_test)
	    PAAD_CRC_pr = len(PAAD_pred[PAAD_CRC])/len(y_PAAD_test)
	    PAAD_HBC_pr = len(PAAD_pred[PAAD_HBC])/len(y_PAAD_test)
	    
	    CRC_acc.append(CRC_accuracy)
	    CRC_HBC_p.append(CRC_HBC_pr)
	    CRC_PAAD_p.append(CRC_PAAD_pr)

	    HBC_acc.append(HBC_accuracy)
	    HBC_CRC_p.append(HBC_CRC_pr)
	    HBC_PAAD_p.append(HBC_PAAD_pr)

	    PAAD_acc.append(PAAD_accuracy)
	    PAAD_CRC_p.append(PAAD_CRC_pr)
	    PAAD_HBC_p.append(PAAD_HBC_pr)
	    

print("Accuracy is : ", mean(accuracy))
print("Recall is : ", mean(recall))

print("CRC Accuracy: ", mean(CRC_acc))
print("CRC to HBC: ", mean(CRC_HBC_p))
print("CRC to PAAD: ", mean(CRC_PAAD_p))

print("HBC Accuracy: ", mean(HBC_acc))
print("HBC to CRC: ", mean(HBC_CRC_p))
print("HBC to PAAD: ", mean(HBC_PAAD_p))

print("PAAD Accuracy: ", mean(PAAD_acc))
print("PAAD to CRC: ", mean(PAAD_CRC_p))
print("PAAD to HBC: ", mean(PAAD_HBC_p))
