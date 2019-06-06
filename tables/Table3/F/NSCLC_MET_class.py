import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn import datasets
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC, NuSVC
from sklearn.naive_bayes import GaussianNB

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

pattern = 'NSCLC_MET_feature.csv'
data_r = pd.read_csv(pattern)

data = biomarker_selection(data_r,10)

MU = data['Group'] == "MET"
WT = data['Group'] == "WT"


y_MU = data[MU] .Group
X_MU = data[MU].drop('Group', axis=1)


y_WT = data[WT] .Group
X_WT = data[WT].drop('Group', axis=1)


accuracy = []
recall = []
MU_acc = []
WT_acc = []
for x in range(1000):
	    X_MU_train, X_MU_test, y_MU_train, y_MU_test = train_test_split(X_MU, y_MU,test_size=0.1)
	    X_WT_train, X_WT_test, y_WT_train, y_WT_test = train_test_split(X_WT, y_WT,test_size=0.1)
	    
	    y_train = pd.concat([y_MU_train, y_WT_train])
	    X_train = pd.concat([X_MU_train, X_WT_train])
	    y_test = pd.concat([y_MU_test, y_WT_test])
	    X_test = pd.concat([X_MU_test, X_WT_test])
	    clf = svm.SVC(kernel='rbf', C=1, gamma = 'auto').fit(X_train, y_train)

	    y_pred = clf.predict(X_test)

	    acc = accuracy_score(y_test, y_pred)

	   
	    rec = recall_score(y_test, y_pred, average='macro')
	    accuracy.append(acc)
	    recall.append(rec)
	    MU_pred = y_pred[0:len(y_MU_test)]
	    WT_pred = y_pred[len(y_MU_test):len(y_test)]
	    #print(MU_pred)
	    MU_corr = MU_pred == "MET"
	    #print(MU_pred[MU_corr])
	    WT_corr = WT_pred == "WT"
	    MU_accuracy = len(MU_pred[MU_corr])/len(y_MU_test)
	    WT_accuracy = len(WT_pred[WT_corr])/len(y_WT_test)
	    MU_acc.append(MU_accuracy)
	    WT_acc.append(WT_accuracy)
	    #print(y_test)
print("Accuracy: ", mean(accuracy))
print("Recall: ", mean(recall))

print("MET+: ", mean(MU_acc))
print("MET-: ", mean(WT_acc))

