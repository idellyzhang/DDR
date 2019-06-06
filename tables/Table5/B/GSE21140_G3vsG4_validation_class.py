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
from sklearn.multiclass import OneVsRestClassifier

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

pattern = 'GSE37418_G3vsG4_feature.csv'
data = pd.read_csv(pattern)


y = data .Group
X = data.drop('Group', axis=1)

test = 'GSE21140_G3vsG4_validation_feature.csv'
data_test = pd.read_csv(test)
y_test = data_test.Group
X_test = data_test.drop('Group', axis=1)


clf = OneVsRestClassifier(svm.SVC(kernel='rbf', C=1, gamma = 'auto')).fit(X, y)
y_pred = clf.predict(X_test)

G3_pred = y_pred[0:27]
G4_pred = y_pred[27:62]

G3_corr = G3_pred == "G3"
G3_G4 = G3_pred == "G4"
	    #print(MU_pred[MU_corr])
G4_corr = G4_pred == "G4"
G4_G3 = G4_pred == "G3"
	    


G3_accuracy = len(G3_pred[G3_corr])/27
G3_G4 = len(G3_pred[G3_G4])/27

G4_accuracy = len(G4_pred[G4_corr])/35
G4_G3 = len(G4_pred[G4_G3])/35


acc = accuracy_score(y_test, y_pred)
rec = recall_score(y_test, y_pred, average='macro')
#print(y_pred)
print("Overall accuracy:", acc)
print("Recall:", rec)

print("G3 accuracy:", G3_accuracy)
print("G4 accuracy:", G4_accuracy)

print("G3 to G4:", G3_G4)

print("G4 to G3:", G4_G3)
