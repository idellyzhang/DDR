import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
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

pattern = 'feature_TCGA_BRCA.csv'
data_r = pd.read_csv(pattern)

data = biomarker_selection(data_r,4)

y = data .Group
X = data.drop('Group', axis=1)

test = 'validation_feature.csv'
data_test_r = pd.read_csv(test)
data_test = biomarker_selection(data_test_r,4)
y_test = data_test.Group
X_test = data_test.drop('Group', axis=1)


clf = svm.SVC(kernel='rbf', C=1, gamma='auto').fit(X, y)

y_pred = clf.predict(X_test)
cancer_pred = y_pred[0:5]
normal_pred = y_pred[5:19]
	    #print(MU_pred)
	    #print(MU_pred[MU_corr])
cancer_corr = cancer_pred == "TN"
cancer_normal = cancer_pred == "OT"

normal_corr = normal_pred == "OT"
normal_cancer = normal_pred == "TN"

	    #print(MU_pred[MU_corr]
cancer_accuracy = len(cancer_pred[cancer_corr])/5
cancer_normal = len(cancer_pred[cancer_normal])/5

normal_accuracy = len(normal_pred[normal_corr])/14
normal_cancer = len(normal_pred[normal_cancer])/14

acc = accuracy_score(y_test, y_pred)
rec = recall_score(y_test, y_pred, average='macro')
prec = precision_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')
#print(y_pred)
print("Overall accuracy:", acc)
print("Recall:", rec)
print("precision accuracy:", prec)
print("F1 Score:", f1)

print("TNBC accuracy:", cancer_accuracy)
print("OTHER accuracy:", normal_accuracy)

print("TNBC to OTHER:", cancer_normal)
print("OTHER to TNBC:", normal_cancer)

