import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn import datasets
from sklearn import svm
from sklearn.metrics import accuracy_score, precision_score, recall_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.multiclass import OneVsRestClassifier
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

pattern = 'GSE37418_20_feature.csv'
data_r = pd.read_csv(pattern)

data = biomarker_selection(data_r,20)

WNT = data['Group'] == "WNT"
SHH = data['Group'] == "SHH"
G3 = data['Group'] == "G3"
G4 = data['Group'] == "G4"


y_WNT = data[WNT] .Group
X_WNT = data[WNT].drop('Group', axis=1)


y_SHH = data[SHH] .Group
X_SHH = data[SHH].drop('Group', axis=1)

y_G3 = data[G3] .Group
X_G3 = data[G3].drop('Group', axis=1)


y_G4 = data[G4] .Group
X_G4 = data[G4].drop('Group', axis=1)


accuracy = []
recall = []
WNT_acc = []
WNT_SHH_pred = []
WNT_G3_pred = []
WNT_G4_pred = []
SHH_acc = []
SHH_WNT_pred = []
SHH_G3_pred = []
SHH_G4_pred = []
G3_acc = []
G3_WNT_pred = []
G3_SHH_pred = []
G3_G4_pred = []
G4_acc = []
G4_WNT_pred = []
G4_SHH_pred = []
G4_G3_pred = []
for x in range(1000):
	    X_WNT_train, X_WNT_test, y_WNT_train, y_WNT_test = train_test_split(X_WNT, y_WNT,test_size=0.2)
	    X_SHH_train, X_SHH_test, y_SHH_train, y_SHH_test = train_test_split(X_SHH, y_SHH,test_size=0.2)
	    X_G3_train, X_G3_test, y_G3_train, y_G3_test = train_test_split(X_G3, y_G3,test_size=0.2)
	    X_G4_train, X_G4_test, y_G4_train, y_G4_test = train_test_split(X_G4, y_G4,test_size=0.2)
	    
	    y_train = pd.concat([y_WNT_train, y_SHH_train, y_G3_train, y_G4_train])
	    X_train = pd.concat([X_WNT_train, X_SHH_train, X_G3_train, X_G4_train])
	    y_test = pd.concat([y_WNT_test, y_SHH_test, y_G3_test, y_G4_test])
	    X_test = pd.concat([X_WNT_test, X_SHH_test, X_G3_test, X_G4_test])
	    clf = OneVsRestClassifier(svm.SVC(kernel='rbf', C=1, gamma='auto')).fit(X_train, y_train)

	    y_pred = clf.predict(X_test)

	    acc = accuracy_score(y_test, y_pred)

	   
	    rec = recall_score(y_test, y_pred, average='macro')
	    accuracy.append(acc)
	    recall.append(rec)
	    WNT_pred = y_pred[0:len(y_WNT_test)]
	    SHH_pred = y_pred[len(y_WNT_test):len(y_WNT_test)+len(y_SHH_test)]
	    G3_pred = y_pred[len(y_WNT_test)+len(y_SHH_test):len(y_WNT_test)+len(y_SHH_test)+len(y_G3_test)]
	    G4_pred = y_pred[len(y_WNT_test)+len(y_SHH_test)+len(y_G3_test):len(y_test)]
	    #print(MU_pred)
	    WNT_corr = WNT_pred == "WNT"
	    WNT_SHH = WNT_pred == "SHH"
	    WNT_G3 = WNT_pred == "G3"
	    WNT_G4 = WNT_pred == "G4"
	    #print(MU_pred[MU_corr])
	    SHH_corr = SHH_pred == "SHH"
	    SHH_WNT = SHH_pred == "WNT"
	    SHH_G3 = SHH_pred == "G3"
	    SHH_G4 = SHH_pred == "G4"
	    
	    G3_corr = G3_pred == "G3"
	    G3_WNT = G3_pred == "WNT"
	    G3_SHH = G3_pred == "SHH"
	    G3_G4 = G3_pred == "G4"
	    #print(MU_pred[MU_corr])
	    G4_corr = G4_pred == "G4"
	    G4_WNT = G4_pred == "WNT"
	    G4_SHH = G4_pred == "SHH"
	    G4_G3 = G4_pred == "G3"
	    
	    WNT_accuracy = len(WNT_pred[WNT_corr])/len(y_WNT_test)
	    WNT_SHH_p = len(WNT_pred[WNT_SHH])/len(y_WNT_test)
	    WNT_G3_p = len(WNT_pred[WNT_G3])/len(y_WNT_test)
	    WNT_G4_p = len(WNT_pred[WNT_G4])/len(y_WNT_test)
	    
	    SHH_accuracy = len(SHH_pred[SHH_corr])/len(y_SHH_test)
	    SHH_WNT_p = len(SHH_pred[SHH_WNT])/len(y_SHH_test)
	    SHH_G3_p = len(SHH_pred[SHH_G3])/len(y_SHH_test)
	    SHH_G4_p = len(SHH_pred[SHH_G4])/len(y_SHH_test)
	    
	    G3_accuracy = len(G3_pred[G3_corr])/len(y_G3_test)
	    G3_WNT_p = len(G3_pred[G3_WNT])/len(y_G3_test)
	    G3_SHH_p = len(G3_pred[G3_SHH])/len(y_G3_test)
	    G3_G4_p = len(G3_pred[G3_G4])/len(y_G3_test)
	    
	    G4_accuracy = len(G4_pred[G4_corr])/len(y_G4_test)
	    G4_WNT_p = len(G4_pred[G4_WNT])/len(y_G4_test)
	    G4_SHH_p = len(G4_pred[G4_SHH])/len(y_G4_test)
	    G4_G3_p = len(G4_pred[G4_G3])/len(y_G4_test)
	    
	    WNT_acc.append(WNT_accuracy)
	    WNT_SHH_pred.append(WNT_SHH_p)
	    WNT_G3_pred.append(WNT_G3_p)
	    WNT_G4_pred.append(WNT_G4_p)
	    
	    SHH_acc.append(SHH_accuracy)
	    SHH_WNT_pred.append(SHH_WNT_p)
	    SHH_G3_pred.append(SHH_G3_p)
	    SHH_G4_pred.append(SHH_G4_p)
	    
	    G3_acc.append(G3_accuracy)
	    G3_WNT_pred.append(G3_WNT_p)
	    G3_SHH_pred.append(G3_SHH_p)
	    G3_G4_pred.append(G3_G4_p)
	    
	    G4_acc.append(G4_accuracy)
	    G4_WNT_pred.append(G4_WNT_p)
	    G4_SHH_pred.append(G4_SHH_p)
	    G4_G3_pred.append(G4_G3_p)
	    #print(y_test)
print("Accuracy: ", mean(accuracy))
print("Recall: ", mean(recall))

print("WNT accuracy:", mean(WNT_acc))
print("SHH accuracy: ",mean(SHH_acc))
print("G3 accuracy:", mean(G3_acc))
print("G4 accuracy: ",mean(G4_acc))

print("WNT to SHH:", mean(WNT_SHH_pred))
print("WNT to G3:", mean(WNT_G3_pred))
print("WNT to G4:", mean(WNT_G4_pred))

print("SHH to WNT:", mean(SHH_WNT_pred))
print("SHH to G3:", mean(SHH_G3_pred))
print("SHH to G4:", mean(SHH_G4_pred))

print("G3 to WNT:", mean(G3_WNT_pred))
print("G3 to SHH:", mean(G3_SHH_pred))
print("G3 to G4:", mean(G3_G4_pred))

print("G4 to WNT:", mean(G4_WNT_pred))
print("G4 to SHH:", mean(G4_SHH_pred))
print("G4 to G3:", mean(G4_G3_pred))

