#####################################################################################################################################################################################################################
# Warnings
import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)
warnings.filterwarnings("ignore",category=FutureWarning)

# Data Processing
import pandas as pd
import numpy as np

# Modelling
from sklearn.ensemble import RandomForestClassifier
#from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
#from sklearn.model_selection import RandomizedSearchCV, train_test_split
#from scipy.stats import randint
#from sklearn.model_selection import cross_val_score
#from sklearn import svm
#from sklearn import linear_model
#from sklearn.tree import DecisionTreeClassifier
#from sklearn import svm

# Tree Visualisation
#from sklearn.tree import export_graphviz
#from IPython.display import Image
#import graphviz
#import matplotlib
#import matplotlib.pyplot as plt
#matplotlib.use('TkAgg') # or "TkAgg"
#import seaborn as sns
#from sklearn.metrics import roc_curve, roc_auc_score

# Model Preservation
import pickle

# System libraries
import os
import time
import sys
orig_stdout = sys.stdout
#####################################################################################################################################################################################################################
### load the model from disk
vector = str(sys.argv[1])
savefile = sys.argv[2]
possition = sys.argv[3]
clustermethod = sys.argv[4]

print(clustermethod)

if "RW" in clustermethod:
	filename = '/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/MachineLearning_Test/Scripts/20230720/RW/3_fix_cluster_evidence/20230720_model_proteorizer_randomforest_improved.sav'
	loaded_model = pickle.load(open(filename, 'rb'))
elif "HClust" in clustermethod:
	filename = '/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/MachineLearning_Test/Scripts/20230720/HClust/2_fix_cluster_evidence/20230720_model_proteorizer_randomforest_improved_HClust.sav'
	loaded_model = pickle.load(open(filename, 'rb'))
else:	### defaults to random walk
	filename = '/net/home.isilon/ag-russell/bq_tschmenger/PhD/MechismoScanner/PERTURBED_INTERFACES/EnrichmentProbability/Hereditary_Cancer/3D_Clustering_For_Any_Variant/MachineLearning_Test/Scripts/20230720/RW/3_fix_cluster_evidence/20230720_model_proteorizer_randomforest_improved.sav'
	loaded_model = pickle.load(open(filename, 'rb'))

vector = vector.replace("[","").replace("]","")
vect = []
for item in vector.split(","):
	vect.append(int(item))
prediction_result 	= loaded_model.predict(np.array(vect).reshape(1, -1))
prediction_proba 	= loaded_model.predict_proba(np.array(vect).reshape(1, -1))

outfile = open(savefile,"a")
sys.stdout = outfile
print(vect,"\t",str(prediction_result).replace("[","").replace("]",""),"\t",str(prediction_proba[0:1,0:1]).replace("[","").replace("]",""),"\t",str(prediction_proba[0:1,1:2]).replace("[","").replace("]",""),"\t",possition)		
sys.stdout = orig_stdout

