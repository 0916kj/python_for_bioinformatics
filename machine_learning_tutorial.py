"""
ML Practice
https://medium.freecodecamp.org/the-hitchhikers-guide-to-
machine-learning-algorithms-in-python-bfad66adb378

Practice on ML basics with Scikit-Learn
Regressions, Decision Trees, Support Vector Machines, K-Nearest Neighbors,
Visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pydotplus
import collections
from pandas.plotting import parallel_coordinates

from sklearn import linear_model
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn import tree
from sklearn.cross_validation import train_test_split
from sklearn.datasets import load_iris
from sklearn import tree
from sklearn import svm

url = "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"
names = ['sepal-length', 'sepal-width', 'petal-length', 'petal-width', 'class']
dataset = pd.read_csv(url, names=names)
##import Iris dataset, apply names, convert to DataFrame

### scatter plot matrix
pd.plotting.scatter_matrix(dataset)
plt.show()
### see which variables have linear correlations
### data is stored as five Series; extract 2 Series to make a linear model

### plot points and linear regression
sns.set_context('notebook', font_scale=1.5)
sns.set_style('ticks')
sns.lmplot('sepal-length','petal-length',data=dataset,hue="class")
### make "class" (3 different flower species) a subset, visualize
### linear regression on those 3 classes
plt.ylabel('Response')
plt.xlabel('Explanatory')

### Now implement ML linear regression ###
linear = linear_model.LinearRegression()
trainX = np.asarray(dataset["sepal-length"][20:len(dataset["sepal-length"])]).reshape(-1, 1)
trainY = np.asarray(dataset["petal-length"][20:len(dataset["petal-length"])]).reshape(-1, 1)
testX = np.asarray(dataset["sepal-length"][:20]).reshape(-1, 1)
testY = np.asarray(dataset["petal-length"][:20]).reshape(-1, 1)
linear.fit(trainX, trainY)
linear.score(trainX, trainY)
print('Coefficient: \n', linear.coef_)
print('Intercept: \n', linear.intercept_)
print('R² Value: \n', linear.score(trainX, trainY))
predicted = linear.predict(testX)
### Logistic Regression ###
## Supervised classification algorithm - predicts values between 0 and 1
##use petal-length/10 as outcome variable
### Visualization
dataset["reg_peta_length"] = (dataset["petal-length"]/10).round(0).astype(int)
### convert outcome variable to integer so logistic regression is possible

sns.set_context("notebook", font_scale=1.1)
sns.set_style("ticks")
sns.lmplot('sepal-length',"reg_peta_length",data=dataset,logistic=True)
plt.ylabel("Probability")
plt.xlabel("Explanatory")

### Implementation
logistic = LogisticRegression()
X = (np.asarray(dataset['sepal-length'])).reshape(-1, 1)
Y = (np.asarray(dataset["reg_peta_length"])).ravel()
logistic.fit(X, Y)
logistic.score(X, Y)
print('Coefficient: \n', logistic.coef_)
print('Intercept: \n', logistic.intercept_)
print('R² Value: \n', logistic.score(X, Y))

### Decision Trees ###
decision = tree.DecisionTreeClassifier(criterion="gini")
X = dataset.values[:, 0:4] #all numerical parameters (array)
Y = dataset.values[:, 4] #classes (list)
trainX, testX, trainY, testY = train_test_split( X, Y, test_size = 0.3)
decision.fit(trainX, trainY)
print("Accuracy: \n", decision.score(testX, testY))

### Visualization
### Training
X = dataset.values[:, 0:4] #all numerical parameters (array)
Y = dataset.values[:, 4] #classes (list)
clf = tree.DecisionTreeClassifier()
clf = clf.fit(X,Y)
### Visualize data ###
dot_data = tree.export_graphviz(clf,
                                out_file=None,
                                feature_names=dataset.columns[0:4],  
                         class_names=list(set(Y)),
                                filled=True,
                                rounded=True)
graph = pydotplus.graph_from_dot_data(dot_data)

colors = ('turquoise', 'orange')
edges = collections.defaultdict(list)

for edge in graph.get_edge_list():
    edges[edge.get_source()].append(int(edge.get_destination()))

for edge in edges:
    edges[edge].sort()    
    for i in range(2):
        dest = graph.get_node(str(edges[edge][i]))[0]
        dest.set_fillcolor(colors[i])

graph.write_png('tree.png') ##visualization saved to directory
#### Support Vector Machines (SVM) ####
### Training 
### Use 2 variables and class
support = svm.SVC(kernel='linear', C = 1.0)
X = dataset.values[:, 0:2] #all numerical parameters (array)
Y = dataset.values[:, 4] #classes (list)
trainX, testX, trainY, testY = train_test_split( X, Y, test_size = 0.3)
support.fit(trainX, trainY)
print("Accuracy: \n", support.score(testX, testY))
pred = support.predict(testX)
### Scatter plot: X1 and X2 plotted, classes by color
sns.set_context("notebook", font_scale=1.1)
sns.set_style("ticks")
sns.lmplot("sepal-length","sepal-width", scatter=True, 
           fit_reg=False, data=dataset, hue="class")
plt.ylabel("sepal-length")
plt.xlabel("sepal-width")

### K-Nearest Neighbors (KNN)
### Implementation
#### Use same training values as SVM example
neighbors = KNeighborsClassifier(n_neighbors=5)

trainX, testX, trainY, testY = train_test_split(X, Y, test_size = 0.3)
neighbors.fit(trainX, trainY)
print("Accuracy: \n", neighbors.score(testX, testY))
pred = neighbors.predict(testX)
### Visualization ###
"""Parallel coordinates is a plotting technique for plotting multivariate data. 
It allows one to see clusters in data and to estimate other statistics visually. 
Using parallel coordinates points are represented as connected line segments. 
Each vertical line represents one attribute. One set of connected line 
represents one data point. Points that tend to cluster will appear 
closer together."""
plt.figure(figsize=(15,10))
parallel_coordinates(dataset,"class")
plt.title('Parallel Coordinates Plot', fontsize=20, fontweight='bold')
plt.xlabel('Features', fontsize=15)
plt.ylabel('Features values', fontsize=15)
plt.legend(loc=1, prop={'size': 15}, frameon=True,shadow=True, 
           facecolor="white", edgecolor="black")
plt.show()
