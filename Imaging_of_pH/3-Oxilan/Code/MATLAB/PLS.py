


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from scipy.stats import randint as sp_randint
from sklearn.model_selection import RandomizedSearchCV
get_ipython().magic('matplotlib inline')


# In[2]:

zspectra = pd.read_csv('centered_cest.csv', header = None).values
conc = pd.read_csv('concentration.csv', header = None).values
pH = pd.read_csv('pH.csv', header = None).values


# In[5]:

def mymetric(yexp, ypred):
    d = np.sum((yexp - ypred)**2 )
    d = d / ypred.shape[0]
    d = np.sqrt(d)
    d = d / np.mean(yexp)
    d = 100 * d
    return d


# In[10]:

X = zspectra[:, 0:101:2]
Y = pH
num_components = 20
Error = np.zeros((num_components -1,1))

for idx,K in enumerate(np.arange(1,num_components)):
    X_train, X_test, y_train, y_test = train_test_split( X, Y, test_size=0.05, random_state=42)
    pls = PLSRegression(n_components = K, scale = False)
    pls.fit(X_train, y_train)
    y_hat = pls.predict(X_test)
    Error[idx] = mymetric(y_test , y_hat)

plt.plot( np.arange(1,num_components), Error ,'o-')


# In[12]:

steps = [1,4,8]
labels = list()

for step in steps:
    X = zspectra[:, 0:101:step]
    labels.append(int(X.shape[1]))
    Y = pH
    num_components = 10
    Error = np.zeros((num_components -1,1))

    for idx,K in enumerate(np.arange(1,num_components)):
        X_train, X_test, y_train, y_test = train_test_split( X, Y, test_size=0.50, random_state=42)
        pls = PLSRegression(n_components = K, scale = False)
        pls.fit(X_train, y_train)
        y_hat = pls.predict(X_test)
        Error[idx] = mymetric(y_test , y_hat)

    plt.plot( np.arange(1,num_components), Error ,'o-')
    plt.legend(labels)
