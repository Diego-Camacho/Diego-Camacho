#############################
#### ML and DL exercises ####
#### Summer 2020         ####
#############################

import os
import sys
from tqdm.autonotebook import tqdm
# scientific python stack
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ML/DL
import sklearn
import sklearn.model_selection
from sklearn.model_selection import train_test_split

import tensorflow as tf
import tensorflow.keras as tfk
import tensorflow.keras.layers as tkl
print('Tensorflow:{}'.format(tf.__version__))
print('Keras:{}'.format(tfk.__version__))

etiquetas = ["T-shirt/top",
             "Trouser",    
             "Pullover",   
             "Dress",      
             "Coat",       
             "Sandal",      
             "Shirt",       
             "Sneaker",     
             "Bag",         
             "Ankle boot"]

from tensorflow.keras.datasets import fashion_mnist

(x_train, y_train), (x_test, y_test) = fashion_mnist.load_data()

x_train = x_train/255.0 * - 0.5
x_test = x_test/255.0 * - 0.5
x_train = x_train[:,:,:,np.newaxis].astype(np.float32)
x_test = x_test[:,:,:,np.newaxis].astype(np.float32)

def plot_img(img):
    plt.imshow(img[:,:,0]* 127.5 + 127.5, cmap='Greys_r')
    plt.axis('off')
    plt.show()
def plot_imglist(imglist):
    n=imglist.shape[0]
    for i in range(n):
        plt.subplot(4,max(4,n/4), i+1)
        plt.imshow(imglist[i,:,:,0]* 127.5 + 127.5, cmap='Greys_r')
        plt.axis('off')
    plt.show()

    
plot_img(x_train[0])
plot_imglist(x_train[:80])

#Model implamantation and fiting

model = tfk.Sequential([
    tkl.Flatten(input_shape=(28, 28,1)),
    tkl.Dense(128, activation='selu'),
    tkl.Dense(10, activation='softmax')
])
model.summary()
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
model.fit(x_train,y_train, epochs=15)

#Visualization of classification

emb_model = tfk.Model(inputs=model.inputs, outputs=model.layers[-2].output)
emb = emb_model.predict(x_test)
print(emb.shape)

###Using now a CNN model

model = tfk.Sequential()
model.add(tkl.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
                                 input_shape=[28, 28, 1]))
model.add(tkl.LeakyReLU())
model.add(tkl.Dropout(0.3))
model.add(tkl.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
model.add(tkl.LeakyReLU())
model.add(tkl.Dropout(0.3))
model.add(tkl.Conv2D(128, (3, 3), strides=(2, 2), padding='same'))
model.add(tkl.LeakyReLU())
model.add(tkl.Dropout(0.3))
model.add(tkl.Flatten())
model.add(tkl.Dense(10, activation='softmax'))

model.summary()
model.compile(optimizer='adam',
              loss='sparse_categorical_crossentropy',
              metrics=['accuracy'])
model.fit(x_train,y_train, epochs=10)


emb_model = tfk.Model(inputs=model.inputs, outputs=model.layers[-2].output)
emb = emb_model.predict(x_test)
print(emb.shape)
from sklearn.decomposition import PCA
pca_model = PCA(2)
x_pca=pca_model.fit_transform(emb)
plt.scatter(x_pca[:,0],x_pca[:,1],s=1,c=y_test,cmap='Dark2')
plt.legend()
plt.show()


