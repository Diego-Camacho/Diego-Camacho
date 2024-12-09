from __future__ import absolute_import, division, print_function, unicode_literals

import tensorflow as tf
from tensorflow.keras import datasets, layers, models
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import h5py
from sklearn.model_selection import train_test_split


device_name = tf.test.gpu_device_name()

print(device_name)


#Local analysis of hitopathology images
DATADIR = "/BaseCancer"
CATEGORIES = ["benign", "malignant"]
training_data = []
label_data=[]


def create_training_data():
  
    datasource = os.path.join(DATADIR, 'cancer_histo.hdf5')
    with h5py.File(datasource, 'r') as f:
      data = f['default']
      database = np.array(data)
    
    labessource = os.path.join(DATADIR, 'label_cancer.hdf5')
    with h5py.File(labessource, 'r') as f:
      data = f['default']
      labels = np.array(data)
    
    print(labels.shape)
      
    xTrain, xTest, yTrain, yTest = train_test_split(database, labels, test_size = 0.25, random_state = 0)
    xTrain = np.array(xTrain)
    xTest = np.array(xTest)
    yTrain = np.array(yTrain)
    yTest = np.array(yTest)
    return (xTrain, yTrain), (xTest, yTest)


def returnOneE(value):
  
    if value == 0:
      label = 'benign'
    else:
      label = 'malignat'
      
    return label
  

st = time.time()
(x,y),(z,w) = create_training_data()
elapsed_time = time.time() - st
elapsed_time_minutes = int(int(elapsed_time) / 60)
elapsed_time_seconds = int(elapsed_time) % 60
print('loaded in %s [min] %s [s]' % (elapsed_time_minutes, elapsed_time_seconds))

img_rev1 = x[0,:,:,:]
img_rev2 = x[1,:,:,:]
img_rev3 = x[2,:,:,:]
img_rev4 = x[3,:,:,:]

plt.figure(figsize=(20,10))

plt.subplot(1, 4, 1)
plt.imshow(img_rev1)
plt.title(returnOneE(y[0]))

plt.subplot(1, 4, 2)
plt.imshow(img_rev2)
plt.title(returnOneE(y[1]))

plt.subplot(1, 4, 3)
plt.imshow(img_rev3)
plt.title(returnOneE(y[3]))

plt.subplot(1, 4, 4)
plt.imshow(img_rev4)
plt.title(returnOneE(y[4]))

plt.show()



def normalize(x_train, x_test, y_train, y_test):

    x_train = np.array(x_train, dtype=np.float32)
    x_test = np.array(x_test, dtype=np.float32)
    y_train = np.array(y_train, dtype=np.float32)
    y_test = np.array(y_test, dtype=np.float32)


    (a, b, c, d) = x_train.shape
    for i in range(a):
        x_train[i, :, :, 0] /= 255
        x_train[i, :, :, 1] /= 255
        x_train[i, :, :, 2] /= 255

    (a, b, c, d) = x_test.shape
    for i in range(a):
        x_test[i, :, :, 0] /= 255
        x_test[i, :, :, 1] /= 255
        x_test[i, :, :, 2] /= 255

    return (x_train,y_train),(x_test,y_test)



def main():

    st = time.time()
    (train_images, train_labels), (test_images, test_labels) = create_training_data()


    print("loading data......")
    elapsed_time = time.time() - st
    elapsed_time_minutes = int(int(elapsed_time) / 60)
    elapsed_time_seconds = int(elapsed_time) % 60
    print('loaded in %s [min] %s [s]' % (elapsed_time_minutes, elapsed_time_seconds))


    # Normalize pixel values to be between 0 and 1
    st = time.time()
    (train_images, train_labels), (test_images, test_labels) = normalize(train_images,test_images,train_labels,test_labels)
    elapsed_time = time.time() - st
    elapsed_time_minutes = int(int(elapsed_time) / 60)
    elapsed_time_seconds = int(elapsed_time) % 60
    print('normalize data in %s [min] %s [s]' % (elapsed_time_minutes, elapsed_time_seconds))

    # Model CNN with dropout
    model = models.Sequential()
    #aqui cambiar numero de filtros donde dice 16,32,64, y en la capa densa poner el ultimo valor (64) y la ultima densa queda en 2s
    model.add(layers.Conv2D(16, (3, 3), padding='same', activation='relu', input_shape=(256,256,3)))
    model.add(layers.BatchNormalization())
    model.add(layers.MaxPooling2D((2, 2), strides=2))
    model.add(layers.Conv2D(32, (3, 3), padding='same', activation='relu'))
    model.add(layers.BatchNormalization())
    model.add(layers.MaxPooling2D((2, 2), strides=2))
    model.add(layers.Conv2D(64, (3, 3), padding='same', activation='relu'))
    model.add(layers.BatchNormalization())

    model.add(layers.Flatten())
    model.add(layers.Dropout(0.5))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.BatchNormalization())
    model.add(layers.Dense(2, activation='softmax'))

    model.summary()

    #training process
    model.compile(optimizer='sgd',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])
    
   

    #testing process
    history = model.fit(train_images, train_labels, epochs=5, batch_size=32, validation_data=(test_images, test_labels))
    print('\nhistory dict:', history.history)
    test_loss, test_acc = model.evaluate(test_images, test_labels)

    print(history)
    acc = history.history['acc']
    val_acc = history.history['val_acc']
    loss = history.history['loss']
    val_loss =  history.history['val_loss']

    print (acc)

    epochs = range(len(acc))

    plt.plot(epochs, acc, 'bo', label='Training acc')
    plt.plot(epochs, val_acc, 'b', label='Validation acc')
    plt.title('Training and validation accuracy')
    plt.legend()
    plt.grid()

    plt.figure()

    plt.plot(epochs, loss, 'bo', label='Training loss')
    plt.plot(epochs, val_loss, 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    
    plt.grid()
    plt.show()
    
    # Generate predictions (probabilities -- the output of the last layer)
    # on new data using `predict`
    print('\n# Generate predictions for 3 samples')
    predictions = model.predict(test_images[:3])
    print('predictions shape:', predictions.shape)
    print(predictions)
    
    
    
    
    # summarize filter shapes 
    for layer in model.layers:
	  # check for convolutional layer
      if 'conv' not in layer.name:
        continue
	    # get filter weights
      filters, biases = layer.get_weights()
      print(layer.name, filters.shape)
    



main()
