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


### GAN

BUFFER_SIZE = 60000
BATCH_SIZE = 256
train_dataset = tf.data.Dataset.from_tensor_slices(x_train).shuffle(BUFFER_SIZE).batch(BATCH_SIZE)

def make_generator_model():
    model = tfk.Sequential()
    model.add(tkl.Dense(7*7*256, use_bias=False, input_shape=(100,)))
    model.add(tkl.BatchNormalization())
    model.add(tkl.LeakyReLU())

    model.add(tkl.Reshape((7, 7, 256)))
    assert model.output_shape == (None, 7, 7, 256) # Note: None is the batch size

    model.add(tkl.Conv2DTranspose(128, (5, 5), strides=(1, 1), padding='same', use_bias=False))
    assert model.output_shape == (None, 7, 7, 128)
    model.add(tkl.BatchNormalization())
    model.add(tkl.LeakyReLU())

    model.add(tkl.Conv2DTranspose(64, (5, 5), strides=(2, 2), padding='same', use_bias=False))
    assert model.output_shape == (None, 14, 14, 64)
    model.add(tkl.BatchNormalization())
    model.add(tkl.LeakyReLU())

    model.add(tkl.Conv2DTranspose(1, (5, 5), strides=(2, 2), padding='same', use_bias=False, activation='tanh'))
    assert model.output_shape == (None, 28, 28, 1)

    return model


renerator = make_generator_model()

noise = tf.random.normal([16, 100])
gen_img = generator(noise, training=False)

plot_imglist(gen_img)

def make_discriminator_model():
    model = tfk.Sequential()
    model.add(tkl.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
                                     input_shape=[28, 28, 1]))
    model.add(tkl.LeakyReLU())
    model.add(tkl.Dropout(0.3))

    model.add(tkl.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
    model.add(tkl.LeakyReLU())
    model.add(tkl.Dropout(0.3))

    model.add(tkl.Flatten())
    model.add(tkl.Dense(1))

    return model

#Testing model 

discriminator = make_discriminator_model()
decision = discriminator(gen_img)
print (decision)


#Cross entropy loss caclulation

cross_entropy = tfk.losses.BinaryCrossentropy(from_logits=True)


def discriminator_loss(real_output, fake_output):
    real_loss = cross_entropy(tf.ones_like(real_output), real_output)
    fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss

def generator_loss(fake_output):
    return cross_entropy(tf.ones_like(fake_output), fake_output)

#optimizing
generator_optimizer = tf.keras.optimizers.Adam(1e-4)
discriminator_optimizer = tf.keras.optimizers.Adam(1e-4)


#Strating conditions

EPOCHS = 10
noise_dim = 100
num_examples = 16


seed = tf.random.normal([num_examples, noise_dim])

@tf.function
def train_step(images):
    noise = tf.random.normal([BATCH_SIZE, noise_dim])

    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
        generated_images = generator(noise, training=True)

        real_output = discriminator(images, training=True)
        fake_output = discriminator(generated_images, training=True)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

    gradients_of_generator = gen_tape.gradient(gen_loss, generator.trainable_variables)
    gradients_of_discriminator = disc_tape.gradient(disc_loss, discriminator.trainable_variables)

    generator_optimizer.apply_gradients(zip(gradients_of_generator, generator.trainable_variables))
    discriminator_optimizer.apply_gradients(zip(gradients_of_discriminator, discriminator.trainable_variables))


##Training loops

def train(dataset, epochs):
    for epoch in tqdm(range(epochs)):

        for image_batch in dataset:
            train_step(image_batch)

        # Produce images for the GIF as we go
       # display.clear_output(wait=True)
        generate_and_save_images(generator,
                                 epoch + 1,
                                 seed)



        # Generate after the final epoch
        #display.clear_output(wait=True)
        #generate_and_save_images(generator,
        #                       epochs,
        #                      seed)

def generate_and_save_images(model, epoch, test_input):

    predictions = model(test_input, training=False)

    fig = plt.figure(figsize=(4,4))

    for i in range(predictions.shape[0]):
        plt.subplot(4, 4, i+1)
        plt.imshow(predictions[i, :, :, 0] * 127.5 + 127.5, cmap='gray')
        plt.axis('off')

    #plt.savefig('image_at_epoch_{:04d}.png'.format(epoch))
    plt.show()

train(train_dataset, 100)
