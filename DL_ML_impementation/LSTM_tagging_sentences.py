import tensorflow.keras as keras
import nltk
nltk.download("treebank")
nltk.download("universal_tagset")

from nltk.corpus import treebank

sentences = treebank.tagged_sents(tagset="universal")
print(sentences[0])

sentences_data = list()
tags_data = list()

#word2vec for vectorizing

for sentence in sentences:
    s, t = zip(*sentence)
    sentences_data.append(s)
    tags_data.append(t)

from sklearn.model_selection import train_test_split

#Dataset divided in 3 sets: train, dev and test

train_sent, test_sent, train_tags, test_tags = train_test_split(sentences_data, tags_data, test_size=0.2)
train_sent, dev_sent, train_tags, dev_tags = train_test_split(train_sent, train_tags, test_size=0.1)


# Initialize a list of tags and words seen in the training set

uniq_tags = list()
uniq_words = list()

#Training data

for sent,tags in zip(train_sent, train_tags):
    for word, tag in zip(sent, tags):
        uniq_tags.append(tag)
        uniq_words.append(word)

uniq_tags = list(set(uniq_tags))
uniq_words = list(set(uniq_words))

word_to_indx = {"PAD":0,"OOV":1, "BOS":2, "EOS":3}
indx_to_word = {0:"PAD",1:"OOV", 2:"BOS", 3:"EOS"}


tag_to_indx = {"PAD":0, "OOV":1, "BOS":2, "EOS":3}
indx_to_tag = {0:"PAD", 1:"OOV", 2:"BOS", 3:"EOS"}

for line in train_sent:
    for word in line:
        if not word in word_to_indx.keys():
            word_to_indx[word] = len(word_to_indx)
            indx_to_word[len(word_to_indx)-1] = word
    
for line in train_tags:
    for tag in line:
        if not tag in tag_to_indx.keys():
            tag_to_indx[tag] = len(tag_to_indx)
            indx_to_tag[len(tag_to_indx)-1] = tag

def sent_to_int(sent):
    int_sent = list()
    for word in sent:
        if word in word_to_indx.keys():
            int_sent.append(word_to_indx[word])
        else:
            int_sent.append(word_to_indx["OOV"])
    return int_sent

def tag_to_int(sent):
    int_sent = list()
    for tag in sent:
        if tag in tag_to_indx.keys():
            int_sent.append(tag_to_indx[tag])
        else:
            int_sent.append(tag_to_indx["OOV"])
    return int_sent

#Cast the setences into values

train_X, dev_X, train_Y, dev_Y, test_X, test_Y = list(), list(), list(), list(), list(), list()

for lineX, lineY in zip(train_sent, train_tags):
    train_X.append(sent_to_int(lineX))
    train_Y.append(tag_to_int(lineY))

for lineX, lineY in zip(dev_sent, dev_tags):
    dev_X.append(sent_to_int(lineX))
    dev_Y.append(tag_to_int(lineY))

for lineX, lineY in zip(test_sent, test_tags):
    test_X.append(sent_to_int(lineX))
    test_Y.append(tag_to_int(lineY))


from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.utils import to_categorical
#Model and training


MAX_LENGTH = len(max(train_X, key=len))
train_X = pad_sequences(maxlen=MAX_LENGTH, sequences=train_X, padding="post")
train_Y = pad_sequences(maxlen=MAX_LENGTH, sequences=train_Y, padding="post")
dev_X = pad_sequences(maxlen=MAX_LENGTH, sequences=dev_X, padding="post")
dev_Y = pad_sequences(maxlen=MAX_LENGTH, sequences=dev_Y, padding="post")
test_X = pad_sequences(maxlen=MAX_LENGTH, sequences=test_X, padding="post")
test_Y = pad_sequences(maxlen=MAX_LENGTH, sequences=test_Y, padding="post")

print(len(train_Y))
print(len(train_X[0]))
print(train_X[0])


train_cat_Y = [to_categorical(i, num_classes=len(tag_to_indx)) for i in train_Y]
dev_cat_Y = [to_categorical(i, num_classes=len(tag_to_indx)) for i in dev_Y]
test_cat_Y = [to_categorical(i, num_classes=len(tag_to_indx)) for i in test_Y]

from tensorflow.keras.models import Model
from tensorflow.keras.layers import LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional


input = keras.Input(shape=(MAX_LENGTH,))
model = Embedding(input_dim=len(word_to_indx), output_dim=50, input_length=MAX_LENGTH)(input)
model = Dropout(0.1)(model)
model = Bidirectional(LSTM(units=100, return_sequences=True, recurrent_dropout=0.1))(model)
out = TimeDistributed(Dense(len(tag_to_indx), activation="softmax"))(model)  # softmax output layer

model = Model(input, out)
model.compile(optimizer="rmsprop", loss="categorical_crossentropy", metrics=["accuracy"])
history = model.fit(train_X, np.array(train_cat_Y), batch_size=32, epochs=5, validation_data=[dev_X, np.array(dev_cat_Y)], verbose=1)


#SAving trained model

model.save("Models/pos_tagging.h5")



