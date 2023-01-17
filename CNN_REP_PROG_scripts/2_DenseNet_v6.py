#!/usr/bin/env python3.7

"""
# ---------------------------------------------
# Programme: DenseNet.py
# Auteur DO
# Date 12/10/2022
#
# Construction dun CNN permettant la classification sur le canceroplol
# PROG (PFS < 6 mois) et REP (PFS >36 mois)
# -------------------------------------------------------------
"""

# On importe tous les modules dont on a besoin
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import datasets, layers, models
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D
from tensorflow.keras.layers import Flatten, Dense, Dropout
from tensorflow.keras.layers import Input
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers
from tensorflow.keras.utils import to_categorical
import keras.backend as K
import pandas as pd
from tensorflow.keras.callbacks import LearningRateScheduler
from tensorflow.keras.callbacks import ReduceLROnPlateau
from tensorflow.keras.applications import DenseNet121
from sklearn.utils import class_weight
from keras.utils import np_utils
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from collections import Counter
import time, os, sys, math, pathlib, functools
start = time.time()
import numpy as np
from skimage import io
from PIL import Image
from tensorflow.keras.preprocessing import image_dataset_from_directory
from tensorflow.keras.layers import Flatten, Dense, Dropout, GlobalAveragePooling2D

sys.path.append("/work/shared/ptbc/PETACC08_JFE/ColoClass_DL/Analyses_stats/Programmes/package_generator/")
import ImageDataGenerator_perso as IDG

print("\n\n+++++ Running CNN DenseNet+++++\n\n")

##### Tuto : CNN with Tensorflow|Keras for Fashion MNIST dans Kaggle et TensorFlow CNN, Data Augmentation: Prostate Cancer pour l'augmentation des données
##### https://penseeartificielle.fr/tp-reseau-de-neurones-convolutifs/

# define some params
SIZE = 402
batch_size = 64#32
EPOCHS=10#20
myData = '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/besancon_macenko_DB_survie/'
output = '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/color_augmentation/densenet_V6_comb_datagen/'

# remove augmented data
def remove_aug_data():
    for data in ['Train', 'Val']:
        root = f'{myData}{data}/REP/'
        for file in os.listdir(root):
            if file.startswith('aug'):
                #print(file)
                os.remove(f'{root}{file}')

# AUGMENTÉ LA CLASSE MINORITAIRE (2 fois)
mino_datagen = ImageDataGenerator(
        rotation_range=45,     #Random rotation between 0 and 45
        width_shift_range=0.2,   #% shift
        height_shift_range=0.2,
        shear_range=0.2,
        zoom_range=0.2,
        horizontal_flip=True,
        fill_mode='constant', cval=0)    #Also try nearest, constant, reflect, wrap

##### FOR TRAIN AND VAL
def minority_class_augmentation(factor):
    for data in ['Train', 'Val']:
        image_directory = f'{myData}{data}/REP/'
        dataset = []
        for i, image_name in enumerate(os.listdir(image_directory)):
            if (image_name.split('.')[1] == 'tif'):
                image = io.imread(image_directory + image_name)
                image = Image.fromarray(image, 'RGB')
                image = image.resize((SIZE,SIZE))
                dataset.append(np.array(image))
        x = np.array(dataset)

        print(f'----------------> {data}_REP = {x.shape}<----------------------')
        i = 1 # une img donnera une aug
        for batch in mino_datagen.flow(x, batch_size=x.shape[0],save_to_dir=image_directory,save_prefix='aug',save_format='jpg'):
            i += 1
            if i > factor: # une img donnera une aug
                break  # otherwise the generator would loop indefinitely

# check data SIZE
#train test val destribution
def check_data_SIZE(myData):
    for fol in os.listdir(myData):
        for sub_fol in os.listdir(f'{myData}{fol}'):
            print(f'--> {fol}/{sub_fol}')
            print(len(os.listdir(f'{myData}{fol}/{sub_fol}')))

# MAIN
remove_aug_data()
minority_class_augmentation(1)
check_data_SIZE(myData)

# On charge les donnees
test_datagen = ImageDataGenerator()

'''train_datagen = IDG.ImageDataGenerator_perso(
    rescale = 1.0 / 255,
    rotation_range = 40,
    width_shift_range = 0.2,
    height_shift_range = 0.2,
    shear_range = 0.2,
    zoom_range = 0.2,
    horizontal_flip = True,
    fill_mode = "nearest",
    color_augmentation = True
)'''

# combinaison ImageDataGenerator_perso et ImageDataGenerator_perso
train_datagen = IDG.ImageDataGenerator_perso(
    rescale = 1.0 / 255,
    rotation_range = 20,
    width_shift_range = 0.05,
    height_shift_range = 0.05,
    shear_range = 0.05,
    zoom_range = 0.05,
    horizontal_flip = True,
    fill_mode = "nearest",
    color_augmentation = True
)

'''train_datagen = ImageDataGenerator(
        rotation_range=20,
        zoom_range=0.05,
        width_shift_range=0.05,
        height_shift_range=0.05,
        shear_range=0.05,
        horizontal_flip=True,
        fill_mode="nearest"
)'''

train_dir = f"{myData}Train"
val_dir = f"{myData}Val"
test_dir = f"{myData}Test"

# Générateurs
train_generator = train_datagen.flow_from_directory(
    directory=train_dir,
    target_size=(402, 402),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=True,
    seed=42, subset = 'training'
)
val_generator = train_datagen.flow_from_directory(
    directory=val_dir,
    target_size=(402, 402),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=True,
    seed=42
)
test_generator = test_datagen.flow_from_directory(
    directory=test_dir,
    target_size=(402, 402),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=True,
    seed=42
)

# On calcule le poids de chacune des deux classes car elles ne sont pas donnees avec le meme poids dans l entrainement
counter = Counter(train_generator.classes)
max_val = float(max(counter.values()))
class_weights = {class_id : max_val/num_images for class_id, num_images in counter.items()}

# instanciez un modèle de base avec des poids pré-entraînés.
base_model = DenseNet121(input_shape= (SIZE,SIZE,3),
                                      include_top=False, # Do not include the ImageNet classifier at the top.
                                      weights='imagenet')# Load weights pre-trained on ImageNet.

#Ensuite, gelez le modèle de base.
base_model.trainable = False
preprocess_input = tf.keras.applications.densenet.preprocess_input

# MY MODEL
inputs = preprocess_input(tf.keras.Input(shape=(SIZE, SIZE, 3)))
# We make sure that the base_model is running in inference mode here, by passing `training=False`. This is important for fine-tuning, as you will learn in a few paragraphs.
x = base_model(inputs, training=False)
# Convert features of shape `base_model.output_shape[1:]` to vectors
x = GlobalAveragePooling2D()(x)
#x = keras.layers.Dropout(0.5)(x)
# A Dense classifier with a single unit (binary classification)
outputs=Dense(2,activation='softmax',kernel_regularizer=tf.keras.regularizers.L1(0.01), activity_regularizer=tf.keras.regularizers.L2(0.01))(x)
model = Model(inputs, outputs)

# COMPILE
model.compile(loss=tf.keras.losses.CategoricalCrossentropy(from_logits=True),optimizer=tf.keras.optimizers.Adam(),metrics = [keras.metrics.BinaryAccuracy()])

# Entrainement
model.fit(train_generator,epochs=EPOCHS, class_weight = class_weights,validation_data=val_generator, verbose = 2)


#Une fois que votre modèle a convergé sur les nouvelles données,
#vous pouvez essayer de dégeler tout ou partie du modèle de base et
#recycler l'ensemble du modèle de bout en bout avec un taux d'apprentissage très faible.

# Unfreeze the base model
base_model.trainable = True
model.compile(optimizer=keras.optimizers.Adam(1e-5),  # Very low learning rate
              loss=keras.losses.BinaryCrossentropy(from_logits=True),
              metrics=[keras.metrics.BinaryAccuracy()])
# Entrainement
history=model.fit(train_generator,epochs=EPOCHS, class_weight = class_weights,validation_data=val_generator, verbose = 2)

# On enregistre les poids du modele et le modele
model.save(f'{output}densenet')
model.save_weights(f'{output}densenet.h5')
hist_df = pd.DataFrame(history.history)
hist_csv_file = f'{output}prog_rep_history.csv'
with open(hist_csv_file, mode='w') as f:
    hist_df.to_csv(f)

remove_aug_data()
check_data_SIZE(myData)

# Predictions
train_generator = test_datagen.flow_from_directory(
    directory=train_dir,
    target_size=(SIZE, SIZE),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=False,
    seed=42
)
val_generator = test_datagen.flow_from_directory(
    directory=val_dir,
    target_size=(SIZE, SIZE),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=False,
    seed=42
)
test_generator = test_datagen.flow_from_directory(
    directory=test_dir,
    target_size=(SIZE, SIZE),
    color_mode="rgb",
    batch_size=batch_size,
    class_mode='categorical',
    shuffle=False,
    seed=42
)

myName = train_generator.filenames
test_predictions_baseline = model.predict(train_generator)
myPrediction = test_predictions_baseline
myPrediction = pd.DataFrame(myPrediction, index = myName)
myPrediction.to_csv(output+"myPrediction_train.csv", sep=',', encoding='utf-8',
               index=True, header = None)

myName = val_generator.filenames
test_predictions_baseline = model.predict(val_generator)
myPrediction = test_predictions_baseline
myPrediction = pd.DataFrame(myPrediction, index = myName)
myPrediction.to_csv(output+"myPrediction_val.csv", sep=',', encoding='utf-8',
               index=True, header = None)

myName = test_generator.filenames
test_predictions_baseline = model.predict(test_generator)
myPrediction = test_predictions_baseline
myPrediction = pd.DataFrame(myPrediction, index = myName)
myPrediction.to_csv(output+"myPrediction_test.csv", sep=',', encoding='utf-8',
               index=True, header = None)


print(f'-'*100)
print('DONE ! ')
