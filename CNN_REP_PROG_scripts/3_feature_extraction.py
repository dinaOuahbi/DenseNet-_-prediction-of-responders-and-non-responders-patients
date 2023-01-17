#!/usr/bin/env python3.7

"""
# ---------------------------------------------
# Programme: DenseNet.py
# Auteur DO
# Date 12/10/2022
#
# Extraction features a partir du modele qui est entraine a predire PROG REP
# -------------------------------------------------------------
"""

# On importe tous les modules dont on a besoin
import tensorflow as tf
from tensorflow.keras import layers, models
import keras
from tensorflow.keras.models import Model
import keras.backend as K
import numpy as np
import pandas as pd
import sys,os, math, pathlib, functools
from keras.preprocessing.image import ImageDataGenerator

arg = sys.argv[1]
model_opt = sys.argv[2]

#arg="TCGA-2J-AAB6-01Z-00-DX1.2FC4D66F-BFBB-48FA-AFCB-ABBC12F37E14"

# params
test_datagen = ImageDataGenerator()
SIZE=402
batch_size = 32

cp_dir='/work/shared/ptbc/CNN_Pancreas_V2/database/canceropol/NT_tiles/'

# Générateurs / model
cp_generator = test_datagen.flow_from_directory(directory=cp_dir, classes = [arg],target_size =(SIZE, SIZE),batch_size=batch_size,shuffle=False,seed=42)

if model_opt == 'model1':
    model = tf.keras.models.load_model('/work/shared/ptbc/CNN_Pancreas_V2/PROG_REP_results/CNN_REP_PROG_bscdata/results/macenko/V6/densenet/', compile = False) #densenet sur macenko
    output = '/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/Feature_extration_prog_rep/without_color_aug/'

#elif model_opt == 'model2':
#    model = tf.keras.models.load_model('/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/color_augmentation/densenet_V6_comb_datagen/densenet/', compile = False) # densenet combdatagen (dina + julie)
#    output = '/work/shared/ptbc/CNN_Pancreas_V2/results/TCGA/macenko_prog_rep/Feature_extration_prog_rep/densenet_v6_julieDatagen/'

elif model_opt == 'model3':
    model = tf.keras.models.load_model('/work/shared/ptbc/CNN_Pancreas_V2/PROG_REP_results/CNN_REP_PROG_bscdata/color_augmentation/densenet_V6/densenet/', compile = False) # densenet julie datagen
    output = '/work/shared/ptbc/CNN_Pancreas_V2/results/CP/macenko_prog_rep/Feature_extration_prog_rep/with_color_aug/'
else:
    print('MODEL IS REQUIRED !')
    sys.exit()

# On recupère les features qui caracterisent le plus les images (avant derniere couche)
intermediate_layer_model = Model(inputs = model.input, outputs = model.get_layer("global_average_pooling2d").output)
myName = cp_generator.filenames
intermediate_output = intermediate_layer_model.predict(cp_generator)
myPrediction = intermediate_output
myPrediction = pd.DataFrame(myPrediction, index = myName)
myPrediction.to_csv(output+"Features_"+arg+".csv", sep=',', encoding='utf-8',index=True, header = None)
