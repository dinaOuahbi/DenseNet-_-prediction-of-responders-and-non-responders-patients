'''
pour ce script j'ai besoin de :
    models
    tile xy
    scan reinhard (tuile normaliser pour une lame)
ce script prend un ID lame et renvoie un jpeg grad_cam
le output est a modifier
'''


import numpy as np
import sys
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.preprocessing.image import load_img
from tensorflow.keras.utils import img_to_array, array_to_img, save_img
from IPython.display import Image, display
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tensorflow.keras.layers import GlobalAveragePooling2D, Dense
from tensorflow.keras.applications import * #Efficient Net included here
from tensorflow.keras.applications import DenseNet121
import pandas as pd
import os
from PIL import Image

slide = sys.argv[1] #21_0_005_00_199

pfs_features = [151,180,213,423,572,611,648,959,1002]
def modify_tensor(my_grad, feature_list):
    my_l = [0 for i in range(my_grad.shape[0])]
    for i in feature_list:
        my_l[i] = my_grad[i]
    new_tensor = tf.constant(np.array(my_l), dtype=tf.float32)
    return new_tensor

def build_model_gradcam(size):
    #base = EfficientNetB3(weights="imagenet",include_top=False, input_shape=(402,402,3)) # transfert learning
    base = DenseNet121(input_shape= (402,402,3),include_top=False,weights='imagenet')
    x= GlobalAveragePooling2D()(base.output)
    x=Dense(2,activation='softmax')(x)
    model = tf.keras.Model(inputs=base.input, outputs=x)
    opt = tf.optimizers.SGD(learning_rate=1e-3)
    loss = tf.keras.losses.BinaryCrossentropy()
    model.compile(optimizer=opt,loss=loss,metrics=['AUC'])
    return model


# Fonctions venant du tuto : https://keras.io/examples/vision/grad_cam/
def make_gradcam_heatmap(img_array, model, last_conv_layer_name, pred_index=None):
    # First, we create a model that maps the input image to the activations
    # of the last conv layer as well as the output predictions
    grad_model = tf.keras.models.Model(
        [model.inputs], [model.get_layer(last_conv_layer_name).output, model.output]
    )
    # Then, we compute the gradient of the top predicted class for our input image
    # with respect to the activations of the last conv layer
    with tf.GradientTape() as tape:
        last_conv_layer_output, preds = grad_model(np.expand_dims(img_array, axis=0))
        if pred_index is None:
            pred_index = tf.argmax(preds[0])
        class_channel = preds[:, pred_index]
    # This is the gradient of the output neuron (top predicted or chosen)
    # with regard to the output feature map of the last conv layer
    grads = tape.gradient(class_channel, last_conv_layer_output)
    # This is a vector where each entry is the mean intensity of the gradient
    # over a specific feature map channel
    pooled_grads = tf.reduce_mean(grads, axis=(0, 1,2))
    #pooled_grads = modify_tensor(pooled_grads, pfs_features)

    # We multiply each channel in the feature map array
    # by "how important this channel is" with regard to the top predicted class
    # then sum all the channels to obtain the heatmap class activation
    last_conv_layer_output = last_conv_layer_output[0]
    heatmap = last_conv_layer_output @ pooled_grads[..., tf.newaxis]
    heatmap = tf.squeeze(heatmap)
    # For visualization purpose, we will also normalize the heatmap between 0 & 1
    heatmap = tf.maximum(heatmap, 0) / tf.math.reduce_max(heatmap)
    return heatmap.numpy(), preds.numpy()[0]


def display_gradcam(img, heatmap,alpha):
    # Load the original image
    img = keras.preprocessing.image.img_to_array(img)
    # Rescale heatmap to a range 0-255
    heatmap = np.uint8(255 * heatmap)
    # Use jet colormap to colorize heatmap
    jet = cm.get_cmap("Set1")
    # Use RGB values of the colormap
    jet_colors = jet(np.arange(256))[:, :3]
    jet_heatmap = jet_colors[heatmap]
    # Create an image with RGB colorized heatmap
    jet_heatmap = array_to_img(jet_heatmap)
    jet_heatmap = jet_heatmap.resize((img.shape[1], img.shape[0]))
    jet_heatmap = img_to_array(jet_heatmap)
    # Superimpose the heatmap on original image
    superimposed_img = jet_heatmap * alpha + img * (1 - alpha)
    superimposed_img = array_to_img(superimposed_img)
    #keras.preprocessing.image.save_img(output_file,superimposed_img)
    return superimposed_img


if __name__ == "__main__":
    print('==========> CONSTRUIRE LE MODEL GRAD CAM ')
    #model = tf.keras.models.load_model('../model/NvsT_EfficientNetB3.h5', compile = False)
    model = tf.keras.models.load_model('/work/shared/ptbc/CNN_Pancreas_V2/Analyses_stat/Resultats/StromaVSautre_DB10/NvsT_densenet/', compile = False)
    model_gradcam = build_model_gradcam(402)
    model_gradcam.set_weights(model.get_weights()) #  et on lui donne les poids du CNN2 que l'on a entraine prealablement
    last_conv_layer_name = list(filter(lambda x: isinstance(x, keras.layers.Activation), model_gradcam.layers))[-1].name
    print(f"Ma couche d'interet de GRADCAM s'appelle {last_conv_layer_name}")
    print(f'Le shape de ma derniere couche de model_gradcam est : {model_gradcam.get_layer(last_conv_layer_name).output.shape}')

    print('==========> CONSTRUIRE LA TABLE DES TUILES ET LES COORDONNÉES')
    coordonate = f'/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/tiles/{slide}/'
    images = '/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/SCAN_reinhard/'
    xy = pd.read_csv(f'{coordonate}{slide}_tileXY.txt', sep='\t')
    slide_list = os.listdir(f'{images}{slide}')
    xy.reset_index(inplace=True)
    xy['point'] = [i.split(' ')[1] for i in xy['Tile-Point']]
    xy['tile'] = [i.split(' ')[0] for i in xy['Tile-Point']]
    xy = xy[xy['point'] == 'Point1']
    xy.drop(['Tile-Point','point','Nb points'], axis=1, inplace=True)
    xy['index'] = [f'{j}.tif' for j in [i.replace('Tile', f'{slide}_') for i in xy['tile']]]
    slide_list = pd.DataFrame(slide_list).rename(columns={0:'index'})
    final_df = pd.merge(slide_list, xy, on='index').drop('tile', axis=1)
    print(f'shape de la TABLE DES TUILES ET LES COORDONNÉES : {final_df.shape}')

    print('==========> CONSTRUIRE UNE MATRICE RGB DE DIM : MAX(X) X MAX(Y)')
    my_img = Image.new('RGB', (max(final_df['x']), max(final_df['y'])))

    print('==========> PARCOURIR LES TUILES POUR GENERER LES GRAD CAM ')
    for i, file in enumerate(final_df['index']):
        print(f'Image {i} / {final_df.shape[0]}')
        img_path=f"{images}{slide}/{file}"
        img = keras.preprocessing.image.load_img(img_path)
        heatmap, preds = make_gradcam_heatmap(keras.preprocessing.image.img_to_array(img)/255,model_gradcam, last_conv_layer_name, pred_index = 1) #0 : normal / 1 : tumeur
        print(preds)
        superimposed_img = display_gradcam(img, heatmap, 0.5)#, output_file = f"../Resultats/CAM/{img_path.split('/')[-1]}")  on superpose l'image de la tuile et de la heatmap et on lenregistre
        # si le stroma est predit comme de la tumeur:
        #if float(preds[1]) <= 0.98:
        #    superimposed_img = Image.new(mode = "RGB", size = (402, 402),color = (176,224,230)) #GRAY
        my_img.paste(superimposed_img,(final_df.set_index('index').loc[file, 'x'], final_df.set_index('index').loc[file, 'y']))
    final_img = my_img.resize((5000,5000))
    final_img.save(f'/work/shared/ptbc/CNN_Pancreas_V2/Donnees/CNN_REP_PROG_bscdata/results/WSI_GRAD_CAM/{slide}.jpeg')
