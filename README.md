# DenseNet-_-prediction-of-responders-and-non-responders-patients

L’objectif etant de prédire le pronostique des patients à partir d’un réseau de neurone à convolution entrainé sur les tuiles tumoral de patients répondeurs et non répondeur. Mais avant ça, il a fallu optimiser certaine étapes préliminaires à l’entrainement et qui sont la filtration des tuiles et la méthode de normalisation.

## Step by step 
## Filtration et normalisation des tuiles
Selection sur la base de l'ecart type des 5 premieres couleurs les plus dominantes sur ma tuile 
Il s'agit de l'utilisation d'un modéle de classification de type Kmean clustering 
Voici un exemple de clustering de couleurs sur une tuile 
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/tile_example.png)
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/kmean_clustering_tile.png)
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/kmean_clustering_scatter.png)

Normalisation par les deux méthodes : Reinhard et Macenko 




## our model architectur 
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/densenet_arch.png)
![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/densenet_arch_details.png)

