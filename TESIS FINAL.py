# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 03:42:04 2020

@author: abe20
"""

from scipy import integrate
import numpy as np
import sklearn.datasets
import pandas as pd
import umap
import matplotlib.pyplot as plt
import umap.plot
#conda install -c confa-forge hdbscan
import hdbscan 
import sklearn.cluster as cluster 
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource, CategoricalColorMapper
from bokeh.palettes import Spectral10
from bokeh.models import ColorBar, ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure, output_file, show
from bokeh.transform import linear_cmap
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar
from bokeh.io import output_file, show
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource,
                          LinearColorMapper, PrintfTickFormatter,)
from bokeh.palettes import Plasma
from bokeh.palettes import Category10
from bokeh.palettes import Viridis
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from bokeh.plotting import figure, output_file, save

groups=pd.read_excel(r'C:\Users\abe20\OneDrive\Documentos\Python Scripts\groupshist.xlsx')
groups=groups.fillna("NA/NA")


fitval=pd.read_excel(r'C:\Users\abe20\OneDrive\Documentos\Python Scripts\Datmulnumnum.xlsx')
fitval.replace([np.inf,-np.inf],np.nan)
fitval.dropna(inplace=True)
fitval.astype("int")
fitvalT=fitval.transpose()
classes =list (range(530))



#3
#el primer analisi labels normal es sin el numero de componentes 3NN Y 2 NCOMP
#LR= 2
mapper= umap.UMAP(n_neighbors=5,min_dist=0.0,
                  n_components=2,learning_rate=2,
                  random_state=42,
                  metric="manhattan").fit_transform(fitvalT)

#training y test 

labels= hdbscan.HDBSCAN(
    min_samples=None,
    min_cluster_size=19, metric="manhattan").fit_predict(mapper)
#, metric="canberra"

clustered= (labels>=0)
plt.scatter(mapper[~clustered,0],
            mapper[~clustered,1],
            c=(0.5,0.5,0.5),
            s=10,
            alpha =0.5)

plt.scatter(mapper[clustered,0],
            mapper[clustered,1],
            c=labels[clustered],
            s=10,
            cmap="Spectral")

#grafica con todos los puntos juntos 
plt.scatter(mapper[:,0],mapper[:,1], c=labels, s=10, cmap="Spectral")

cbar=plt.colorbar(boundaries=np.arange(3)-0.5)
cbar.set_ticks(np.arange(7))
cbar.set_ticklabels

#solo falta poner el cuantitative assesment
#.shape es para ver la cantidad de datos osea 529
#se divide entre los clustered 

np.sum(clustered)/fitvalT.shape[0]




mapper2= umap.UMAP(n_neighbors=5,min_dist=0.0,
                  n_components=2,learning_rate=2,
                  random_state=42,
                  metric="manhattan").fit(fitvalT)

hover_data=pd.DataFrame({'index':fitvalT.index,
                         'label':labels,
                         'groups':groups.grouphist})

umap.plot.output_notebook()

#label

colors= Plasma[7]
p=umap.plot.interactive(mapper2,labels=hover_data.groups,
                        hover_data=hover_data, point_size=10,
                        color_key_cmap="plasma",background="black", )

color_mapper = (LinearColorMapper(palette=colors))


color_bar=ColorBar(color_mapper=color_mapper,
                   ticker=BasicTicker(desired_num_ticks=len(hover_data.groups)),
                   formatter=PrintfTickFormatter(format="%d%%"),
                   label_standoff=6,
                   major_label_text_font_size="7px",
                   location=(0,0))

p.add_layout(color_bar,'left')
output_file("umapgroups.html")
umap.plot.show(p)



#label


colors= Plasma[3]
p2=umap.plot.interactive(mapper2,labels=labels,
                        hover_data=hover_data, point_size=10,
                        color_key_cmap="plasma",background="black", )

color_mapper = (LinearColorMapper(palette=colors))


color_bar=ColorBar(color_mapper=color_mapper,
                   ticker=BasicTicker(desired_num_ticks=len(labels)),
                   formatter=PrintfTickFormatter(format="%d%%"),
                   label_standoff=6,
                   major_label_text_font_size="7px",
                   location=(0,0))

p.add_layout(color_bar,'left')
output_file("umaplabels.html")
umap.plot.show(p2)

save(p)
save(p2)


export=pd.DataFrame(labels)

export.to_csv('C:\\Users\\abe20\\OneDrive\\Documentos\\R\\TESIS\\labels5.csv')