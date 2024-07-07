#-*- coding utf-8 -*-
import ee
import gee 
import os
import sys
import random
from datetime import date
import copy
import math
import json
import collections
collections.Callable = collections.abc.Callable
try:
  ee.Initialize()
  print('The Earth Engine package initialized successfully!')
except ee.EEException as e:
  print('The Earth Engine package failed to initialize!')
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise


param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  # 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY",
    'bandsS2' : ["B2","B3","B4","B8","B9","B10","QA10"],
    'bandVis': ["B2","B3","B4","B8","B9","B10"],

}
geomOrbTile = ee.FeatureCollection(param['gradeS2Corr'])
lstOrbTiles = geomOrbTile.reduceColumns(ee.Reducer.toList(2), ['SENSING_ORBIT_NUMBER','MGRS_TILE']).get('list').getInfo()

mdictOrb = {}

for cc, obrTile in enumerate(lstOrbTiles):
    print(obrTile)
    lstOrb = [kk for kk in mdictOrb.keys()]

    if str(obrTile[0]) not in lstOrb:
        mdictOrb[str(obrTile[0])] = [obrTile[1]]
    else:
       lsttmp = mdictOrb[str(obrTile[0])]
       lsttmp.append(obrTile[1])
       mdictOrb[str(obrTile[0])] = lsttmp

for kk, vv in mdictOrb.items():
    print(kk, vv)

with open("dict_Tile_Orbita_2023_Sentinel2.json", "w") as file:
    json.dump(mdictOrb, file)
print("the dictionary was saved ")