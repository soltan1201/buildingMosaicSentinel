#-*- coding utf-8 -*-
import ee
import os
import sys
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


poligons =  [
        [[-40.231804078876614, -13.173196663547518],
          [-40.231804078876614, -14.772325765997135],
          [-38.248771852314114, -14.772325765997135],
          [-38.248771852314114, -13.173196663547518]]
    ];
geometry =  ee.Geometry.Polygon(poligons, None, False)


def exportImg (clipIm, nnewNome, nbuffer):
    param ={
        'image': clipIm.toInt8(), #//.toUint8(),
        'description': nnewNome,
        'folder': "MAPBIOMAS-EXPORT",
        'maxPixels':1e13,
        'scale': 30,
        'region': nbuffer
    }
    task = ee.batch.Export.image.toDrive(**param)    
    task.start()
        
    print (f"salvando {nnewNome}... !")


anos = [kk for kk in range(1985,2022)]
print("lista de anos ", anos);
assetCol8 = "projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1";
collectionx = ee.Image(assetCol8);
print(collectionx, 'collection')
classMapbiomas = [1,3,4,5,6,49,10,11,12,32,29,50,13,14,15,18,19,39,20,40,62,41,36,46,47,35,48, 9, 21, 22, 23, 24, 30, 25, 26, 33, 31, 27];
reclass = [1,1,1,1,1, 2, 1, 1, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3,  4,  4,  5,  4,  4,  6,  6,  6, 27];
imgRemap = ee.Image().byte();

for ano in anos:
    bandAct = 'classification_'+ str(ano)
    mosaico1 = collectionx.select(bandAct).clip(geometry);
    mosaico= ee.Image(mosaico1).remap(classMapbiomas, reclass);
    imgRemap = imgRemap.addBands(mosaico.rename(bandAct));    


collection = collectionx;
coordlist =[
        ['APA_1',-38.999000 , -13.708000],['APA_2',-39.019000 , -13.705000],
        ['APA_3',-39.023000 , -13.721000],['APA_4',-39.027000 , -13.736000],
        ['APA_5',-39.028000 , -13.708000],['APA_6',-39.042000 , -13.695000],
        ['APA_7',-39.055000 , -13.690000],['APA_8',-39.073000 , -13.678000],
        ['APA_9',-39.239014 , -13.845203],['APA_10',-39.449997 , -13.892000],
        ['REM_1',-39.218008 , -13.813992],['REM_2',-39.239014 , -13.845203],
        ['APA_11',-39.460000 , -13.906000],['APA_12',-39.479914 , -13.947525],
        ['REM_3',-39.168961 , -13.782947],['REM_4',-39.164000 , -13.777000],
        ['APA_13',-39.477344 , -13.937703],['APA_14',-39.481250 , -13.943169],
        ['APA_15',-39.310000 , -13.785000],['APA_16',-39.309000 , -13.789000],
        ['APA_17',-39.471000 , -13.936000],['APA_18',-39.469000 , -13.931000],
        ['REM_5',-39.178253 , -13.793558],['REM_6',-39.174878 , -13.785056],
        ['REM_7',-39.233000 , -13.854000],['REM_8',-39.234000 , -13.861000],
        ['REM_9',-39.239908 , -13.854033],['REM_10',-39.236197 , -13.868000],
        ['REM_11',-39.240056 , -13.859661],['REM_12',-39.248781 , -13.862581],
        ['REM_13',-39.208036 , -13.826569],['REM_14',-39.164000 , -13.777000],
        ['REM_15',-39.230000 , -13.850000],['REM_16',-39.233756 , -13.832536],
        ['REM_17',-39.221000 , -13.826000],['REM_18',-39.174878 , -13.785056]
    ];

print("points list ", coordlist);
tambuffer = ee.List.sequence(500, 2500, 500).getInfo();
for ano in anos:
    nomeBand = "classification_" + str(ano);
    imgTemp = imgRemap.select(nomeBand);
    for buftam in tambuffer:          
        for cc, mPoint in enumerate(coordlist):            
            print(f"{cc} Processing coordenadas {mPoint} in buffer {buftam}")
            buffer = ee.Geometry.Point([mPoint[1], mPoint[2]]).buffer(buftam , 10);
            clipa = imgTemp.clip(buffer)
            clipa = clipa.set('banda', nomeBand,'ano', ano, 'buffer_size', buftam);
            newNome = "Itu_Igrap_Buffer_clipPoint_" + mPoint[0] + "_" + nomeBand + "_" + str(buftam);

            exportImg(clipa, newNome, buffer);
                     

