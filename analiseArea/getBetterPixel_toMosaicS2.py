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
sys.setrecursionlimit(1000000000)

class ClassCalcIndicesSpectral(object):

    # default options
    bndInd = []  
    options = {        
        "bandas": ['B2','B3', 'B4', 'B8', 'B9', 'B11','B12', 'MSK_CLDPRB'],        
        "Feature" :[
                'evi','ndvi','ndwi','lai','gcvi',
                'osavi', 'isoil', 'gv',  'ndfia',
                'cvi','iia','npv','msi','soil','spri',
                'B2','B3', 'B4', 'B8', 'B11', 'B12'
        ],
        'bandasExtas' : [],
        "featAdd": [],     
        "FeatureG": [],
        "feat_Calcular": [], 
        'bandsFraction': ['gv','npv','soil','cloud','shade'] 
    }  

    CLOUD_FILTER_TOLER = 1
    CLD_PRB_THRESH = 30
    NIR_DRK_THRESH = 0.15
    B9_SHAD_THRESH = 230
    CLD_PRJ_DIST = 1
    dictClassifRef = {}
    CLOUD_FILTER = 70
    limiarProb = 0.45
    numRuido = None
    geomet = None
    imgUltima = None
    # dict_propCloudshadow = None
    maskVegetation = None
    maskRuido = None

    # cloudsMascara = True
    imgColClouds = None
    imgMosaic = None

    periodos = {
            # 'wet': ['-01-01', '-07-01'], 
            'dry': ['-06-01', '-01-01']
        }    

    def __init__(self, parametros, sOrbita, sTile, withgeomet):
        self.options = parametros
        self.sOrbita = sOrbita
        self.sTile = sTile
        imgColS2 = ee.ImageCollection(parametros['asset_S2harmonic']).filter(
                                    ee.Filter.eq('SENSING_ORBIT_NUMBER', int(sOrbita))).filter(
                                        ee.Filter.eq('MGRS_TILE', sTile)).filter(
                                            ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', self.CLOUD_FILTER))
        
        # get how many image will be processing and cleaning clouds 
        lstCloudsPercent = imgColS2.reduceColumns(ee.Reducer.toList(), ['CLOUDY_PIXEL_PERCENTAGE']).get('list').getInfo()
        numImgLoad = imgColS2.size().getInfo()
        print(f"we have {numImgLoad} images with band clouds ")
        numberimgMaskCloud = [kk for kk in lstCloudsPercent if kk > self.CLOUD_FILTER_TOLER]
        numberimgMaskCloud = len(numberimgMaskCloud) / 2

        self.footPrintimgCloudFull = imgColS2.sort('CLOUDY_PIXEL_PERCENTAGE', False).first().get('system:footprint')

        # get geometry  of reference 
        if withgeomet:
            self.geomet = ee.FeatureCollection(parametros['gradeS2Corr']).filter(
                                ee.Filter.eq('SENSING_ORBIT_NUMBER', int(kOrb))) .filter(
                                    ee.Filter.eq('MGRS_TILE', ntile)).geometry()
            
        else:
            self.geomet = imgColS2.sort('CLOUDY_PIXEL_PERCENTAGE', False).first().geometry()

        # Import and filter s2cloudless.
        s2_cloudless_col = (ee.ImageCollection(parametros['asset_S2Mask'])
                            .filterBounds(self.geomet))
        
        # Join asset of imagens Sentinel S2 and asset of s2cloudless
        imgColS2join = ee.ImageCollection(
                                ee.Join.saveFirst('s2cloudless').apply(**{
                                    'primary': imgColS2,
                                    'secondary': s2_cloudless_col,
                                    'condition': ee.Filter.equals(**{
                                        'leftField': 'system:index',
                                        'rightField': 'system:index'
                                    })
                                }))

        if parametros['showinfoImgCol']:            
            print(f"bands of the imageCollection {imgColS2join.first().bandNames().getInfo()} ")
            print("band CLouds ", ee.Image(imgColS2join.first().get('s2cloudless')).bandNames().getInfo())

        # add all bands of clouds 
        self.imColwithCLD = imgColS2join.map(lambda imS2 : self.add_cld_shdw_mask(imS2))

        self.imDarkPixel = self.imColwithCLD.select('dark_pixels').sum()
        self.imDarkPixel = self.imDarkPixel.gte(numberimgMaskCloud).rename('dark_pixels_est')


    def get_features (self, feat_lista):

        ls_features = []
        ls_feat_mosaic = []
        ls_all_feat = []
        for variavel in feat_lista:
            # print("variavel == ", variavel)
            # deixando a variavel sem sufixo
            var_limpa = variavel.replace("_d", "")          

            # selecionando os índices que interatuam na diferencias 
            # e precisam sere calculados no mosaico 

            if '_d' in variavel and variavel not in ls_features:
                ls_features.append(variavel.replace("_d", ""))            
            #  índices do mosaico que entram na classificação            
            # selecionando todos os índices que serão calculados na última imagem 
            if (var_limpa not in ls_all_feat and 
                var_limpa not in self.options["bandas"]):                
                ls_all_feat.append(var_limpa)
        
        self.options["Feature_Dif"] = copy.deepcopy(ls_features)        
        self.options["feat_Calcular"] = copy.deepcopy(ls_all_feat)

        for variavel in self.options["featAdd"]:
            var_limpa = variavel.replace("_m", "")
            ls_feat_mosaic.append(var_limpa)

        self.options['bandasExtas'] = ls_feat_mosaic
        print("serão selcionado os siguente indices para diferencia ")
        print(ls_features)
        print(" do mosaico se extrairão as siguentes features ") 
        print(self.options['bandasExtas'])       
        print(self.options["featAdd"])
        print(" as features a calcular serão ")
        print(ls_all_feat)
        
    
    def add_cloud_bands(self, img):

        pmt_focal = {
            'radius': 2,
            'kernelType': 'square',
            'iterations': 1
        }

        # Get s2cloudless image, subset the probability band.
        cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

        # Condition s2cloudless by the probability threshold value.
        is_cloud = cld_prb.gt(self.CLD_PRB_THRESH).rename('clouds')

        qa = img.select('QA60');
        # Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = 1 << 10;
        cirrusBitMask = 1 << 11;
        # Both flags should be set to zero, indicating clear conditions.
        maskCloudsQ60 = qa.bitwiseAnd(cloudBitMask).eq(1).And(
                    qa.bitwiseAnd(cirrusBitMask).eq(1)).rename('clouds_Q60');

        is_cloud = is_cloud.add(maskCloudsQ60).gt(0).rename('clouds')
        is_cloud = is_cloud.focalMax(** pmt_focal)
        
        # Add the cloud probability layer and cloud mask as image bands.
        return img.addBands(ee.Image(is_cloud))


    def add_shadow_bands(self, img):
        pmt_focal = {
            'radius': 1,
            'kernelType': 'square',
            'iterations': 1
        }
        # buiilding mask clouds 
        imWithCLD = self.add_cloud_bands(img)
        imWithCLD  =imWithCLD.select(['clouds'])
        imWithCLDInv = imWithCLD.eq(0).rename("cloudMaskInv")
        
        # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
        SR_BAND_SCALE = 1e4
        # .multiply(not_water)
        dark_pixels = img.select('B8').lt(self.NIR_DRK_THRESH * SR_BAND_SCALE).rename('dark_pixels')
        shadownMask = img.select("B9").lt(self.B9_SHAD_THRESH).rename('shadowsB9')
        shadownMask = shadownMask.focalMax(** pmt_focal)
        dark_pixels = dark_pixels.focalMax(** pmt_focal);
 
        ##########################################################################
        # Determine the direction to project cloud shadow from clouds           ##
        # we're working in a UTM projection.                                    ##
        ##########################################################################
        shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

        # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
        cld_proj = (imWithCLD.select('clouds').directionalDistanceTransform(shadow_azimuth, self.CLD_PRJ_DIST * 10)
                        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
                            .select('distance').mask().rename('cloud_transform'))

        # Identify the intersection of dark pixels with cloud shadow projection.
        shadows = cld_proj.multiply(dark_pixels).add(shadownMask).gt(0).rename('shadows')

        return img.addBands(ee.Image([imWithCLD, imWithCLDInv, dark_pixels, shadows]))


    # https://code.earthengine.google.com/9e24373e43b239e4310ecd1affc3190e
    def add_cld_shdw_mask(self, img):
        # Add cloud component bands.
        bandConst = img.select('B2').gt(1e6)
        percentCloud = ee.Number(img.get('CLOUDY_PIXEL_PERCENTAGE'))    
        # Add cloud shadow component bands.
        img_cloud_shadow = ee.Algorithms.If(
                                ee.Algorithms.IsEqual(percentCloud.gt(self.CLOUD_FILTER_TOLER), 1),
                                self.add_shadow_bands(img),
                                img.addBands(bandConst.rename('clouds')).addBands(
                                    bandConst.rename('shadows')).addBands(
                                        bandConst.rename('dark_pixels')).addBands(
                                            bandConst.eq(0).rename('cloudMaskInv'))
                            )
        img_cloud_shadow = ee.Image(img_cloud_shadow)
        
        return img_cloud_shadow #.updateMask(maskPixelClean)

    def apply_cloud_shadow_mask(self, img):
        bndCloudQA = ["cloudMaskInv"];
        maskClouds = img.select('clouds');
        maskshadows = img.select('shadows');
        maskdarkPixel = img.select('dark_pixels');
        maskshadows = maskClouds.add(maskshadows).add(maskdarkPixel).subtract(self.imDarkPixel).gt(0);
        maskPixelClean = maskshadows.eq(0);  #.add(is_cld_shdw)
        maskPixelClean = maskPixelClean.focalMax(2);

        return img.select(self.options['bandVis']).updateMask(maskPixelClean).addBands(img.select(bndCloudQA));
    

    def applicateCoeficientes (self, image):
        intercept = [-0.004, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012];
        slope = [0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003];

        bndCloudQA = ["cloudMaskInv"];
        imgsCloudQA = image.select(bndCloudQA);
  
        imagecoor = image.select(self.options['bandVis']).multiply(slope);
        imagecoor = imagecoor.add(intercept);
    
        return imgsCloudQA.addBands(imagecoor);

    def building_mosaic_andExport(self, nyear):
        # 'bandVis': ["B2","B3","B4","B8","B11","B12"],
        bndName = ['blue','green','red','nir','swir1','swir2']

        imgBase = ee.Image().toUint16()
        bandasBase = []
        for station, period in self.periodos.items():
            print(" --------- Processing station = {} and period = {} ------------- ".format(station, period))
            date_start = str(nyear) + period[0]
            date_end = str(nyear) + period[1]
            if station == 'dry':
                date_end = str(nyear + 1) + period[1]
            

            imgColPeriod = self.imColwithCLD.filterDate(date_start, date_end)             
            
            imgColPeriod = imgColPeriod.map(self.apply_cloud_shadow_mask)
            imgColPeriodCorr = imgColPeriod.map(self.applicateCoeficientes)

            imgPeriodQM = imgColPeriodCorr.qualityMosaic("cloudMaskInv")
            # print("know the bands imgPeriodQM ", imgPeriodQM.bandNames().getInfo())
            bndQual = [kk + "_QM_" + station for kk in bndName]
            bandasBase += bndQual
            imgPeriodQM = imgPeriodQM.select(self.options['bandVis'], bndQual).toUint16()

            imgPeriodMedian = imgColPeriodCorr.median()
            # print("know the bands imgPeriodMedian ", imgPeriodMedian.bandNames().getInfo())
            # bndMedian = [kk + "_median" for kk in bndName]
            bndMedianPer = [kk + "_median_" + station for kk in bndName]
            bandasBase += bndMedianPer
            imgPeriodMedian = imgPeriodMedian.select(self.options['bandVis'], bndMedianPer).toUint16()

            # imgPeriodMin = imgColPeriodCorr.reduce(ee.Reducer.min(), 4)
            # # print("know the bands imgPeriodMin ", imgPeriodMin.bandNames().getInfo())
            # bndMin = [kk + "_min"  for kk in self.options['bandVis']]
            # bndMinPer = [kk + "_min_" + station for kk in bndName]
            # bandasBase += bndMinPer
            # imgPeriodMin = imgPeriodMin.select(bndMin, bndMinPer).toUint16()

            # imgPeriodMax = imgColPeriodCorr.reduce(ee.Reducer.max(), 4)
            # # print("know the bands imgPeriodMax ", imgPeriodMax.bandNames().getInfo())
            # bndMax = [kk + "_max" for kk in self.options['bandVis']]
            # bndMaxPer = [kk + "_max_" + station for kk in bndName]
            # bandasBase += bndMaxPer
            # imgPeriodMax = imgPeriodMax.select(bndMax, bndMaxPer).toUint16()

            # imgPeriodstdDev = imgColPeriodCorr.reduce(ee.Reducer.stdDev(), 4)
            # # print("know the bands imgPeriodstdDev ", imgPeriodstdDev.bandNames().getInfo())
            # bndstdDev = [kk + "_stdDev" for kk in self.options['bandVis']] 
            # bndstdDevPer = [kk + "_stdDev_" + station for kk in bndName] 
            # bandasBase += bndstdDevPer 
            # imgPeriodstdDev = imgPeriodstdDev.select(bndstdDev, bndstdDevPer).toUint16()

            imgBase = imgBase.addBands(imgPeriodQM).addBands(imgPeriodMedian)
                            # .addBands(imgPeriodMin
                            # ).addBands(imgPeriodMax).addBands(imgPeriodstdDev)
                        

            # sys.exit()

        imgBase = imgBase.set(
                                'system:footprint', self.footPrintimgCloudFull,
                                'SENSING_ORBIT_NUMBER', self.sOrbita,
                                'MGRS_TILE', self.sTile,
                                'year', nyear,
                                'version', '4',
                                'biome', 'CAATINGA',
                                'satellite', 'S2_HARMONIZED',
                                'source', 'geodatin'
                )
        imgBase =imgBase.select(bandasBase)
        nameImg = 'CAATINGA-' + str(self.sOrbita) + '-' + str(self.sTile) + '-' + str(nyear) + '-' + 'S2-4'
        self.exportarImagem(imgBase, nameImg)

    def exportarImagem(self, imgU, nameAl):

        IdAsset = self.options['output'] + "/" + nameAl           
        optExp = {
            'image': imgU,
            'description': nameAl, 
            'assetId':IdAsset, 
            'pyramidingPolicy': {".default": "mode"},  
            'region': self.geomet.getInfo()['coordinates'],
            'scale': 10,
            'maxPixels': 1e13 
        }

        task = ee.batch.Export.image.toAsset(**optExp)    
        task.start()
        
        print ("salvando ... !" + nameAl)

    
"""
* Function to mask clouds using the Sentinel-2 QA band
* @param {ee.Image} image Sentinel-2 image
* @return {ee.Image} cloud masked Sentinel-2 image
"""

#========================METODOS=============================
def gerenciador(cont, paramet):
    #0, 18, 36, 54]
    #=====================================
    # gerenciador de contas para controlar 
    # processos task no gee
    # cada conta vai rodar com 18 cartas X 3 anos
    #=====================================
    numberofChange = [kk for kk in paramet['conta'].keys()]
    
    if str(cont) in numberofChange:

        print("conta ativa >> {} <<".format(paramet['conta'][str(cont)]))        
        gee.switch_user(paramet['conta'][str(cont)])
        gee.init()        
        gee.tasks(n= paramet['numeroTask'], return_list= True)        
    
    elif cont > paramet['numeroLimit']:
        cont = -1
    
    cont += 1    
            
    return cont


###########################################################################################
##    Computes spectral indices of cloudyness and take the minimum of them.              ##
##                                                                                       ##                               
##    Each spectral index is fairly lenient because the group minimum                    ##
##    is a somewhat stringent comparison policy. side note -> this                       ##
##    seems like a job for machine learning :)                                           ##
##    originally written by Matt Hancher for Landsat imagery                             ##
##    adapted to Sentinel by Chris Hewig and Ian Housman                                 ##
###########################################################################################


param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'output': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    # 'output': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/mosaics-CAATINGA-4',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  # 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY",
    'bandVis': ["B2","B3","B4","B8","B11","B12"],
    'showinfoImgCol': False,
    'numeroTask': 6,
    'numeroLimit': 100,
    'conta' : {
        '0': 'caatinga01',
        '7': 'caatinga02',
        '14': 'caatinga03',
        '21': 'caatinga04',
        '28': 'caatinga05',        
        '35': 'solkan1201',
        # '36': 'rodrigo',
        '42': 'diegoGmail',    
        '50': 'superconta'
    },
}
# 266 Tiles Orbitas
geomOrbTile = ee.FeatureCollection(param['gradeS2Corr'])
contAuth = 50
contAuth = gerenciador(contAuth, param)
vjson = open('dict_Tile_Orbita_2023_Sentinel2.json')
dictTileOrb = json.load(vjson)
cc = 0           
for kOrb, myLst in dictTileOrb.items():    
    for ntile in myLst:               
        cc += 1
        if cc > 130:
            print("================================================================")
            mgeomet = geomOrbTile.filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', int(kOrb))) .filter(
                                        ee.Filter.eq('MGRS_TILE', ntile))        
            # print(f"{cc} ---PROCESSING ->  Orbita  {kOrb}   Tile  {ntile} --- { mgeomet.size().getInfo()}")
            print(f"{cc} ---PROCESSING ->  Orbita  {kOrb}   Tile  {ntile} --- ")

            ClassCalcIndicesSpectralMosiac =  ClassCalcIndicesSpectral(param, kOrb, ntile, True)
            for yy in range(2015, 2016):
                ClassCalcIndicesSpectralMosiac.building_mosaic_andExport(yy)
                contAuth = gerenciador(contAuth, param)

        