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
        "FeatureG" :[ ],
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
            'wet': ['-01-01', '-07-01'], 
            'dry': ['-07-01', '-01-01']
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
        numberimgMaskCloud = [kk for kk in lstCloudsPercent if kk > self.CLOUD_FILTER_TOLER]
        numberimgMaskCloud = len(numberimgMaskCloud) / 2

        self.footPrintimgCloudFull = imgColS2.sort('CLOUDY_PIXEL_PERCENTAGE', False).first().get('system:footprint')

        # get geometry  of reference 
        if withgeomet:
            self.geomet = ee.FeatureCollection(parametros['gradeS2Corr']).filter(
                                ee.Filter.eq('SENSING_ORBIT_NUMBER', int(kOrb))) .filter(
                                    ee.Filter.eq('MGRS_TILE', ntile)).geometry()
            
        else:
            self.geomet = imgColS2.first().geometry()

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
            print(f"we have {numImgLoad} images with band clouds ")
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
        
        percentCloud = ee.Number(img.get('CLOUDY_PIXEL_PERCENTAGE'))    
        # Add cloud shadow component bands.
        img_cloud_shadow = ee.Algorithms.If(
                                ee.Algorithms.IsEqual(percentCloud.gt(self.CLOUD_FILTER_TOLER), 1),
                                self.add_shadow_bands(img),
                                img.addBands(ee.Image.constant(0).rename('clouds')).addBands(
                                    ee.Image.constant(0).rename('shadows')).addBands(
                                        ee.Image.constant(0).rename('dark_pixels')).addBands(
                                        ee.Image.constant(1).rename('cloudMaskInv'))
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
        # self.sOrbita = sOrbita
        # self.sTile = sTile#
        # periodos = {
        #     'wet': ['-01-01', '-07-01'], 
        #     'dry': ['-07-01', '-01-01']
        # }
        imgBase = ee.Image().toUint16()

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
            bndQual = [kk + "_QM_" + station for kk in bndName]
            imgPeriodQM = imgPeriodQM.select(self.options['bandVis'], bndQual) 

            imgPeriodMedian = imgColPeriodCorr.median()
            bndMedian = [kk + "_median_" + station for kk in bndName]
            imgPeriodMedian = imgPeriodMedian.select(self.options['bandVis'], bndMedian)

            imgPeriodMin = imgColPeriodCorr.reduce(ee.Reducer.min(), 4)
            bndMin = [kk + "_min_" + station for kk in bndName]
            imgPeriodMin = imgPeriodMin.select(self.options['bandVis'], bndMin)

            imgPeriodMax = imgColPeriodCorr.reduce(ee.Reducer.max(), 4)
            bndMax = [kk + "_max_" + station for kk in bndName]
            imgPeriodMax = imgPeriodMax.select(self.options['bandVis'], bndMax)

            imgPeriodstdDev = imgColPeriodCorr.reduce(ee.Reducer.stdDev(), 4)
            bndstdDev = [kk + "_stdDev_" + station for kk in bndName]  
            imgPeriodstdDev = imgPeriodstdDev.select(self.options['bandVis'], bndstdDev)

            imgBase = imgBase.addBands(imgPeriodQM).addBands(imgPeriodMedian
                        ).addBands(imgPeriodMin).addBands(imgPeriodMax).addBands(imgPeriodstdDev)

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
        nameImg = 'CAATINGA-' + str(self.sOrbita) + '-' + str(self.sTile) + '-' + str(nyear) + '-' + 'S2-4'


    def agregateBandsIndexNDVI(self, img):
    
        ndviImg = img.expression("float(b('B8') - b('B12')) / (b('B8') + b('B12'))")\
            .rename(['ndvi'])       

        return img.addBands(ndviImg)

    def agregateBandsIndexWater(self, img):
    
        ndwiImg = img.expression("float(b('B8') - b('B12')) / (b('B8') + b('B12'))")\
            .rename(['ndwi'])       

        return img.addBands(ndwiImg)

    def agregateBandsIndexEVI(self, img):
            
        eviImg = img.expression(
            "float(2.4 * (b('B8') - b('B4')) / (1 + b('B8') + b('B4')))")\
                .rename(['evi'])     
        
        return img.addBands(eviImg)
    
    def IndiceIndicadorAgua(self, img):
    
        iiaImg = img.expression(
                            "float((b('B3') - 4 *  b('B8')) / (b('B3') + 4 *  b('B8')))"
                        ).rename("iia")
        
        return img.addBands(iiaImg)
    
    def agregateBandsIndexLAI(self, img):
    
        laiImg = img.expression(
            "(3.618 * float(b('evi') - 0.118))")\
                .rename(['lai'])     
    
        return img.addBands(laiImg)
    
    def agregateBandsIndexGCVI(self, img):
    
        gcviImgA = img.expression(
            "float(b('B8')) / (b('B3')) - 1")\
                .rename(['gcvi'])        
        
        return img.addBands(gcviImgA)
    
    def agregateBandsIndexOSAVI(self,img):
    
        osaviImg = img.expression(
            "float(b('B8') - b('B4')) / (0.16 + b('B8') + b('B4'))")\
                .rename(['osavi'])        
        
        return img.addBands(osaviImg)
    
    def agregateBandsIndexSoil(self, img):
        
        soilImg = img.expression(
            "float(b('B8') - b('B3')) / (b('B8') + b('B3'))")\
                .rename(['isoil'])       
        
        return img.addBands(soilImg)    
    
    # Moisture Stress Index (MSI)
    def agregateBandsIndexMSI(self, img):
    
        msiImg = img.expression(
            "float( b('B11') / b('B8'))").rename(['msi']) 
        
        return img.addBands(msiImg)
    
    # char soil Index
    def agregateBandsIndexCSI(self, img):
    
        csiImg = img.expression(
            "float( b('B8') / b('B11'))").rename(['csi']) 
        
        return img.addBands(csiImg)

    def agregateBandsIndexBAI(self, img):
    
        baiImg = img.expression(
            "float(1) / ((0.1 - b('B4'))**2 + (0.06 - b('B8'))**2)")\
                .rename(['bai']) 
        
        return img.addBands(baiImg)
    
    def agregateBandsIndexsPRI(self, img):
        
        priImg = img.expression(
                                "float((b('B3') - b('B2')) / (b('B3') + b('B2')))"
                            ).rename(['pri'])   
        spriImg =  priImg.expression(
                                "float((b('pri') + 1) / 2)").rename(['spri'])  
    
        return img.addBands(spriImg)
    
    # Chlorophyll vegetation index
    def agregateBandsIndexCVI(self, img):
    
        cviImgA = img.expression(
            "float(b('B8') * (b('B3') / (b('B2') * b('B2'))))").rename(['cvi'])        
        
        return img.addBands(cviImgA)

    def AutomatedWaterExtractionIndex(self, img):
    
        awei = img.expression(
                            "float(4 * (b('B3') - b('B12')) - (0.25 * b('B8') + 2.75 * b('B11')))"
                        ).rename("awei")          
        
        return img.addBands(awei)
    
     # Tasselled Cap - brightness 
    
    def agregateBandsIndexBrightness(self, img):
    
        tasselledCapImg = img.expression(
            "float(0.3037 * b('B2') + 0.2793 * b('B3') + 0.4743 * b('B4')  + 0.5585 * b('B8') + 0.5082 * b('B11') +  0.1863 * b('B12'))")\
                .rename(['brightness']) 
        
        return img.addBands(tasselledCapImg)
    
    # Tasselled Cap - wetness 
    def agregateBandsIndexwetness(self, img):
    
        tasselledCapImg = img.expression(
            "float(0.1509 * b('B2') + 0.1973 * b('B3') + 0.3279 * b('B4')  + 0.3406 * b('B8') + 0.7112 * b('B11') +  0.4572 * b('B12'))")\
                .rename(['wetness']) 
        
        return img.addBands(tasselledCapImg)
    
    def agregateBandsIndexGVMI(self, img):
        
        gvmiImg = img.expression(
                        "float ((b('B8')  + 0.1) - (b('B11') + 0.02)) / ((b('B8') + 0.1) + (b('B11') + 0.02))" 
                    ).rename(['gvmi'])     
    
        return img.addBands(gvmiImg)
    
    def agregateBandsIndexRATIO(self, img):
    
        ratioImg = img.expression("float(b('B8') / b('B4'))")\
                            .rename(['ratio'])      

        return img.addBands(ratioImg)
    
    def agregateBandsIndexRVI(self, img):
    
        rviImg = img.expression("float(b('B4') / b('B8'))")\
                        .rename(['rvi'])       

        return img.addBands(rviImg)

    def agregateBandsTexturasGLCM(self, img):
        
        img = img.toInt()
                
        textura2 = img.select('B8').glcmTexture(3)            
        contrast = textura2.select('B8_contrast').divide(1000).rename('contrast')    

        return contrast
        
    def agregateBandsgetFractions(self, img):

       # Define endmembers
        endmembers =  [
            [ 119.0,  475.0,  169.0, 6250.0, 2399.0,  675.0], #/*gv*/
            [1514.0, 1597.0, 1421.0, 3053.0, 7707.0, 1975.0], #/*npv*/
            [1799.0, 2479.0, 3158.0, 5437.0, 7707.0, 6646.0], #/*soil*/
            [4031.0, 8714.0, 7900.0, 8989.0, 7002.0, 6607.0], #/*cloud*/
            [   0.0,    0.0,    0.0,    0.0,    0.0,    0.0]  #/*Shade*/
        ]

        # Uminxing data
        bandas6 = ['B2','B3', 'B4', 'B8', 'B11','B12']
        fractions = ee.Image(img).select(bandas6)\
                        .unmix(endmembers= endmembers).float()
        
        fractions = fractions.select([0,1,2,3,4], self.options['bandsFraction'])

        ndfia = fractions.expression(
            "float(b('gv') - b('soil')) / float( b('gv') + b('npv') + b('soil'))")
        
        ndfia = ndfia.rename('ndfia')
        
        return fractions.select(['gv', 'npv','soil']).addBands(ndfia)
    

    def CalculateIndice(self, imageW, propCloud):         
        # print("presentes bandas da imagem == ", imageW.bandNames().getInfo()) 
        if propCloud:
            imageW = self.maskS2clouds(imageW)
        print("olhando as bandas da imagem == ", ee.Image(imageW).bandNames().getInfo())
        # self.cloudsMascara == False     
        imageW = ee.Image(imageW).set('system:footprint', self.geomet)
        # por causa do bucket  a imagem sai com [0, 10.000]
        
        if 'contrast' in self.options["FeatureG"]:
            imageT = self.agregateBandsTexturasGLCM(imageW)  
             
                
        # imagem em Int16 com valores inteiros ate 10000     
        if ('gv' in self.options["feat_Calcular"] or 
                'npv' in self.options["feat_Calcular"] or 
                    'soil' in self.options["feat_Calcular"] or  
                        'ndfia' in self.options["feat_Calcular"]):
            
            imageF = self.agregateBandsgetFractions(imageW)        
        
            print("calculadas as 🔰 FRAÇÕES 🔰 ")
        
        imageW = imageW.divide(10000).toFloat()   
        if 'evi' in self.options["feat_Calcular"] or 'lai' in self.options["feat_Calcular"]:
            imageW = self.agregateBandsIndexEVI(imageW)            
        if 'awei' in self.options["feat_Calcular"] :
            imageW = self.AutomatedWaterExtractionIndex(imageW)
        if 'brightness' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexBrightness(imageW)        
        if 'cvi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexCVI(imageW)       
        if 'gcvi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexGCVI(imageW)   
        if 'gvmi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexGVMI(imageW)
        if 'iia' in self.options["feat_Calcular"] :
            imageW = self.IndiceIndicadorAgua(imageW)
        if 'isoil' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexSoil(imageW)
        if 'lai' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexLAI(imageW)
        if 'msi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexMSI(imageW)
        if 'ndvi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexNDVI(imageW)
        if 'ndwi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexWater(imageW)
        if 'osavi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexOSAVI(imageW)
        if 'ratio' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexRATIO(imageW)
        if 'rvi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexRVI(imageW)
        if 'csi' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexCSI(imageW)
        if 'spri' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexsPRI(imageW)
        if 'wetness' in self.options["feat_Calcular"] :
            imageW = self.agregateBandsIndexwetness(imageW)         
        
        imageW = imageW.set('system:footprint', self.geomet)
        imageW = imageW.addBands(imageF)
        # print("Indices calculados na Imagem: \n ====>", imageW.bandNames().getInfo())

        if 'contrast' in self.options["FeatureG"]:
            return imageW.addBands(imageT)

        return imageW
    
param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  # 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY",
    'bandVis': ["B2","B3","B4","B8",,"B11","B12"],
    'showinfoImgCol': True
}
geomOrbTile = ee.FeatureCollection(param['gradeS2Corr'])

vjson = open('dict_Tile_Orbita_2023_Sentinel2.json')
dictTileOrb = json.load(vjson)
cc = 0           
for kOrb, myLst in dictTileOrb.items():    
    for ntile in myLst:               
        cc += 1
        print("================================================================")
        mgeomet = geomOrbTile.filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', int(kOrb))) .filter(
                                    ee.Filter.eq('MGRS_TILE', ntile))        
        print(f"{cc} ---PROCESSING ->  Orbita  {kOrb}   Tile  {ntile} --- { mgeomet.size().getInfo()}")

        ClassCalcIndicesSpectralMosiac =  ClassCalcIndicesSpectral(param, kOrb, ntile, True)
        for yy in range(2016, 2017):
            ClassCalcIndicesSpectralMosiac.building_mosaic_andExport(yy)

        sys.exit()