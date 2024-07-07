
// https://code.earthengine.google.com/09e10a3e046c6b338ee7dd2c80157616
var vis = {
    mosaico: {
        min:30, 
        max: 3500,
        bands: ["red_median","green_median","blue_median"]
    },
    imgS2: {
        min:20, 
        max: 2500,
        bands: ["B4","B3","B2"]
    }, 
    maskCC: {
        min: 0, 
        max: 1,
        palette: "000000,00FF00,FFFA00"
    }
}
// 'bandVis': ["B2","B3","B4","B8","B9","B10"],
var intercept = [-0.004, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012];
var slope = [0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003];
slope = ee.Image.constant(slope);
intercept = ee.Image.constant(intercept);

var addNDVI = function(image) {
    var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
    return image.addBands(ndvi);
}; 
var add_shadow_bands= function(img){
    var pmt_focal = {
        'radius': 1,
        'kernelType': 'square',
        'iterations': 1
    }

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4
    var NIR_DRK_THRESH = 0.15
    var CLD_PRJ_DIST = 1
    print("valor SR_BAND_SCALE", SR_BAND_SCALE)
    // .multiply(not_water)
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH * SR_BAND_SCALE).rename('dark_pixels');
 
    var shadownMask = img.select("B9").lt(230).rename('shadowsB9')
    shadownMask = shadownMask.focalMax(pmt_focal);
    // pmt_focal['radius'] = 6
    dark_pixels = dark_pixels.focalMax(pmt_focal);

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('cloudMask').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
                    .reproject({'crs': img.select(0).projection(), 'scale': 100})
                        .select('distance').mask().rename('cloud_transform'))

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).add(shadownMask).gt(0).rename('shadows')
    // shadows = shadows.add(shadownMask).gt(0)

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([shadows, dark_pixels]))
}

var maskingCoverCloud = function(image){
    var pmt_focal = {
        'radius': 1,
        'kernelType': 'square',
        'iterations': 1
    }
    var imgCloud = ee.Image(image.select('probability')).gt(30).rename('cloudMask');
    var imgCloudInv = imgCloud.eq(0).rename('cloudMaskInv');
    var imgCloudShadow = add_shadow_bands(image.addBands(imgCloud));
    // var maskPixelClean = imgCloudShadow.select('shadows').add(imgCloudShadow.select('shadowsB9')
    //                             ).add(imgCloud).focalMax(pmt_focal).eq(0) 
    return image.addBands(imgCloudInv).addBands(imgCloudShadow.select(['cloudMask','shadows', 'dark_pixels']));  //.updateMask(maskPixelClean)
};
var apply_mask_coverCloudsShadows = function (imcolCS , imgDarkPix){
      var pmt_focal = {
          'radius': 1,
          'kernelType': 'square',
          'iterations': 1
      }
      var bndCloudQA = ["QA60", "cloudMask",'cloudMaskInv',"shadows", "dark_pixels"];
      var imgColClean = imcolCS.map(function(image){
            var maskPixelClean = image.select('shadows').add(image.select('dark_pixels')
                                ).add(image.select('cloudMask')).focalMax(pmt_focal).gt(0)
            maskPixelClean = maskPixelClean.subtract(imgDarkPix).eq(0);
            maskPixelClean = maskPixelClean.focalMax(2);
            return image.select(param.bandVis).updateMask(maskPixelClean).addBands(image.select(bndCloudQA));
      })
      
      return imgColClean;
}

// Join two collections on their 'system:index' property.
// The propertyName parameter is the name of the property
// that references the joined image.
function indexJoin(collectionA, collectionB, propertyName) {
    var joined = ee.ImageCollection(ee.Join.saveFirst(propertyName).apply({
        primary: collectionA,
        secondary: collectionB,
        condition: ee.Filter.equals({
            leftField: 'system:index',
            rightField: 'system:index'})
    }));
    // Merge the bands of the joined image.
    return joined.map(function(image) {
        return image.addBands(ee.Image(image.get(propertyName)));
    });
}

var applicateCeof = function(image){
    var bndCloudQA = ["QA60", "cloudMask","cloudMaskInv","shadows", "dark_pixels"];
    var imgsCloudQA = image.select(bndCloudQA);
  
    var imagecoor = image.select(param.bandVis).multiply(slope);
    imagecoor = imagecoor.add(intercept);
    
    return imgsCloudQA.addBands(imagecoor);
};


var param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY",
    'bandVis': ["B2","B3","B4","B8","B11","B12"],
}
var sOrbita = 138;
var sTile = '23KMB';
var year = 2016;
var periodos = {
    'wet': ['-01-01', '-07-01'], 
    'dry': ['-07-01', '-01-01']
};
var date_start = year.toString() + periodos['wet'][0];
var date_end =  year.toString() + periodos['wet'][1];


var geomet = ee.FeatureCollection(param['gradeS2Corr']).filter(
                ee.Filter.eq('SENSING_ORBIT_NUMBER', sOrbita)) .filter(
                    ee.Filter.eq('MGRS_TILE', sTile)).geometry();

var s2c = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
                    .filterBounds(geomet);
print("show the first images ", s2c.first())                    
var imColS2 = ee.ImageCollection(param.asset_S2harmonic).filter(
                ee.Filter.eq('SENSING_ORBIT_NUMBER', sOrbita)).filter(
                    ee.Filter.eq('MGRS_TILE', sTile)).filter(
                        ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 70)).filterDate(
                            date_start, date_end);
                            
var numberImg = imColS2.size().getInfo();
print("know numbers of images ", imColS2.sort('CLOUDY_PIXEL_PERCENTAGE', false));

// Join the cloud probability dataset to surface reflectance.
var withCloudProbability = indexJoin(imColS2, s2c, 'cloud_probability');
withCloudProbability = withCloudProbability.map(maskingCoverCloud);
print("know numbers of bands in images withCloudProbability ", withCloudProbability);

var imgdarkPixelest = withCloudProbability.select('dark_pixels').sum()
print(numberImg - 2);
// 
imgdarkPixelest = imgdarkPixelest.gt(numberImg/2 ); // 
withCloudProbability = apply_mask_coverCloudsShadows(withCloudProbability, imgdarkPixelest)
print("know numbers of bands in images withCloudClean", withCloudProbability);

var imgCleanedCorr = withCloudProbability.map(applicateCeof);
print("know numbers of bands in images of imgCleanedCorr ", imgCleanedCorr);
var imgCleanedCorrInd = imgCleanedCorr.map(addNDVI);
// Make a "greenest" pixel composite.
var greenest = imgCleanedCorrInd.qualityMosaic('NDVI');
var cloudsnest = imgCleanedCorr.qualityMosaic("cloudMask");
var cloudsInvsnest = imgCleanedCorr.qualityMosaic("cloudMaskInv");
print("image quality mosaic ", cloudsInvsnest)

// Display the result.
// Map.centerObject(geomet, 9);
// var ndviParams = {min: -1, max: 1, palette: ['blue', 'white', 'green']};

Map.addLayer(withCloudProbability.mosaic(), vis.imgS2, "Mosaic free Cloud");
Map.addLayer(imgCleanedCorr.mosaic(), vis.imgS2, "Mosaic Corregido");
Map.addLayer(imgCleanedCorr.median(), vis.imgS2, "Mosaic Median");
Map.addLayer(greenest, vis.imgS2, "Mosaic qualityMosaic NDVI", false);
Map.addLayer(cloudsnest, vis.imgS2, "Mosaic qualçityMosaic ", false);
Map.addLayer(cloudsInvsnest, vis.imgS2, "Mosaic qualçityMosaic CloudsInv");
Map.addLayer(imgdarkPixelest.selfMask(), {min:0, max: 1, palette: 'red'}, 'darkPixel estavel', false)

