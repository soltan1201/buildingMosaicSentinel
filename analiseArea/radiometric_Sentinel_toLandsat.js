// //------------------------------------     blue          green        red          NIR        SWIR-1    SWIR-2 
// 	CString Landsat8_Bands[] = { "band02" ,"band03","band04","band05" ,"band06","band07" };
// 	CString Sentinel2_Bands[] = { "B02" ,       "B03",       "B04",      "B8A" ,       "B11",      "B12" };
// 	CString Landsat8_outputBands[] = { "blue" ,"green","red","NIR" ,"SWIR1","SWIR2" };
// 	double paraBandAdjust_slope_V14A[6] = {  };
// 	double paraBandAdjust_intercept_V14A[6] = {  };
// 	int iTotalCount = 365 * 4 + 121; // from 20150101 to 20190430
// 	int32 start[2] = { 0, 0 };
// 	int32 end[2] = { 3660, 3660 };  // Total size

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

var param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY",
    'bandsS2' : ["B2","B3","B4","B8","B9","B10","QA10"],
    'bandVis': ["B2","B3","B4","B8","B9","B10"],
}

var intercept = [-0.004, -0.0009, 0.0009, -0.0001, -0.0011, -0.0012];
var slope = [0.9778, 1.0053, 0.9765, 0.9983, 0.9987, 1.003];
slope = ee.Image.constant(slope);
intercept = ee.Image.constant(intercept);

// cloudHeights = ee.List.sequence(200, 10000, 500)
// var shadowSumThresh = 0.5;
// var dilatePixels = 2;
// ver teste removindo uma parte da projeção 
// https://code.earthengine.google.com/5addf9b02520e27491e66fa3df2faba7
var cloudBand = 'cloudMask';
var cloudProject = function (image) {

    // Get the cloud mask
    var cloud = image.select('cloudMask');
    var cirrus = image.select('B10').multiply(0.0001);
    var cdi = image.select('cdi');
    
    // Assume low-to-mid atmospheric clouds to be pixels where probability
    // is greater than 65%, and CDI is less than -0.5. For higher atmosphere
    // cirrus clouds, assume the cirrus band is greater than 0.01.
    // The final cloud mask is one or both of these conditions.
    var isCloud = cloud.and(cdi.lt(-0.5)).or(cirrus.gt(0.01));

    // Reproject is required to perform spatial operations at 20m scale.
    // 20m scale is for speed, and assumes clouds don't require 10m precision.
    // isCloud = isCloud.focal_min(3).focal_max(6);
    isCloud = isCloud.reproject({crs: cdi.projection(), scale: 20});
  
    // Project shadows from clouds we found in the last step. This assumes 
    // we're working in a UTM projection.
    var azimuth_angle = 'MEAN_SOLAR_AZIMUTH_ANGLE';
    var shadowAzimuth = ee.Number(90).subtract(
                          ee.Number(image.get(azimuth_angle)));

    // With the following reproject, the shadows are projected 5km.
    var isCloud30 = isCloud.directionalDistanceTransform(shadowAzimuth, 50);
    isCloud30 = isCloud30.reproject({crs: cdi.projection(), scale: 10});
    print("isCoud30 ", isCloud30);
   
    isCloud30 = isCloud30.select("distance").unmask(0).gt(0)
    isCloud30 = isCloud30.focal_min(2).focal_max(6);
    print("isCoud30 ", isCloud30);
    // Map.addLayer(isCloud30, {min: 0, max: 1}, 'isCloud30');  //, palette: 'blue'
    cloud = cloud.focal_min(2).focal_max(6);
    var ccCloudShadown = isCloud30.add(cloud).gt(0).rename('mcloud_shadow'); //
    
    // return image; //.not()
    return image.addBands(ccCloudShadown).select(['mcloud_shadow']); //
};

var add_shadow_bands= function(img){
    var pmt_focal = {
        'radius': 1,
        'kernelType': 'square',
        'iterations': 1
    }
    // Identify water pixels from the SCL band.
    // not_water = img.select('SCL').neq(6)

    // # shadownMask = img.select('SCL').eq(2).multiply(
    // #             img.select('SCL').eq(7)).add(
    // #                 img.select('SCL').eq(8)).add(
    // #                     img.select('SCL').eq(9))

    // Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4
    var NIR_DRK_THRESH = 0.15
    var CLD_PRJ_DIST = 1
    print("valor SR_BAND_SCALE", SR_BAND_SCALE)
    // .multiply(not_water)
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH * SR_BAND_SCALE).rename('dark_pixels')
    var shadownMask = img.select("B9").lt(230).rename('shadowsB9')
    shadownMask = shadownMask.focalMax(pmt_focal);

    // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('cloudMask').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
                    .reproject({'crs': img.select(0).projection(), 'scale': 100})
                        .select('distance').mask().rename('cloud_transform'))

    // Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows')
    // shadows = shadows.add(shadownMask).gt(0)

    var mascara = img.select('clouds').add(shadownMask).gt(0)

    // Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([shadows, shadownMask]))
}

var applicateCeof = function(image){
                image = image.multiply(slope);
                image = image.add(intercept);
                return image;
          };
          
var maskingCoverCloud = function(image, idImage){
            var imgCloud = ee.Image(param.asset_S2Mask + idImage).gt(30).rename('cloudMask');
            var imgCloudShadow = add_shadow_bands(image.addBands(imgCloud));
            return image.addBands(imgCloudShadow.select(['cloudMask','shadows', 'shadowsB9']));//; .updateMask(imgCloudShadow.eq(0))
    };
          
          
var idIm1 = '/20160407T131202_20160407T193256_T23KMB';
// var idIm1 = '/20160424T130240_20160424T181611_T24LTM';
var idIm2 = '/20160424T130242_20160424T130240_T24LTM';
var idIm3 = '/20160623T130251_20160623T193957_T24LTM';
var Image1 = ee.Image(param.asset_S2harmonic + idIm1);//
var cdi = ee.Algorithms.Sentinel2.CDI(Image1);
Image1 = Image1.addBands(cdi.rename('cdi'));

// var Image2 = ee.Image(param.asset_S2harmonic + idIm2).select(param.bandVis);
// var Image3 = ee.Image(param.asset_S2harmonic + idIm3).select(param.bandVis);

print("show properties of image # 1", Image1);
// print("show properties of image # 2", Image2);
// print("show properties of image # 3", Image3);
Map.addLayer(Image1, vis.imgS2, "image # 1");

Image1 = maskingCoverCloud(Image1, idIm1);
// Image2 = removeCloud(Image2, idIm2);
// Image3 = removeCloud(Image3, idIm3);
print("show variaveis of image ", Image1);

var cloudMask = Image1.select("cloudMask")
var maskShadow = Image1.select("shadows");
var cloudShadown = Image1.select('shadowsB9')
// Map.addLayer(cloudShadown, {min: 0, max: 1, palette: 'green, yellow'}, 'shadownCloud'); //.selfMask()
Map.addLayer(cloudMask.selfMask(), {min: 0, max: 1, palette: 'blue'}, 'cloud');
Map.addLayer(maskShadow.selfMask(), {min: 0, max: 1}, 'shadown');
Map.addLayer(cloudShadown.selfMask(), {min: 0, max: 1, palette: 'yellow'}, 'shadownB9');

// Map.addLayer(Image2, vis.imgS2, "image # 2");
// Map.addLayer(Image3, vis.imgS2, "image # 3");

Image1 = Image1.select(param.bandVis);
var img1corr = applicateCeof(Image1);
// // var img2corr = applicateCeof(Image2);
// // var img3corr = applicateCeof(Image3);
// print("img correlate ", img1corr);

Map.addLayer(img1corr, vis.imgS2, "imgCorr # 1");
// Map.addLayer(img2corr, vis.imgS2, "imgCorr # 2");
// Map.addLayer(img3corr, vis.imgS2, "imgCorr # 3");





