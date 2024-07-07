
/*/////////////////////////////////////////////////////////////////////////////
DESCRIPTION:
This script accomplishes several tasks:
 1) Loads study region 
 2) Gathers Landsat 4,5,7,8 imagery for a specified season (see below)
 3) Masks shadows and clouds 
 4) Exports composite to asset
 
*//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// STEP 0: Define user inputs:

// 1. Specify study area: Study area
var StudyRegion = omo;

// 2. Update the startJulian and endJulian variables to indicate your seasonal constraints. 
// This supports wrapping for tropics and southern hemisphere.
// startJulian: Starting Julian date 
// endJulian: Ending Julian date
var startJulian = 306;
var endJulian = 61;     

// 3. Specify the year for which you would like to make a composite.
// Add in an annual buffer to include imagery from the same season 
// timeframe from the prior and following year. 
var year = 2015;
var timebuffer = 2;

// 4. Set up Names for the export
var SRname = 'omo_fr'; 
var exportPathRoot = 'users/Kolawoleayobami72';

// Metadata tag
var versionNumber = 1;

// 5. Cloud and cloud shadow masking parameters.
// cloudThresh: If using the cloudScoreTDOMShift method-Threshold for cloud 
//    masking (lower number masks more clouds.  Between 10 and 30 generally 
//    works best)
var cloudThresh = 20;

// dilatePixels: Number of pixels to buffer clouds and cloud 
//    shadows by (1 or 2 generally is sufficient)
var dilatePixels = 2;
    
// cloudHeights: Height of clouds to use to project cloud shadows
var cloudHeights = ee.List.sequence(200,5000,500);

// zScoreThresh: Threshold for cloud shadow masking- lower number masks out 
//    less.  Between -0.8 and -1.2 generally works well
var zScoreThresh = -0.8;
    
// shadowSumThresh: Sum of IR bands to include as shadows within TDOM and the 
//    shadow shift method (lower number masks out less)
var shadowSumThresh = 0.35;
var useCloudProject = true;
    
// metadataCloudCoverMax: Cloud cover percentage in image metadata threshold.
//    Will not include images with cloud cover > this number. Set to 100 if 
//    using metadata cloud cover is not wanted
var metadataCloudCoverMax = 100;

// 6. Adjust the map visualization parameters.
// vizParams: Options for map visualization of exported data set.
var vizParams = {
  bands: ['nir', 'swir1', 'red'],
  min: 500,
  max: 5000,
  gamma: [1.6]
  };
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// End User defined variables.
// Change items below with caution.
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// Prepare dates
var startYear = year - timebuffer;
var endYear = year + timebuffer;

// Prepare dates, cont.
if(startJulian > endJulian){endJulian = endJulian + 365}
var startDate = ee.Date.fromYMD(startYear,1,1).advance(startJulian-1,'day');
var endDate = ee.Date.fromYMD(endYear,1,1).advance(endJulian-1,'day');
print('Start and end dates:',startDate,endDate);

// set up export name and path
var exportName = SRname + 'Landsat_TOA_' + startJulian + '_' + endJulian + 
  '_' + startYear + '_' + endYear;
var exportPath = exportPathRoot + '/' + exportName;
print('Write down the Asset ID:', exportPath);

// STEP 1: Load Mekong study region (Myanmar, Thailand, Laos, Vietnam, Cambodia)
// var mekongBuffer = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ');
// var mekongRegion = mekongBuffer.geometry();
var studyArea = StudyRegion.geometry();

// STEP 2: Get Landsat 4,5,8,7 Image Collections
var lsTOA = getImageCollection(studyArea,startDate,endDate,startJulian,endJulian);

//Map.addLayer(ls.median(),vizParams,'Landsat before masking composite',false);

// STEP 3: Compute a cloud score and mask clouds
lsTOA = lsTOA.map(landsatCloudScore);

// Find and mask out dark outliers
lsTOA = simpleTDOM2(lsTOA,zScoreThresh,shadowSumThresh,dilatePixels);

if (useCloudProject) {
  // Run cloud project to get the final cloud and shadow masked image collection
  lsTOA = lsTOA.map(function(img){return cloudProject(img,shadowSumThresh,dilatePixels,
       cloudHeights,'SUN_AZIMUTH','SUN_ELEVATION')});
} else {
  // Just apply the TDOM Mask
  lsTOA = lsTOA.map(function(img){return img.updateMask(img.select(['TDOMMask']))});
}

// STEP 4: Add additional bands
// Add common spectral indices
lsTOA = lsTOA.map(addIndices);

// Add tasseled cap transformation, tasseled cap angles, and NDSV
var tcInputBands = ee.List(['blue','green','red','nir','swir1','swir2']);
lsTOA = lsTOA.map(function(img){
  img = getTasseledCap(img,tcInputBands);
  img = addTCAngles(img);
  return img;
});

// Get post cloud-masking counts
// var countTOA_post = lsTOA.select(['blue']).count().rename('count');

// STEP 5: Prepare composite for export
// Get median bands
///////////////////////////////////////////////////////////////////////////////
var medianBands = ee.List(['blue','green','red','nir','swir1','swir2' //,'temp'
  //, 'ND_blue_green','ND_blue_red','ND_blue_nir','ND_blue_swir1','ND_blue_swir2'
  //, 'ND_green_red','ND_green_nir','ND_green_swir1','ND_green_swir2','ND_red_swir1'
  //, 'ND_red_swir2','ND_nir_red','ND_nir_swir1','ND_nir_swir2','ND_swir1_swir2'
  //, 'R_swir1_nir','R_red_swir1','EVI','SAVI','IBI','brightness','greenness'
  //, 'wetness','fourth', 'fifth', 'sixth','tcAngleBG','tcDistBG','tcAngleGW'
  //, 'tcDistGW','tcAngleBW','tcDistBW'
  ]);
var medianCompositeTOA = lsTOA.select(medianBands).median();

// Get standard deviation bands
var stdDevBands = ee.List(['blue','green','red','nir','swir1','swir2'] 
  //, 'temp', 'ND_nir_red', 'ND_nir_swir2', 'ND_green_swir1'
  );
var stdDevCompositeTOA = lsTOA.select(stdDevBands).reduce(ee.Reducer.stdDev());

// Combine all bands with mask and count bands
var compositeTOA = medianCompositeTOA.addBands(stdDevCompositeTOA);//.addBands(countTOA_post);

// STEP 6: Add data layers into the map window. 
// Display the Landsat Composite.
///////////////////////////////////////////////////////////////////////////////
// vizParams: Options for map visualization
var vizParams = {
  bands: ['nir', 'swir1', 'red'],
  min: 500,
  max: 5000,
  gamma: [1.6]
  };
    
Map.addLayer(compositeTOA, vizParams,'Landsat Composite TOA', true);

// Load the study region, with a blue outline.
// Create an empty image into which to paint the features, cast to byte.
// Paint all the polygon edges with the same number and width, display.
var empty = ee.Image().byte();
var outline = empty.paint({
  featureCollection: studyArea,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: '0000FF'}, "Study Area", true);

// Center the map window on the study region.
Map.centerObject(studyArea, 6);

// STEP 7: Reformat data for export. 
///////////////////////////////////////////////////////////////////////////////

var compositeBands = compositeTOA.bandNames();
var nonDivideBands = ee.List(['temp','temp_stdDev','count','SR_mask','R_swir1_nir','R_red_swir1']);
var composite10000 = compositeTOA.select(compositeBands.removeAll(nonDivideBands)).multiply(10000);
// var composite100 = compositeTOA.select(['temp_stdDev','R_swir1_nir','R_red_swir1'],
//   ['thermal_stdDev','R_swir1_nir','R_red_swir1']).multiply(100);
// var composite10 = compositeTOA.select(['temp'],['thermal']).multiply(10);
// var composite1 = compositeTOA.select(['count']);
compositeTOA = composite10000.int16();//.addBands(composite100).addBands(composite10).addBands(composite1)

// Add metadata, cast to integer, and export composite
compositeTOA = compositeTOA.set({
  'system:time_start': ee.Date.fromYMD(year,6,1).millis(),
  'date': ee.Date.fromYMD(year,6,1),
  'source':'TOA',
  'version': versionNumber});

print(compositeTOA);

// STEP 7: Export the composite. 
///////////////////////////////////////////////////////////////////////////////

asyncExportToAssetWrapper(compositeTOA,exportName,exportPath,'mean',
  studyArea,30,'EPSG:4326');

// STEP 8: Display re-formatted Landsat imagery. 
///////////////////////////////////////////////////////////////////////////////

Map.addLayer(compositeTOA, vizParams2,'Exported Landsat Composite TOA', false);

///////////////////////////////////////////////////////////////////////////////
//FUNCTIONS
///////////////////////////////////////////////////////////////////////////////
//Function for acquiring Landsat TOA image collection
function getImageCollection(studyArea,startDate,endDate,startJulian,endJulian){
  var ls;var l4TOAs;var l5TOAs;var l7TOAs;var l8TOAs;var out;
  
  var sensorBandDictLandsatTOA =ee.Dictionary({L8 : ee.List([1,2,3,4,5,9,6]),
                        L7 : ee.List([0,1,2,3,4,5,7]),
                        L5 : ee.List([0,1,2,3,4,5,6]),
                        L4 : ee.List([0,1,2,3,4,5,6])
  });
  var bandNamesLandsatTOA = ee.List(['blue','green','red','nir','swir1','temp',
      'swir2']);

  l4TOAs = ee.ImageCollection('LANDSAT/LT4_L1T_TOA')
      .filterDate(startDate,endDate)
      .filter(ee.Filter.calendarRange(startJulian,endJulian))
      .filterBounds(studyArea)
      .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
      .select(sensorBandDictLandsatTOA.get('L4'),bandNamesLandsatTOA);
  
  l5TOAs = ee.ImageCollection('LANDSAT/LT5_L1T_TOA')
      .filterDate(startDate,endDate)
      .filter(ee.Filter.calendarRange(startJulian,endJulian))
      .filterBounds(studyArea)
      .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
      .select(sensorBandDictLandsatTOA.get('L5'),bandNamesLandsatTOA);
  
  l8TOAs = ee.ImageCollection('LANDSAT/LC8_L1T_TOA')
      .filterDate(startDate,endDate)
      .filter(ee.Filter.calendarRange(startJulian,endJulian))
      .filterBounds(studyArea)
      .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
      .select(sensorBandDictLandsatTOA.get('L8'),bandNamesLandsatTOA);

  l7TOAs = ee.ImageCollection('LANDSAT/LE7_L1T_TOA')
      .filterDate(startDate,endDate)
      .filter(ee.Filter.calendarRange(startJulian,endJulian))
      .filterBounds(studyArea)
      .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
      .select(sensorBandDictLandsatTOA.get('L7'),bandNamesLandsatTOA);
  
  ls = ee.ImageCollection(l4TOAs.merge(l5TOAs).merge(l7TOAs).merge(l8TOAs));
  out = ls;
  return out;
}

///////////////////////////////////////////////////////////////////////////////
// Function to merge images from the same date but different collections
function joinCollections(c1,c2){
  // Define an inner join.
  var innerJoin = ee.Join.inner();
  
  // Specify an equals filter for image timestamps.
  var filterTimeEq = ee.Filter.equals({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });
  
  // Apply the join.
  var innerJoined = innerJoin.apply(c1, c2, filterTimeEq);
  
  
  // Map a function to merge the results in the output FeatureCollection.
  var joined = innerJoined.map(function(feature) {
    return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
  });
  return ee.ImageCollection(joined);
}

///////////////////////////////////////////////////////////////////////////////
// A helper to apply an expression and linearly rescale the output.
// Used in the landsatCloudScore function below.
function rescale(img, exp, thresholds) {
  return img.expression(exp, {img: img})
      .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
}

///////////////////////////////////////////////////////////////////////////////
// Compute a cloud score and adds a band that represents the cloud mask.  
// This expects the input image to have the common band names: 
// ["red", "blue", etc], so it can work across sensors.
function landsatCloudScore(img) {
  // Compute several indicators of cloudiness and take the minimum of them.
  var score = ee.Image(1.0);
  // Clouds are reasonably bright in the blue band.
  score = score.min(rescale(img, 'img.blue', [0.1, 0.3]));
 
  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]));
   
  // Clouds are reasonably bright in all infrared bands.
  score = score.min(
      rescale(img, 'img.nir + img.swir1 + img.swir2', [0.3, 0.8]));

  // Clouds are reasonably cool in temperature.
  score = score.min(rescale(img,'img.temp', [300, 290]));

  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['green', 'swir1']);
  score =  score.min(rescale(ndsi, 'img', [0.8, 0.6])).multiply(100).byte();
  score = score.lt(cloudThresh).rename('cloudMask');
  img = img.updateMask(score);
  return img.addBands(score);
}

///////////////////////////////////////////////////////////////////////////////
//Function for finding dark outliers in time series.
//Original concept written by Carson Stam and adapted by Ian Housman.
//Adds a band that is a mask of pixels that are dark, and dark outliers.
function simpleTDOM2(collection,zScoreThresh,shadowSumThresh,dilatePixels){
  var shadowSumBands = ['nir','swir1'];
  
  //Get some pixel-wise stats for the time series
  var irStdDev = collection.select(shadowSumBands).reduce(ee.Reducer.stdDev());
  var irMean = collection.select(shadowSumBands).mean();
  
  //Mask out dark dark outliers
  collection = collection.map(function(img){
    var zScore = img.select(shadowSumBands).subtract(irMean).divide(irStdDev);
    var irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum());
    var TDOMMask = zScore.lt(zScoreThresh).reduce(ee.Reducer.sum()).eq(2)
        .and(irSum.lt(shadowSumThresh)).not();
    TDOMMask = TDOMMask.focal_min(dilatePixels);
    return img.addBands(TDOMMask.rename('TDOMMask'));
  });
  
  return collection;
}

///////////////////////////////////////////////////////////////////////////////
//Function for wrapping cloud and shadow masking together.
//Assumes image has cloud mask band called "cloudMask" and a TDOM mask called 
//"TDOMMask".
function cloudProject(img,shadowSumThresh,dilatePixels,cloudHeights,
    azimuthField,zenithField){
    
    //Get the cloud mask
    var cloud = img.select('cloudMask').not();
    cloud = cloud.focal_max(dilatePixels);
    cloud = cloud.updateMask(cloud);
    
    //Get TDOM mask
    var TDOMMask = img.select(['TDOMMask']).not();
    
    //Project the shadow finding pixels inside the TDOM mask that are dark and 
    //inside the expected area given the solar geometry
    //Find dark pixels
    var darkPixels = img.select(['nir','swir1','swir2'])
        .reduce(ee.Reducer.sum()).lt(shadowSumThresh);//.gte(1);
    
    //Get scale of image
    var nominalScale = cloud.projection().nominalScale();

    //Find where cloud shadows should be based on solar geometry
    //Convert to radians
    var meanAzimuth = img.get(azimuthField);
    var meanZenith = img.get(zenithField);
    var azR = ee.Number(meanAzimuth).multiply(Math.PI).divide(180.0)
        .add(ee.Number(0.5).multiply(Math.PI ));
    var zenR = ee.Number(0.5).multiply(Math.PI )
        .subtract(ee.Number(meanZenith).multiply(Math.PI).divide(180.0));
    
    //Find the shadows
    var shadows = cloudHeights.map(function(cloudHeight){
        cloudHeight = ee.Number(cloudHeight);
        var shadowCastedDistance = zenR.tan()
            .multiply(cloudHeight);//Distance shadow is cast
        var x = azR.cos().multiply(shadowCastedDistance)
            .divide(nominalScale).round();//X distance of shadow
        var y = azR.sin().multiply(shadowCastedDistance)
            .divide(nominalScale).round();//Y distance of shadow
        return cloud.changeProj(cloud.projection(), cloud.projection()
            .translate(x, y));
    });

    var shadow = ee.ImageCollection.fromImages(shadows).max();
 
  //Create shadow mask
  shadow = shadow.updateMask(cloud.mask().not());
  shadow = shadow.focal_max(dilatePixels);
  shadow = shadow.updateMask(darkPixels.and(TDOMMask));

  //Combine the cloud and shadow masks
  var combinedMask = cloud.mask().or(shadow.mask()).eq(0);
  
  //Update the image's mask and return the image
  img = img.updateMask(combinedMask);
  img = img.addBands(combinedMask.rename(['cloudShadowMask']));
  return img;
}

///////////////////////////////////////////////////////////////////////////////
// Function to add common (and less common) spectral indices to an image.
// Includes the Normalized Difference Spectral Vector from (Angiuli and Trianni, 2014)
function addIndices(img){
  // Add Normalized Difference Spectral Vector (NDSV)
  img = img.addBands(img.normalizedDifference(['blue','green']).rename('ND_blue_green'));
  img = img.addBands(img.normalizedDifference(['blue','red']).rename('ND_blue_red'));
  img = img.addBands(img.normalizedDifference(['blue','nir']).rename('ND_blue_nir'));
  img = img.addBands(img.normalizedDifference(['blue','swir1']).rename('ND_blue_swir1'));
  img = img.addBands(img.normalizedDifference(['blue','swir2']).rename('ND_blue_swir2'));

  img = img.addBands(img.normalizedDifference(['green','red']).rename('ND_green_red'));
  img = img.addBands(img.normalizedDifference(['green','nir']).rename('ND_green_nir')); //NDWBI
  img = img.addBands(img.normalizedDifference(['green','swir1']).rename('ND_green_swir1')); //NDSI, MNDWI
  img = img.addBands(img.normalizedDifference(['green','swir2']).rename('ND_green_swir2'));

  img = img.addBands(img.normalizedDifference(['red','swir1']).rename('ND_red_swir1'));
  img = img.addBands(img.normalizedDifference(['red','swir2']).rename('ND_red_swir2'));

  img = img.addBands(img.normalizedDifference(['nir','red']).rename('ND_nir_red')); //NDVI
  img = img.addBands(img.normalizedDifference(['nir','swir1']).rename('ND_nir_swir1')); //NDWI, LSWI, -NDBI
  img = img.addBands(img.normalizedDifference(['nir','swir2']).rename('ND_nir_swir2')); //NBR, MNDVI

  img = img.addBands(img.normalizedDifference(['swir1','swir2']).rename('ND_swir1_swir2'));
  
  // Add ratios
  img = img.addBands(img.select('swir1').divide(img.select('nir')).rename('R_swir1_nir')); //ratio 5/4
  img = img.addBands(img.select('red').divide(img.select('swir1')).rename('R_red_swir1')); // ratio 3/5

  // Add Enhanced Vegetation Index (EVI)
  var evi = img.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': img.select('nir'),
      'RED': img.select('red'),
      'BLUE': img.select('blue')
  }).float();
  img = img.addBands(evi.rename('EVI'));
  
  // Add Soil Adjust Vegetation Index (SAVI)
  // using L = 0.5;
  var savi = img.expression(
    '(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
      'NIR': img.select('nir'),
      'RED': img.select('red')
  }).float();
  img = img.addBands(savi.rename('SAVI'));
  
  // Add Index-Based Built-Up Index (IBI)
  var ibi_a = img.expression(
    '2*SWIR1/(SWIR1 + NIR)', {
      'SWIR1': img.select('swir1'),
      'NIR': img.select('nir')
    }).rename('IBI_A');
  var ibi_b = img.expression(
    '(NIR/(NIR + RED)) + (GREEN/(GREEN + SWIR1))', {
      'NIR': img.select('nir'),
      'RED': img.select('red'),
      'GREEN': img.select('green'),
      'SWIR1': img.select('swir1')
    }).rename('IBI_B');
  ibi_a = ibi_a.addBands(ibi_b);
  var ibi = ibi_a.normalizedDifference(['IBI_A','IBI_B']);
  img = img.addBands(ibi.rename('IBI'));
  
  return img;
}

///////////////////////////////////////////////////////////////////////////////
// Function to compute the Tasseled Cap transformation and return an image
// with the following bands added: ['brightness', 'greenness', 'wetness', 
// 'fourth', 'fifth', 'sixth']
function getTasseledCap(image,bands) {
  // Kauth-Thomas coefficients for Thematic Mapper data
  var coefficients = ee.Array([
    [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
    [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
    [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
    [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
    [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
    [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
  ]);
  // Make an Array Image, with a 1-D Array per pixel.
  var arrayImage1D = image.select(bands).toArray();
  
  // Make an Array Image with a 2-D Array per pixel, 6x1.
  var arrayImage2D = arrayImage1D.toArray(1);
  
  var componentsImage = ee.Image(coefficients)
    .matrixMultiply(arrayImage2D)
    // Get rid of the extra dimensions.
    .arrayProject([0])
    // Get a multi-band image with TC-named bands.
    .arrayFlatten(
      [['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']])
    .float();
  
  return image.addBands(componentsImage);
}

// Function to add Tasseled Cap angles and distances to an image.
// Assumes image has bands: 'brightness', 'greenness', and 'wetness'.
function addTCAngles(image){
  // Select brightness, greenness, and wetness bands
  var brightness = image.select(['brightness']);
  var greenness = image.select(['greenness']);
  var wetness = image.select(['wetness']);
  
  // Calculate Tasseled Cap angles and distances
  var tcAngleBG = brightness.atan2(greenness).divide(Math.PI).rename('tcAngleBG');
  var tcAngleGW = greenness.atan2(wetness).divide(Math.PI).rename('tcAngleGW');
  var tcAngleBW = brightness.atan2(wetness).divide(Math.PI).rename('tcAngleBW');
  var tcDistBG = brightness.hypot(greenness).rename('tcDistBG');
  var tcDistGW = greenness.hypot(wetness).rename('tcDistGW');
  var tcDistBW = brightness.hypot(wetness).rename('tcDistBW');
  image = image.addBands(tcAngleBG).addBands(tcAngleGW)
    .addBands(tcAngleBW).addBands(tcDistBG).addBands(tcDistGW)
    .addBands(tcDistBW);
  return image;
}

// Function to create a quality mosaic based on absolute difference from a percentile
function addAbsDiff (inCollection, qualityBand, percentile, sign){
  var bestQuality = inCollection.select([qualityBand]).reduce(ee.Reducer.percentile([percentile]));
  var out = inCollection.map(function(image) {
    var delta = image.select([qualityBand]).subtract(bestQuality).abs().multiply(sign);
    return image.addBands(delta.select([0], ['delta']));
  });
  return out;
}

function customQualityMosaic(inCollection,qualityBand,percentile){
  var inCollectionDelta = addAbsDiff(inCollection, qualityBand, percentile,-1);
  return inCollectionDelta.qualityMosaic('delta');
}

///////////////////////////////////////////////////////////////////////////////
// Function to export a provided image to an EE asset
function asyncExportToAssetWrapper(
  imageForExport,assetName,assetPath,pyramidingPolicy,roi,scale,crs){
  //Make sure image is clipped to roi in case it's a multi-part polygon
  imageForExport = imageForExport.clip(roi);
  assetName = assetName.replace(/\s+/g,'-');//Get rid of any spaces
  
  //Asynchronous approach to gathering converting server-side vectors to 
  //client-side without locking the browser
  roi.evaluate(function(roiInfo){
    var roiType = roiInfo.type.toString();
    //If it is a Polygon geometry...
    if( roiType === 'Polygon'){
      roi.bounds(1000).evaluate(function(polygonInfo){
        var region = polygonInfo.coordinates[0];
        Export.image.toAsset(imageForExport, assetName, assetPath, 
        {'.default':pyramidingPolicy}, null, region, scale, crs, null, 1e13);
      });
    }
    //If it is a MultiPolygon gometry.....
    else if( roiType === 'MultiPolygon'){
      roi.bounds(1000).evaluate(function(multiPolygonInfo){
        var region = multiPolygonInfo.coordinates[0];
        Export.image.toAsset(imageForExport, assetName, assetPath, 
        {'.default':pyramidingPolicy}, null, region, scale, crs, null, 1e13);
      });
    }
    //If it is a FeatureCollection.....
    else if( roiType === 'FeatureCollection'){
      roi.geometry(1000).bounds(1000).evaluate(function(featureCollectionInfo){
        var region = featureCollectionInfo.coordinates[0];
        Export.image.toAsset(imageForExport, assetName, assetPath, 
        {'.default':pyramidingPolicy}, null, region, scale, crs, null, 1e13);
      });
    }
    //Alert user if not supported
    else(
    alert('Type of feature is "'+roiType+ '". This is not handled\nIf a ' + 
    'Feature, can manually cast to featureCollections by using: ' + 
    'ee.FeatureCollection([myFeature])')
    );
  });
}