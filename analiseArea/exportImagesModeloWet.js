var vis = {
    imagem: {
        min:0, 
        max: 1500,
        bands: ["B4","B3","B2"]
    },
    maskCC: {
        min:0, 
        max: 1,
        palette: '000000,FF0000'
    },
    mosaic_Median_Wet: {
        min:10, 
        max: 2500,
        bands: ["red_median_wet","green_median_wet","blue_median_wet"]
    },
};


var exportSHP_poligonos = function(featColAlerts, namesExp){
    var assetIdExp = 'users/mapbiomascaatinga04/poligonos_exc/' + namesExp
    var pmtroExpAsset = {
              "collection": featColAlerts,
              'description': namesExp,
              'assetId': assetIdExp            
    }
    
    Export.table.toAsset(pmtroExpAsset);
    print("tables SHP  exported ! ");
}
/**
 * Function to mask clouds using the Sentinel-2 QA band
 * @param {ee.Image} image Sentinel-2 image
 * @return {ee.Image} cloud masked Sentinel-2 image
 */
function maskS2clouds(image) {
    var qa = image.select('QA60');
  
    // Bits 10 and 11 are clouds and cirrus, respectively.
    var cloudBitMask = 1 << 10;
    var cirrusBitMask = 1 << 11;
  
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
        .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  
    return image.updateMask(mask).divide(10000);
}

var dictTiles = {
    '138': [
        "23KPA","23KPV","23KQB","23LNH","23LNJ","23LPJ","23LRH","23MQQ","23MRQ",
        "23MRR","23MRS","24MTA","24MTB","24MTV","24MUA","24MUB","24MUU","24MUV"
    ], 
    '38': [
        "23LLF","23LLG","23LLH","23LMF","23LMG","23LMJ","23LMK","23LNH","23LNJ",
        "23LNK","23LNL","23MPR","23MQS"
    ], 
    '52': [
        "24KUG","24KVG","24LUH","24LVH","24LVJ","24LVK","24LVL","24LVM","24LWL",
        "24LWM","24LWN","24LWP","24LWQ","24LXM","24LXN","24LXP","24LXQ","24LYP",
        "24LYQ","24LYR","24LZR","24MZS","24MZT","24MZU","24MZV","25MBP","25MBQ"
    ],
    '81': ["23MLM"],
    '95': [
        "23KPV","23KQA","23KQB","23KQV","23KRA","23KRB","23LQD","23LQE","23LRD",
        "23LRE","23LRF","23LRG","24KTG","24LTH","24LTJ","24LTK","24LTL","24LTM",
        "24LTN","24LTP","24LUH","24LUJ","24LUK","24LUL","24LUM","24LUN","24LUP",
        "24LVM","24LVN","24LVP","24MUS","24MUT","24MUU","24MVA","24MVB","24MVS",
        "24MVT","24MVU","24MWA","24MWB","24MWT","24MWU","24MWV","24MXA","24MXV"
    ],
    '9': [
        "24LZQ","24LZR","24MZS","24MZT","25LBL","25MBM","25MBN","25MBP"
    ]
}
var listOrb = [9,52,81,95,138];

var param = {
    'asset_geodatin_work': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/mosaics-CAATINGA-4',
    'asset_geodatin_next': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
};
var bandsS2 = ['B2','B3', 'B4', 'B8', 'B9', 'B11','B12']
var limiarCC = 40;
var orbitSelect = 95;
var tileSelect = 9;
var yearShow = 2016;
var dataInic = '-01-01';
var yearPos = yearShow + 1;
var data_inic = yearShow.toString() + dataInic;
var data_end = yearShow.toString() + '-06-01';
print(" Processing desde " + data_inic + " até " + data_end);

var geoGrade = ee.FeatureCollection(param['gradeS2Corr'])
                    .filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbitSelect))
                    .filter(ee.Filter.eq('MGRS_TILE', dictTiles[orbitSelect.toString()][tileSelect]))

print("numero de grades selectionadas ", geoGrade.size())

var imgColGeodatin = ee.ImageCollection(param.asset_geodatin_next).merge(
                          ee.ImageCollection(param.asset_geodatin_work))
                              .filter(ee.Filter.eq('year', yearShow))
                              .filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbitSelect.toString()))
                              .filter(ee.Filter.eq('MGRS_TILE', dictTiles[orbitSelect.toString()][tileSelect])); 
print("imagem Geodatin ",imgColGeodatin);
var datasetS2 = ee.ImageCollection(param.asset_S2harmonic)
                  .filterDate(data_inic , data_end)
                  .filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbitSelect))
                  .filter(ee.Filter.eq('MGRS_TILE', dictTiles[orbitSelect.toString()][tileSelect]))
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
                  .sort('CLOUDY_PIXEL_PERCENTAGE');
                  
var datasetMask = ee.ImageCollection(param.asset_S2Mask)
                  .filterDate(data_inic , data_end)
                  .filterBounds(geoGrade.geometry().centroid());

Map.addLayer(geoGrade, {color: 'green'}, "tile", false);
Map.addLayer(imgColGeodatin.first(), vis.mosaic_Median_Wet, "geodatin Median Wet");
var lstSystem = datasetS2.reduceColumns(
                                ee.Reducer.toList(2), 
                                ['system:index', 'CLOUDY_PIXEL_PERCENTAGE']
                            ).get('list').getInfo();
print("lista de parametros ", lstSystem);
var imgsSelect = ee.List([]);
var contador = 0;
lstSystem.forEach(function(pmtrosCC){
        var imS2 = datasetS2.filter(ee.Filter.eq('system:index', pmtrosCC[0])).first();
        if (contador < 1){
            print("imagem " + pmtrosCC[0] , imS2);
        }
        
        var maskImS2  = datasetMask.filter(
                        ee.Filter.eq('system:index', pmtrosCC[0])
                    ).first();
           
        imS2 = imS2.select(bandsS2)
        Map.addLayer(imS2, vis.imagem, pmtrosCC[0], false);
        Map.addLayer(ee.Image(maskImS2).lt(limiarCC).eq(0).selfMask(), vis.maskCC, 'mascara', false); //
        contador = contador + 1
})


// Paint all the polygon edges with the same number and width, display.
var outline = ee.Image().toByte().paint({
  featureCollection: geoGrade,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: 'FF0000'}, 'bordeGrade');

exportSHP_poligonos(poligono, '20161021T130242_20161021T192022_T23KRB');