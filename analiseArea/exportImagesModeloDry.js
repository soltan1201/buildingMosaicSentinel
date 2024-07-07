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
    mosaic_Median_Dry: {
        min:10, 
        max: 2500,
        bands: ["red_median_dry","green_median_dry","blue_median_dry"]
    },
    planet: {
        'min':14,
        'max':3454,
        // 'gamma':1.8,
        'bands': ['R','G','B']
    }
};
// ver aqui  https://code.earthengine.google.com/ec95c96167e74c1055af0c8aa2a642f4

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
// tiles seleccionados por aqui 
//https://code.earthengine.google.com/1c3e9dd507bb2d904c8c212329353a9b
var dictTiles = {
      '95': [
              "23KQA","23KQB","23KQV","23KRA","23KRB","23LQC","23LRC",
              "23LRD","23LRF","23LRG","23LRK","24KTG","24LTH","24LTJ",
              "24LTL","24LTM","24LTN","24LTQ","24LUH","24LUJ","24LUK",
              "24LUM","24LUN","24LUQ","24LVJ","24LVN","24LVP","24LVQ",
              "24LVR","24LWR","24MTS","24MUA","24MUB","24MVA","24MVB",
              "24MWA","24MWV","24LUR","24LVR","24LWR"
            ],
            '138': [
                "23KPA","23KPV","23LND","23LQJ","23LRK","23MQS","23MRS",
                "24LTQ","24MTB","24MTS","24MTT","24MUA","24MUB","24MUU",
                "24MUV"
            ], 
            '38': ["23MNN","23MPN"],
            '52': [
                "24KUG","24KVG","24LUJ","24LVH","24LVJ","24LVK","24LVL",
                "24LVM","24LWH","24LWK","24LWL","24LWM","24LWQ","24LWR",
                "24LXM","24LXQ","24LYP","24MXV","24MYU","24MYV","24MZU",
                "24MZV","25MBQ","24LXP","24LXQ","24LYP","24LYQ"
            ],
            '9':[
                "24LZQ","24LZR","24MZS","24MZT","24MZU","24MZV","25LBL",
                "25MBM","25MBN","25MBP","25MBQ","24LZP"
            ]
}
var listOrb = [9,52,28,95,138];

var param = {
    'asset_geodatin_work': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/mosaics-CAATINGA-4',
    'asset_geodatin_next': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
     'assetNICFI': 'projects/planet-nicfi/assets/basemaps/americas',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
};
var bandsS2 = ['B2','B3', 'B4', 'B8', 'B9', 'B11','B12']
var limiarCC = 40;
var orbitSelect = 95;
var tileSelect = 3;
var yearShow = 2017;
var dataInic = '-06-01';
var yearPos = yearShow + 1;
var data_inic = yearShow.toString() + dataInic;
var data_end = yearPos.toString() + '-01-01';
print(" Processing desde " + data_inic + " até " + data_end);

var geoGrade = ee.FeatureCollection(param['gradeS2Corr'])

var knowCruza = geoGrade.filterBounds(ponto)
print(knowCruza)

geoGrade = geoGrade.filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbitSelect))
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
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
                  .sort('CLOUDY_PIXEL_PERCENTAGE');
                  
var datasetMask = ee.ImageCollection(param.asset_S2Mask)
                  .filterDate(data_inic , data_end)
                  .filterBounds(geoGrade.geometry().centroid());
                  
// This collection is not publicly accessible. To sign up for access,
// please see https://developers.planet.com/docs/integrations/gee/nicfi
var nicfi = ee.ImageCollection(param.assetNICFI);
// Filter basemaps by date and get the first image from filtered results
var basemap2016 = nicfi.filter(ee.Filter.date('2017-08-01','2018-01-01')).first().clip(geoGrade);
                  
                  
                  

Map.addLayer(geoGrade, {color: 'green'}, "tile", false);
Map.addLayer(basemap2016, vis.planet, 'Planet 2016');
Map.addLayer(imgColGeodatin.first(), vis.mosaic_Median_Dry, "geodatin Median Dry");
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

exportSHP_poligonos(poligono, "20160809T125311_20160809T192726_T24LWM");