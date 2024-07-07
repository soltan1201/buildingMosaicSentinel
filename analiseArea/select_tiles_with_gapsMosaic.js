//https://code.earthengine.google.com/081493d9ee38e50a4f229540a6ca1eed
var vis = {
    mosaicQM_Wet: {
        min:0, 
        max: 1500,
        bands: ["red_QM_wet","green_QM_wet","blue_QM_wet"]
    },
    mosaic_QM_Dry: {
        min:0, 
        max: 1500,
        bands: ["red_QM_dry","green_QM_dry","blue_QM_dry"]
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
    }
};
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
var listTiles = [
    "23KPA","23KPV","23KQB","23LNH","23LNJ","23LPJ","23LRH","23MQQ","23MRQ",
    "23MRR","23MRS","24MTA","24MTB","24MTV","24MUA","24MUB","24MUU","24MUV",
    "23LLF","23LLG","23LLH","23LMF","23LMG","23LMJ","23LMK","23LNH","23LNJ",
    "23LNK","23LNL","23MPR","23MQS",
    "24KUG","24KVG","24LUH","24LVH","24LVJ","24LVK","24LVL","24LVM","24LWL",
    "24LWM","24LWN","24LWP","24LWQ","24LXM","24LXN","24LXP","24LXQ","24LYP",
    "24LYQ","24LYR","24LZR","24MZS","24MZT","24MZU","24MZV","25MBP","25MBQ",
    "23MLM",
    "23KPV","23KQA","23KQB","23KQV","23KRA","23KRB","23LQD","23LQE","23LRD",
    "23LRE","23LRF","23LRG","24KTG","24LTH","24LTJ","24LTK","24LTL","24LTM",
    "24LTN","24LTP","24LUH","24LUJ","24LUK","24LUL","24LUM","24LUN","24LUP",
    "24LVM","24LVN","24LVP","24MUS","24MUT","24MUU","24MVA","24MVB","24MVS",
    "24MVT","24MVU","24MWA","24MWB","24MWT","24MWU","24MWV","24MXA","24MXV",
    "24LZQ","24LZR","24MZS","24MZT","25LBL","25MBM","25MBN","25MBP"
];

var param = {
    'asset_geodatin_work': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/mosaics-CAATINGA-4',
    'asset_geodatin_next': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
};
var geoGrade = ee.FeatureCollection(param['gradeS2Corr'])

var gradesSelect = geoGrade.filterBounds(pointCorregWet);
print("numero de grades selectionadas ", gradesSelect.size())
var lstGradesprop = gradesSelect.reduceColumns(
                            ee.Reducer.toList(2), 
                            ['SENSING_ORBIT_NUMBER', 'MGRS_TILE']).get('list');
print(" lista de propiedades de TILES e ORBITAS ", lstGradesprop);

var limiarCC = 40;
var tileSelect = 95;
var yearShow = 2016;
var dataInic = '-06-01';
var yearPos = yearShow + 1;
var data_inic = yearShow.toString() + dataInic;
var data_end = yearPos.toString() + '-01-01';
print(" Processing desde " + data_inic + " até " + data_end);


var imgColGeodatin = ee.ImageCollection(param.asset_geodatin_next).merge(
                          ee.ImageCollection(param.asset_geodatin_work))
                              .filter(ee.Filter.eq('year', yearShow)); 
print("show imagens geodatin ", imgColGeodatin);

var datasetS2 = ee.ImageCollection(param.asset_S2harmonic)
                  .filterDate(data_inic , data_end)
                  .filter(ee.Filter.inList('SENSING_ORBIT_NUMBER', listOrb))
                  .filter(ee.Filter.inList('MGRS_TILE', listTiles))
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30));
                  
var datasetMask = ee.ImageCollection(param.asset_S2Mask)
                  .filterDate(data_inic , data_end);

Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Wet, "geodatin Median Wet");
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Dry, "geodatin Median Dry");
Map.addLayer(geoGrade, {color: 'green'}, "tile", false);

var lstSystem = datasetS2.reduceColumns(
                                ee.Reducer.toList(2), 
                                ['system:index', 'CLOUDY_PIXEL_PERCENTAGE']
                            ).get('list').getInfo();
print("lista de parametros ", lstSystem);
var imgsSelect = ee.List([]);
lstSystem.forEach(function(pmtrosCC){
        var imS2 = datasetS2.filter(ee.Filter.eq('system:index', pmtrosCC[0])).first();
        if (parseInt(pmtrosCC[1]) > 5){
            var maskImS2  = datasetMask.filter(
                            ee.Filter.eq('system:index', pmtrosCC[0])
                        ).first();
            // print("bandas mask " + pmtrosCC[0] , maskImS2);
            imS2 = imS2.updateMask(maskImS2.lt(limiarCC));
        }
        // Map.addLayer(imS2, vis.imgS2, pmtrosCC[0], true);
        // Map.addLayer(ee.Image(maskImS2).lt(limiarCC).eq(0).selfMask(), vis.maskCC, 'mascara', false); //
})
imgsSelect = ee.ImageCollection.fromImages(imgsSelect)
Map.addLayer(imgsSelect, vis.imgS2, 'imagesSelect', true);

// Paint all the polygon edges with the same number and width, display.
var outline = ee.Image().toByte().paint({
  featureCollection: geoGrade,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: 'FF0000'}, 'bordeGrade');