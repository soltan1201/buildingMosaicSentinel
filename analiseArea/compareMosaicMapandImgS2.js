
var vis = {
    mosaico: {
        min:30, 
        max: 3500,
        bands: ["red_median","green_median","blue_median"]
    },
    imgS2: {
        min:30, 
        max: 3500,
        bands: ["B4","B3","B2"]
    }, 
    maskCC: {
        min: 0, 
        max: 1,
        palette: "000000,00FF00,FFFA00"
    }
}

var dictTiles = {
    '138': [
          "23LPH","23LPJ","23LRK","23MRS","24LTQ","24MUA","24MUB"
    ],
    '52': [
          "24LVJ","24LVK","24LVL","24LWQ","24LWR","24LXM","24LXN",
          "24LXQ","24LXR","24MYT","24MZT"
    ],
    // '95': [
    //       "23LRF","23LRG","24LTL","24LTM","24LTN","24LTQ","24LUL",
    //       "24LUM","24LUN","24LUP","24LVM","24LVN","24LVP","24LWR",
    //       "24MUA","24MUS","24MUT","24MUU","24MUV","24MVA","24MVU",
    //       "24MVV","24MWA","24MWT","24MWU"
    // ],
    '95': ['24LUL', '24LTL', '23LRG', '24LTM', '24LUM'],
    '9' : ["25MBM"]
}

var param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
}
var limiarCC = 30;
var tileSelect = 95;
var yearShow = 2016;
var dataInic = '-06-01';
var yearPos = yearShow + 1;
var data_inic = yearShow.toString() + dataInic;
var data_end = yearPos.toString() + '-01-01';
print(" Processing desde " + data_inic + " até " + data_end);

var limitCaat = ee.FeatureCollection(param.assetCaat).geometry();
var gradeS2 = ee.FeatureCollection(param.gradeS2Corr);
var mosaicMapB = ee.ImageCollection(param.asset_mapbiomas)
                    .filterBounds(limitCaat)
                    .filter(ee.Filter.eq('year', yearShow));
print("show mosaic ", mosaicMapB)                    
var gradesselect =  gradeS2.filterBounds(geometry);
var lstTiles = gradesselect.reduceColumns(
                              ee.Reducer.toList(2), 
                              ['SENSING_ORBIT_NUMBER','MGRS_TILE']).get('list');
                              
print(" list grades afetadas ", ee.List(lstTiles).distinct());

var datasetS2 = ee.ImageCollection(param.asset_S2harmonic)
                  .filterDate(data_inic , data_end)
                  .filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', tileSelect))
                  .filter(ee.Filter.inList('MGRS_TILE', dictTiles[tileSelect.toString()]))
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60));
                  
var datasetMask = ee.ImageCollection(param.asset_S2Mask)
                  .filterDate(data_inic , data_end);

                 
print("DataSet Sentinel 2 ", datasetS2);
print("Dataset Mask Sentinel 2", datasetMask.limit(3));
var maskImS2 = null;
var imS2 = null;


Map.addLayer(gradeS2, {color: 'yellow'}, 'grades', false)
Map.addLayer(limitCaat, {color: 'green'}, 'limitCaat', false);

Map.addLayer(mosaicMapB, vis.mosaico, 'mosaic'); 

var lstSystem = datasetS2.reduceColumns(
                      ee.Reducer.toList(2), 
                      ['system:index', 'CLOUDY_PIXEL_PERCENTAGE']
                    ).get('list').getInfo();
print("lista de parametros ", lstSystem);
lstSystem.forEach(function(pmtrosCC){
    imS2 = datasetS2.filter(ee.Filter.eq('system:index', pmtrosCC[0])).first();
    if (parseInt(pmtrosCC[1]) > 1){
        maskImS2  = datasetMask.filter(
                                    ee.Filter.eq('system:index', pmtrosCC[0])
                                ).first();
        print("bandas mask " + pmtrosCC[0] , maskImS2);
        imS2 = imS2.updateMask(maskImS2.lt(limiarCC));
    }
    Map.addLayer(imS2, vis.imgS2, pmtrosCC[0], false);
    Map.addLayer(maskImS2.lt(limiarCC).eq(0).selfMask(), vis.maskCC, 'mascara', false); //
    
})
// Map.addLayer(maskImS2, {min:0, max:100}, 'mascaraOrin', false); //
// 

