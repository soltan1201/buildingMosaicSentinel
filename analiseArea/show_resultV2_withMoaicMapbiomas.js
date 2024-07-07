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

var param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'asset_geodatin_work': 'projects/mapbiomas-workspace/AMOSTRAS/col8/CAATINGA/mosaics-CAATINGA-4',
    'asset_geodatin_next': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
};
var geoGrade = ee.FeatureCollection(param['gradeS2Corr'])

var yearEst = 2016;
var imgColMapBiomas = ee.ImageCollection(param.asset_mapbiomas)
                        .filter(ee.Filter.eq('biome', 'CAATINGA'))
                        .filter(ee.Filter.eq('year', yearEst)); 
var imgColGeodatin = ee.ImageCollection(param.asset_geodatin_next).merge(
                          ee.ImageCollection(param.asset_geodatin_work))
                              .filter(ee.Filter.eq('year', yearEst)); 
print("show imagens geodatin ", imgColGeodatin)
print("lista de bandas Mapbiomas ", imgColMapBiomas.limit(10));
print("lista de bandas Geodatin ", imgColGeodatin.limit(10));

Map.addLayer(imgColMapBiomas.mosaic(), vis.mosaic_Median_Wet, "mapbiomas median Wet", false);
Map.addLayer(imgColMapBiomas.mosaic(), vis.mosaic_Median_Dry, "mapbiomas median dry", false);
Map.addLayer(geoGrade, {color: 'green'}, "tile", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaicQM_Wet, "geodatin QM Wet", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_QM_Dry, "geodatin QM Dry", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Wet, "geodatin Median Wet");
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Dry, "geodatin Median Dry");

// Paint all the polygon edges with the same number and width, display.
var outline = ee.Image().toByte().paint({
  featureCollection: geoGrade,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: 'FF0000'}, 'bordeGrade');