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
    mosaic_Min_Wet: {
        min:0, 
        max: 1500,
        bands: ["red_min_wet","green_min_wet","blue_min_wet"]
    },
    mosaic_Min_Dry: {
        min:0, 
        max: 1500,
        bands: ["red_min_dry","green_min_dry","blue_min_dry"]
    },
    mosaic_Max_Wet: {
        min:0, 
        max: 2500,
        bands: ["red_max_wet","green_max_wet","blue_max_wet"]
    },
    mosaic_Max_Dry: {
        min:0, 
        max: 2500,
        bands: ["red_max_dry","green_max_dry","blue_max_dry"]
    },
    mosaic_stdDev_Wet: {
        min:0, 
        max: 500,
        bands: ["red_stdDev_wet","green_stdDev_wet","blue_stdDev_wet"]
    },
    mosaic_stdDev_Dry: {
        min:0, 
        max: 500,
        bands: ["red_stdDev_dry","green_stdDev_dry","blue_stdDev_dry"]
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
     mapbiomas_Min: {
        min:0, 
        max: 1500,
        bands: ["red_min","green_min","blue_min"]
    },
    // mapbiomas_Max: {
    //     min:0, 
    //     max: 10000,
    //     bands: ["red_max","green_max","blue_max"]
    // },
    mapbiomas_stdDev: {
        min:0, 
        max: 750,
        bands: ["red_stdDev","green_stdDev","blue_stdDev"]
    },
};

var param = {
    'asset_mapbiomas': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-3',
    'asset_geodatin': 'projects/nexgenmap/MapBiomas2/SENTINEL/mosaics-CAATINGA-4',
    'assetCaat': 'users/CartasSol/shapes/nCaatingaBff3000',
    'gradeS2Corr': 'projects/mapbiomas-arida/ALERTAS/auxiliar/shpGradeSent_Caat_42reg',
    'asset_S2harmonic': 'COPERNICUS/S2_HARMONIZED',  // 2015-06-23T00:00:00
    'asset_S2Mask': "COPERNICUS/S2_CLOUD_PROBABILITY"
};
var yearEst = 2016;
var imgColMapBiomas = ee.ImageCollection(param.asset_mapbiomas)
                        .filter(ee.Filter.eq('biome', 'CAATINGA'))
                        .filter(ee.Filter.eq('year', yearEst)); 
var imgColGeodatin = ee.ImageCollection(param.asset_geodatin);

print("lista de bandas Mapbiomas ", imgColMapBiomas.limit(10));
print("lista de bandas Geodatin ", imgColGeodatin.limit(10));

Map.addLayer(imgColMapBiomas.mosaic(), vis.mosaic_Median_Wet, "mapbiomas median Wet");
Map.addLayer(imgColMapBiomas.mosaic(), vis.mosaic_Median_Dry, "mapbiomas median dry");
// Map.addLayer(imgColMapBiomas.mosaic(), vis.mapbiomas_Max, "mapbiomas Max ");
Map.addLayer(imgColMapBiomas.mosaic(), vis.mapbiomas_Min, "mapbiomas Min ");
Map.addLayer(imgColMapBiomas.mosaic(), vis.mapbiomas_stdDev, "mapbiomas stdDev");


Map.addLayer(imgColGeodatin.mosaic(), vis.mosaicQM_Wet, "geodatin QM Wet", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_QM_Dry, "geodatin QM Dry", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Min_Wet, "geodatin Min Wet", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Min_Dry, "geodatin Min Dry", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Max_Wet, "geodatin Max Wet", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Max_Dry, "geodatin Max Dry", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_stdDev_Wet, "geodatin stdDev Wet", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_stdDev_Dry, "geodatin stdDev Dry", false);
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Wet, "geodatin Median Wet");
Map.addLayer(imgColGeodatin.mosaic(), vis.mosaic_Median_Dry, "geodatin Median Dry");