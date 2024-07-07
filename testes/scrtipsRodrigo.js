var geometry = 
    /* color: #d63000 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[-40.231804078876614, -13.173196663547518],
          [-40.231804078876614, -14.772325765997135],
          [-38.248771852314114, -14.772325765997135],
          [-38.248771852314114, -13.173196663547518]]], null, false);


function exportImg (clipIm, nnewNome, nbuffer){
    var param ={
        image: clipIm.toInt8(), //.toUint8(),
        description: nnewNome,
        folder: "MAPBIOMAS-EXPORT",
        maxPixels:1e13,
        scale: 30,
        region: nbuffer
    }
    Export.image.toDrive(param)
}

var anos = ee.List.sequence(1985,2022, 1).getInfo();
print("lista de anos ", anos);
var assetCol8 = "projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1";
var collectionx = ee.Image(assetCol8);
print(collectionx, 'collection')
var classMapbiomas = [1,3,4,5,6,49,10,11,12,32,29,50,13,14,15,18,19,39,20,40,62,41,36,46,47,35,48, 9, 21, 22, 23, 24, 30, 25, 26, 33, 31, 27];
var reclass = [1,1,1,1,1, 2, 1, 1, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  3,  4,  4,  5,  4,  4,  6,  6,  6, 27];
var imgRemap = ee.Image().byte();

anos.forEach(function(ano) {
    var bandAct = 'classification_'+ ano.toString()
    var mosaico1 = collectionx.select(bandAct).clip(geometry);
    var mosaico= ee.Image(mosaico1).remap(classMapbiomas, reclass);
    imgRemap = imgRemap.addBands(mosaico.rename(bandAct));    
})

var collection = collectionx;
var Palettes = require('users/mapbiomas/modules:Palettes.js');
var palette = Palettes.get('classification2');
palette[13] = 'd1b390';

var visMap = {
    min: 1,
    max:62,
    palette : palette
} 

var coordlist =[
        [-38.999000 , -13.708000],[-39.019000 , -13.705000],[-39.023000 , -13.721000],
        [-39.027000 , -13.736000],[-39.028000 , -13.708000],[-39.042000 , -13.695000],
        [-39.055000 , -13.690000],[-39.073000 , -13.678000],[-39.239014 , -13.845203],
        [-39.218008 , -13.813992],[-39.239014 , -13.845203],[-39.460000 , -13.906000],
        [-39.479914 , -13.947525],[-39.168961 , -13.782947],[-39.164000 , -13.777000],
        [-39.477344 , -13.937703],[-39.481250 , -13.943169],[-39.310000 , -13.785000],
        [-39.309000 , -13.789000],[-39.471000 , -13.936000],[-39.469000 , -13.931000],
        [-39.178253 , -13.793558],[-39.174878 , -13.785056],[-39.233000 , -13.854000],
        [-39.234000 , -13.861000],[-39.239908 , -13.854033],[-39.236197 , -13.868000],
        [-39.240056 , -13.859661],[-39.248781 , -13.862581],[-39.208036 , -13.826569],
        [-39.164000 , -13.777000],[-39.230000 , -13.850000],[-39.233756 , -13.832536],
        [-39.221000 , -13.826000],[-39.174878 , -13.785056]
    ];
var mListCoord = ee.List([]);
coordlist.forEach(
            function(coord){
                mListCoord = mListCoord.add(ee.Geometry.Point(coord));
            }
        )

var points = ee.List(mListCoord);
print(points, "points");
var tambuffer = ee.List.sequence(500, 2500, 500).getInfo();
anos.forEach(function (ano){
    var nomeBand = "classification_" + String(ano);
    var imgTemp = imgRemap.select(nomeBand);
    tambuffer.forEach(function (buftam){            
        points.evaluate(function (coorda){
            coorda.forEach(function (coordlist){
                print("coordenadas ", coordlist)
                var buffer = ee.Geometry(coordlist).buffer(buftam , 10);
                var clipa = imgTemp.clip(buffer)
                clipa = clipa.set('banda', nomeBand,'ano', ano, 'buffer_size', buftam);
                var newNome = "Itu_Igrap_Buffer_clip_"+ nomeBand + "_" + buftam.toString();
                // Map.centerObject(clipa, 16)
                Map.addLayer(clipa, visMap, newNome, false); // visualizar
                exportImg(clipa, newNome, buffer);
                // print("print buffer", clipa);          

            })
        })
    })
})

