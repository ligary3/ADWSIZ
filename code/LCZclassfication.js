var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") 
var yrd = ee.FeatureCollection('users/lastrye00/LCZ/YRD')
var visualization = {
//  bands: ['LCZ_Filter'],
  min: 1,
  max: 17,
  palette: [
    '8c0000','d10000','ff0000','bf4d00','ff6600',
    'ff9955','faee05','bcbcbc','ffccaa','555555',
    '006a00','00aa00','648525','b9db79','000000',
    'fbf7ae','6a6aff'
    ]
};
Map.centerObject(yrd)
Map.setOptions('SATELLITE')
var RectRegion=yrd.union().geometry().bounds()
function prepSrL8(image) {
  // Develop masks for unwanted pixels (fill, cloud, cloud shadow).
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var getFactorImg = function(factorNames) {
    var factorList = image.toDictionary().select(factorNames).values();
    return ee.Image.constant(factorList);
  };
  var scaleImg = getFactorImg([
    'REFLECTANCE_MULT_BAND_.|TEMPERATURE_MULT_BAND_ST_B10']);
  var offsetImg = getFactorImg([
    'REFLECTANCE_ADD_BAND_.|TEMPERATURE_ADD_BAND_ST_B10']);
  var scaled = image.select('SR_B.|ST_B10').multiply(scaleImg).add(offsetImg);

  // Replace original bands with scaled bands and apply masks.
  return image.addBands(scaled, null, true)
    .updateMask(qaMask).updateMask(saturationMask);
}
var bands = ['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5',
             'SR_B6', 'SR_B7'];
var image = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(RectRegion)
  .filterDate('2024-01-01', '2024-12-01')// 
  .map(prepSrL8)
  .select(bands)
  .mean()
var labels=ee.List.sequence(1,17,1)
var nsamples=300
var trainpoints=ee.FeatureCollection(labels.map(function(value){
      return ee.FeatureCollection.randomPoints(classes.filterMetadata("Label","equals",ee.Number(value)),nsamples, 111)
      .map(function(feacol){
        return ee.Feature(feacol.geometry()).set("Label",ee.Number(value)) 
      })})).flatten()
      
var testpoints=ee.FeatureCollection(labels.map(function(value){
      return ee.FeatureCollection.randomPoints(classes.filterMetadata("Label","equals",ee.Number(value)),nsamples, 222)
      .map(function(feacol){
        return ee.Feature(feacol.geometry()).set("Label",ee.Number(value)) 
      })})).flatten()//no include 0
      
    var label='Label'
    var trainsample=image.sampleRegions({
      collection:trainpoints, 
      scale:500,
      tileScale:1, 
      geometries:false,
      properties: [label],
    })
//}
  var testsample=image.sampleRegions({
      collection:testpoints, 
      scale:500,
      tileScale:1, 
      geometries:false,
      properties: [label],
    })
  var trained = ee.Classifier.smileRandomForest({
      numberOfTrees:200,
      seed:1
    }).train(trainsample, "Label",bands);
    // Get a confusion matrix representing resubstitution accuracy.
    var trainAccuracy = trained.confusionMatrix();
    print('Resubstitution error matrix: ', trainAccuracy);
    print('Training overall accuracy: ', trainAccuracy.accuracy());
    print('Training kappa: ', trainAccuracy.kappa());
    print("Trained ConsumersAccuracy",trainAccuracy.consumersAccuracy())
    print("Trained producersAccuracy",trainAccuracy.producersAccuracy())
    // Classify the validation data.
    var validated = testsample.classify(trained);
    var testAccuracy = validated.errorMatrix('Label', 'classification');
    print('Validation error matrix: ', testAccuracy);
    print('Validation overall accuracy: ', testAccuracy.accuracy());
    print('Validation kappa: ', testAccuracy.kappa());
    print("ConsumersAccuracy",testAccuracy.consumersAccuracy())
    print("producersAccuracy",testAccuracy.producersAccuracy())
    
  var classified = image.classify(trained)
 // print(classified)//只有一个band
 
  //Map.addLayer(classified,visualization,"classified")
  Map.addLayer(classified.clip(RectRegion),visualization,"classified2")
  
//   Export.image.toDrive({
// image:classified.toUint8(),
// description:'lcz_yrd2', 
// folder:"LCZ_PredReuslt", 
// crs:'EPSG:32651',
// region:RectRegion,//yrd.geometry().bounds(),
// scale:100,
// fileFormat:'GeoTIFF',
// maxPixels:16000000000
// });