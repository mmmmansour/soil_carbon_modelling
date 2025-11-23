# soil_carbon_modelling

// This script models forest aboveground biomass (AGB) using the Global Ecosystem Dynamics Investigation (GEDI) L4B.
// The GEDI L4B product provides 1 km x 1 km estimates of mean aboveground biomass density (AGBD). 
// More information about GEDI L4B is available at: https://daac.ornl.gov/GEDI/guides/GEDI_L4B_Gridded_Biomass.html.
// Sentinel-1, Sentinel-2, SRTM elevation and slope data are predictor variables.
// We will derive the forest mask from the ESA Global Land Cover dataset (2020).
//var carbone = carbone.multiply(0.5).rename('carbon_sample')
//var carbone = carbone.reproject({crs: 'EPSG:32628', scale: 10});
// Load Sentinel-1 for the post-rainy season.
var S1_PRS = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterDate('2021-09-01', '2022-10-30')
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
    .filterBounds(boundary);

// Prepare inter-quartile range (IQR) 
var S1_PRS_pc = S1_PRS.reduce(ee.Reducer.percentile([25,50,75]));

// Convert to natural units (linear units, which can be averaged)
var S1_PRS_pc = ee.Image(10).pow(S1_PRS_pc.divide(10));

var S1_PRS_pc_Feats = S1_PRS_pc.select(['VH_p50','VV_p50']).clip(boundary);

// Reproject to WGS 84 UTM zone 35s                
var S1_PRS_pc_Feats = S1_PRS_pc_Feats.reproject({crs: 'EPSG:32628',scale: 10}); 
  
// Check projection information
print('Projection, crs, and crs_transform:', S1_PRS_pc_Feats.projection());    

// Calculate inter-quartile range (IQR), a measure of Sentinel-1 backscatter variability
var PRS_VV_iqr = S1_PRS_pc_Feats.addBands((S1_PRS_pc.select('VV_p75').subtract(S1_PRS_pc.select('VV_p25'))).rename('VV_iqr'));
var PRS_VH_iqr = S1_PRS_pc_Feats.addBands((S1_PRS_pc.select('VH_p75').subtract(S1_PRS_pc.select('VH_p25'))).rename('VH_iqr'));

// Print the image to the console
print('Post-rainy Season VV IQR', PRS_VV_iqr);
// Print the image to the console
print('Post-rainy Season VV IQR', PRS_VH_iqr);

// Display S1 inter-quartile range imagery
Map.addLayer(PRS_VV_iqr.clip(boundary), {'bands': 'VV_iqr', min: 0,max: 0.1}, 'Sentinel-1 IW VV');
Map.addLayer( PRS_VH_iqr.clip(boundary), {'bands': 'VH_iqr', min: 0,max: 0.1}, 'Sentinel-1 IW VH');

/////////////////////
// Load Sentinel-2 spectral reflectance data.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR');

// Create a function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = ee.Number(2).pow(11).int();

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Return the masked and scaled data.
  return image.updateMask(mask).divide(10000);
}

// Filter clouds from Sentinel-2 for a given period.
var composite = s2.filterDate('2022-09-01', '2022-10-30')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .select('B2', 'B3','B8A', 'B4','B5','B6','B7','B8','B11', 'B12'); 

// Reproject to WGS 84 UTM zone 35s                  
var S2_composite = composite.median().reproject({crs: 'EPSG:32628', scale: 10});

// Check projection information                 
print('Projection, crs, and crs_transform:', S2_composite.projection());

// Display a composite S2 imagery
Map.addLayer(S2_composite.clip(boundary), {bands: ['B11', 'B8', 'B3'], min: 0, max: 0.3});
////////////////

  var BAI= S2_composite.expression(
                  '1/( (0.1-RED)*(0.1-RED)+(0.06-NIR)*(0.06-NIR) )',{
                   
                  'RED': S2_composite.clip(boundary).select('B4'),
                  'NIR': S2_composite.clip(boundary).select('B8')
                  }).rename('BAI'); 
                
Map.addLayer(BAI, {}, 'BAI')
var NDVI = S2_composite.clip(boundary).normalizedDifference(['B8','B4']).rename('NDVI')
var SAVI = S2_composite.expression(
                  '((NIR - RED) / (NIR + RED+0.5))* 1.5 ', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4')
                  }).rename('SAVI');
                  
var EVI =  S2_composite.expression(
                  '2.5*(NIR - RED) / (NIR + 6*RED-7.5*BLUE+1) ', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4'),
                  'BLUE': S2_composite.clip(boundary).select('B2')
                  }).rename('EVI');
                  

var OSAVI =S2_composite.expression(
                  '(NIR - RED) / (NIR + RED+0.16)', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4')
                  }).rename('OSAVI');
                  
                  
var MSAVI2 = S2_composite.expression(
                  '(0.5)*(2*(NIR +1)-sqrt((2*NIR+1)*(2*NIR+1)-8 *(NIR-RED)))', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4')
                  }).rename('MSAVI2');
                  
                  
var ARVI = S2_composite.expression(
                  '(NIR-(2*RED-BLUE))/(NIR+(2*RED-BLUE))', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4'),
                  'BLUE': S2_composite.clip(boundary).select('B2')
                  }).rename('ARVI');
                 
                 
                 
var IRECI = S2_composite.expression(
                  '(NIR - RED) / ( RED1 - RED2)', {
                  'NIR': S2_composite.clip(boundary).select('B8'),
                  'RED': S2_composite.clip(boundary).select('B4'),
                  'RED1': S2_composite.clip(boundary).select('B5'),
                  'RED2': S2_composite.clip(boundary).select('B6')
                  }).rename('IRECI');
                  
                  
                  
var MRENDVI = S2_composite.expression(
                  '(RED2 - RED1)/(RED2 + RED1 - 2 * BLUE)', {
                  'BLUE': S2_composite.clip(boundary).select('B2'),
                  'RED1': S2_composite.clip(boundary).select('B5'),
                  'RED2': S2_composite.clip(boundary).select('B6')
                  }).rename('MRENDVI');
                  
                  
var RENDVI =  S2_composite.expression(
                  '(RED2 - RED1)/(RED2 + RED1)', {
                  'RED1': S2_composite.clip(boundary).select('B5'),
                  'RED2': S2_composite.clip(boundary).select('B6')
                  }).rename('RENDVI');
                  
                
                
var MRESR  = S2_composite.expression(
                  '(RED2 - BLUE)/(RED1-BLUE)', {
                  'BLUE': S2_composite.clip(boundary).select('B2'),
                  'RED1': S2_composite.clip(boundary).select('B5'),
                  'RED2': S2_composite.clip(boundary).select('B6')
                  }).rename('MRESR');
                                    
//////////////////
// Load SRTM
var SRTM = ee.Image("USGS/SRTMGL1_003");
// Clip Elevation
var elevation = SRTM.clip(boundary);

// Reproject 'elevation' to WGS 84 UTM zone 35s                
var elevation = elevation.reproject({crs: 'EPSG:32628',scale: 10}); 
  
// Check projection information
print('Projection, crs, and crs_transform:', elevation.projection()); 

// Derive slope from the SRTM
var slope = ee.Terrain.slope(SRTM).clip(boundary);

// Reproject 'slope' to WGS 84 UTM zone 35s                
var slope = slope.reproject({crs: 'EPSG:32628',scale: 30}); 
  
// Check projection information
print('Projection, crs, and crs_transform:', slope.projection()); 

/////////////////
var dataset = ee.ImageCollection("ESA/WorldCover/v100").first();

// Clip the land cover to the boundary
var ESA_LC_2020 = dataset.clip(boundary);

// Extract forest areas from the land cover
var forest_mask = ESA_LC_2020.updateMask(
  ESA_LC_2020.eq(10) // Only keep pixels where class equals 2
);

// Display forests only
var visualization = {bands: ['Map'],};

Map.addLayer(forest_mask, visualization, "Trees");
//////////////////
var dataset = ee.ImageCollection('FIRMS').filterDate('2021-01-01', '2021-05-30').filterBounds(boundary).mean();
var dataset = dataset.select('T21').rename('T21');
var mean_fire = dataset.reduce(ee.Reducer.mean())
var fire_roi = mean_fire.clip(boundary)

var firesVis = {
  min: 325.0,
  max: 400.0,
  palette: ['red', 'orange', 'yellow'],
};
//Map.setCenter(-119.086, 47.295, 6);
Map.addLayer(fire_roi, firesVis, 'Fires');

Export.image.toDrive({
  image: fire_roi,
  description: 'fire_roi',
  scale: 20,
  crs: 'EPSG:32628', // EPSG:32735 (WGS 84 UTM Zone 35S)
  maxPixels: 6756353855,
  region: boundary
});
///////////////////
// Merge the predictor variables
var mergedCollection = S2_composite.addBands(NDVI.addBands(MRESR.addBands(IRECI.addBands(MRENDVI.addBands(RENDVI.addBands(SAVI.addBands(OSAVI.addBands(EVI.addBands(MSAVI2.addBands(ARVI.addBands(PRS_VV_iqr.addBands(dataset.addBands(BAI.addBands(PRS_VH_iqr.addBands(elevation.addBands(slope.addBands(forest_mask)))))))))))))))));

// Clip to the output image to Harare study area boundary.
var clippedmergedCollection = mergedCollection.clipToCollection(boundary);
print('clippedmergedCollection: ', clippedmergedCollection);
//Map.addLayer(clippedmergedCollection, {bands: ['B8', 'B4', 'B3'], max: 0.3}, 'mergedCollection');

// Bands to include in the classification
var bands = ['B2','B3','B6', 'B8A','B11', "B5",'B12','B4','elevation', 'B8', 'NDVI','OSAVI', 'MSAVI2', 'EVI', 'ARVI','MRESR','RENDVI','VH_iqr'];
//var bands = ['B11', 'B12', 'B5', 'B8', 'B2','NDVI']
////////////////////
// Prepare training dataset
// More information at https://developers.google.com/earth-engine/datasets/catalog/LARSE_GEDI_GEDI04_B_002

var l4b = ee.Image('LARSE/GEDI/GEDI04_B_002');

var dataset = l4b.select('MU').clip(boundary);
Map.setCenter(28.8713,-18.4492, 12);

// Reproject to WGS 84 UTM zone 35s                  
var dataset = dataset.reproject({crs: 'EPSG:32628', scale: 100});

// Check projection information                 
print('Projection, crs, and crs_transform:', dataset.projection());

// Display the GEDI L4B dataset
Map.addLayer(dataset,
    {min: 10, max: 250, palette: '440154,414387,2a788e,23a884,7ad151,fde725'},
    'Mean Biomass');

// Sample the training points from the dataset
var points = carbone.sample({
   region: boundary,
   scale: 0.5,
   numPixels: 1000000, 
   geometries: true});

// Print and display the points derived from the GEDI L4B dataset
print(points.size());
print(points.limit(10));

Map.addLayer(points);

// Split training data into training and testing sets 
// Add a random column (named random) and specify the seed value for repeatability
var datawithColumn = points.randomColumn('random', 27);

// Use 70% for training, 30% for validation
var split = 0.8; 
var trainingData = datawithColumn.filter(ee.Filter.lt('random', split));
print('training data', trainingData);

var validationData = datawithColumn.filter(ee.Filter.gte('random', split));
print('validation data', validationData);

////////////////////
// Perform random forest regression

// Collect training data
var training = clippedmergedCollection.select(bands).sampleRegions({
  collection: trainingData,
  properties: ['b1'],
  scale: 10 // Need to change the scale of training data to avoid the 'out of memory' problem
  });

// Train a random forest classifier for regression 
var classifier = ee.Classifier.smileRandomForest(50)
  .setOutputMode('REGRESSION')
  .train({
    features: training, 
    classProperty: "b1",
    inputProperties: bands
    });

//Run the classification and clip it to the boundary
var regression = clippedmergedCollection.select(bands).classify(classifier, 'predicted').clip(boundary);

// Load and define a continuous palette
var palettes = require('users/gena/packages:palettes');

// Choose and define a palette
var palette = palettes.kovesi.rainbow_bgyr_35_85_c73[7].reverse();

// Display the input imagery and the regression classification.
  // get dictionaries of min & max predicted value
  var regressionMin = (regression.reduceRegion({
    reducer: ee.Reducer.min(),
    scale: 10, 
    crs: 'EPSG:32628',
    geometry: boundary,
    bestEffort: true,
    tileScale: 5
  }));
  
  var regressionMax = (regression.reduceRegion({
    reducer: ee.Reducer.max(),
    scale: 10, 
    crs: 'EPSG:32628',
    geometry: boundary,
    bestEffort: true,
    tileScale: 5
  }));
  
// Add to map
var viz = {palette: palette, min: regressionMin.getNumber('predicted').getInfo(), max: regressionMax.getNumber('predicted').getInfo()};
Map.addLayer(regression, viz, 'Regression');

// Create the panel for the legend items.
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Create and add the legend title.
var legendTitle = ui.Label({
  value: 'g CO/KG_soil',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});

legend.add(legendTitle);

// create the legend image
var lon = ee.Image.pixelLonLat().select('latitude');
var gradient = lon.multiply((viz.max-viz.min)/100.0).add(viz.min);
var legendImage = gradient.visualize(viz);
 
// create text on top of legend
var panel = ui.Panel({
widgets: [
ui.Label(viz['max'])
],
});
 
legend.add(panel);
 
// create thumbnail from the image
var thumbnail = ui.Thumbnail({
image: legendImage,
params: {bbox:'0,0,10,100', dimensions:'10x200'},
style: {padding: '1px', position: 'bottom-center'}
});
 
// add the thumbnail to the legend
legend.add(thumbnail);
 
// create text on top of legend
var panel = ui.Panel({
widgets: [
ui.Label(viz['min'])
],
});

legend.add(panel);
Map.add(legend);

// Zoom to the regression on the map
Map.centerObject(boundary, 11);

//////////////////////////
// Check model performance
// Get details of classifier
var classifier_details = classifier.explain();

// Explain the classifier with importance values
var importance = ee.Dictionary(classifier_details.get('importance'))
var keys = importance.keys().sort(importance.values()).reverse()
var values = importance.values(keys);
var rows = keys.zip(values).map(function(list) {
  return {c: ee.List(list).map(function(n) { return {v: n}; })}
})

var dataTable = {
  cols: [{id: 'band', label: 'Band', type: 'string'},
         {id: 'importance', label: 'Importance', type: 'number'}],
  rows: rows
};

ee.Dictionary(dataTable).evaluate(function(result) {
  var chart = ui.Chart(result)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Random Forest Band Importance',
      legend: {position: 'none'},
      hAxis: {title: 'Bands'},
      vAxis: {title: 'Importance',
      gridlines: {color: 'transparent'}
      }
      ,
      series: {
      0: {
      
        dataOpacity: 1.5,
        gridlines: {color: 'transparent'},
        color: 'black'
       
      },
    
    },
    });
  print(chart);
})


// Create model assessment statistics
// Get predicted regression points in same location as training data
var predictedTraining = regression.sampleRegions({collection:trainingData, geometries: true});

// Separate the observed (agbd_GEDI) and predicted (regression) properties
var sampleTraining = predictedTraining.select(['b1', 'predicted']);

// Create chart, print it
var chartTraining = ui.Chart.feature.byFeature(sampleTraining, 'b1', 'predicted')
.setChartType('ScatterChart').setOptions({
title: 'Predicted vs Observed - Training data ',
hAxis: {'title': 'observed'},
vAxis: {'title': 'predicted'},
pointSize: 3,
series: {
      0: {
      
        dataOpacity: 1.5,
        type: 'linear',
        color: '3c7a18',
        visibleInLegend: true,
        showR2: true
       
      },
    
    },
trendlines: {  0: {
        opacity: 1,
        type: 'linear',
        color: 'red',
        visibleInLegend: true,
        showR2: true
      } ,
1: {
        opacity: 4,
        type: 'linear',
        color: 'black',
        visibleInLegend: true
      }}});
print(chartTraining);

// Compute Root Mean Squared Error (RMSE)
// Get array of observation and prediction values 
var observationTraining = ee.Array(sampleTraining.aggregate_array('b1'));

var predictionTraining = ee.Array(sampleTraining.aggregate_array('predicted'));

// Compute residuals
var residualsTraining = observationTraining.subtract(predictionTraining);

// Compute RMSE with equation and print the result
var rmseTraining = residualsTraining.pow(2).reduce('mean', [0]).sqrt();
print('Training RMSE', rmseTraining);

/////////////////////
//Perform validation
// Get predicted regression points in same location as validation data
var predictedValidation = regression.sampleRegions({collection:validationData, geometries: true});

// Separate the observed (MU) and predicted (regression) properties
var sampleValidation = predictedValidation.select(['b1', 'predicted']);

// Create chart and print it
var chartValidation = ui.Chart.feature.byFeature(sampleValidation, 'b1', 'predicted')
.setChartType('ScatterChart').setOptions({
title: 'Predicted vs Observed - Validation data',
hAxis: {'title': 'predicted'},
vAxis: {'title': 'observed'},
 series: {
      0: {
      
        dataOpacity: 1.5,
        type: 'linear',
        color: '3c7a18',
        visibleInLegend: true,
        showR2: true
       
      },
    
    },
pointSize: 3,
trendlines: {  0: {
        opacity: 1,
        type: 'linear',
        color: 'black',
        visibleInLegend: true,
        showR2: true
      } ,
1: {
        opacity: 4,
        type: 'linear',
        color: 'black',
        visibleInLegend: true}}});
print(chartValidation);

// Compute RMSE
// Get array of observation and prediction values 
var observationValidation = ee.Array(sampleValidation.aggregate_array('b1'));

var predictionValidation = ee.Array(sampleValidation.aggregate_array('predicted'));

// Compute residuals
var residualsValidation = observationValidation.subtract(predictionValidation);

// Compute RMSE with equation and print it
var rmseValidation = residualsValidation.pow(2).reduce('mean', [0]).sqrt();
print('Validation RMSE', rmseValidation);

//////////////////
// Export the image, specifying scale and region.
Export.image.toDrive({
  image: regression,
  description: 'Muf_AGBD_GEDI_2021',
  scale: 20,
  crs: 'EPSG:32628', // EPSG:32735 (WGS 84 UTM Zone 35S)
  maxPixels: 6756353855,
  region: boundary
});
