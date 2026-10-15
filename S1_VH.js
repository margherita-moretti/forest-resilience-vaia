/****************************************************************************
 * FOREST RECOVERY AFTER VAIA - SELECTED WINDTHROW POLYGONS
 * Areas: Regione Veneto Vaia polygons inside the AOI with > 100 fallen trees
 * Sentinel-2 NDVI + Sentinel-1 VH (terrain-flattened), July-August composites
 * Years: 2017-2026, year-over-year differences
 * Forest mask: CORINE Land Cover 2018 (pre-Vaia)
 ****************************************************************************/


/* ---------------------------------------------------------------------------
 * 0. CONFIGURATION
 * ------------------------------------------------------------------------- */

// Study area.
var AOI = ee.Geometry.Rectangle([11.90, 46.26, 12.19, 46.46], null, false);
Map.centerObject(AOI, 11);

// Regione Veneto layer "Perimetrazione dei danni al patrimonio forestale
// arrecati dalla tempesta Vaia", uploaded as an asset.
var VAIA = ee.FeatureCollection('projects/lunar-spring-501022-c3/assets/c0605101_schiantivaia');
var NAME_FIELD = 'nomcom';         // municipality name of each polygon
var FALLEN_FIELD = 'schianti_c';   // number of fallen trees
var MIN_FALLEN = 100;

print('Fields of the Vaia layer:', VAIA.first().propertyNames());

// Polygons inside the AOI with schianti_c > MIN_FALLEN.
var vaiaInAoi = VAIA.filterBounds(AOI).map(function (f) {
  var v = f.get(FALLEN_FIELD);
  var num = ee.Algorithms.If(
    ee.String(ee.Algorithms.ObjectType(v)).equals('String'),
    ee.Number.parse(v), v);
  return f.set('fallen_num', num);
});
var DAMAGED = vaiaInAoi.filter(ee.Filter.gt('fallen_num', MIN_FALLEN));

// One feature per municipality (all its selected polygons merged).
var comuniNames = DAMAGED.aggregate_array(NAME_FIELD).distinct().sort();
var COMUNI = ee.FeatureCollection(comuniNames.map(function (name) {
  var polys = DAMAGED.filter(ee.Filter.eq(NAME_FIELD, name));
  return ee.Feature(polys.geometry(ee.ErrorMargin(1)),
    {comune: name, n_polygons: polys.size()});
}));

// Processing is limited to the bounding box of the selected polygons.
var STUDY_BOX = DAMAGED.geometry(ee.ErrorMargin(1)).bounds(ee.ErrorMargin(1));
print('Vaia polygons in the AOI:', vaiaInAoi.size());
print('Selected polygons (' + FALLEN_FIELD + ' > ' + MIN_FALLEN + '):', DAMAGED.size());
print('Selected polygons by municipality:',
  DAMAGED.select([NAME_FIELD, FALLEN_FIELD]).sort(FALLEN_FIELD, false));
print('Municipalities (one line each in the charts):',
  COMUNI.aggregate_array('comune'));


var START_YEAR = 2017;
var END_YEAR = 2026;

var SUMMER_MONTHS = [7, 8];
var S1_START_MD = '07-01';
var S1_STOP_MD = '09-01';   
var S1_ORBIT = 'DESCENDING';

var SCALE = 20;                 // m
var EXPORT_CRS = 'EPSG:32632';  // UTM 32N, metric pixels
var EXPORT_FOLDER = 'GEE_ForestResilience';


/* ---------------------------------------------------------------------------
 * 1. MASK: pre-storm forest (CORINE 2018) inside the selected polygons
 * ------------------------------------------------------------------------- */

var FOREST_CLASSES = [311, 312, 313];
var forestMask = ee.Image('COPERNICUS/CORINE/V20/100m/2018')
  .select('landcover')
  .remap(FOREST_CLASSES, ee.List.repeat(1, FOREST_CLASSES.length), 0)
  .eq(1);

var studyMask = forestMask
  .and(ee.Image(0).paint(DAMAGED, 1))
  .clip(AOI).selfMask();

function applyMask(image) {
  return image.updateMask(studyMask).clip(AOI);
}

var areaHa = ee.Number(ee.Image.pixelArea().updateMask(studyMask).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: STUDY_BOX, scale: SCALE,
  crs: EXPORT_CRS, maxPixels: 1e10
}).get('area')).divide(1e4);
print('Analysed forest area (ha):', areaHa);

function safeReduce(col, bandName, reduceFn) {
  var fallback = ee.Image.constant(0).rename(bandName).selfMask();
  return ee.Image(ee.Algorithms.If(col.size().gt(0), reduceFn(col), fallback))
    .toFloat();
}


/* ---------------------------------------------------------------------------
 * 2. SENTINEL-2 NDVI
 * ------------------------------------------------------------------------- */

// SCL classes kept: 4 vegetation, 5 bare soil.
// Excluded: 7 unclassified, 11 snow, 6 water.
function maskS2(image) {
  var scl = image.select('SCL');
  var valid = scl.eq(4).or(scl.eq(5));
  return image.updateMask(valid).divide(10000)
    .copyProperties(image, ['system:time_start']);
}

var s2Ndvi = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(STUDY_BOX)
  .filterDate(ee.Date.fromYMD(START_YEAR, 1, 1), ee.Date.fromYMD(END_YEAR + 1, 1, 1))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 70))
  .select(['B4', 'B8', 'SCL'])
  .map(maskS2)
  .map(function (img) {
    return img.normalizedDifference(['B8', 'B4']).rename('NDVI')
      .copyProperties(img, ['system:time_start']);
  });

// Summer NDVI = mean of the July and August monthly medians.
function ndviSummer(year) {
  var monthly = SUMMER_MONTHS.map(function (m) {
    var start = ee.Date.fromYMD(year, m, 1);
    return safeReduce(s2Ndvi.filterDate(start, start.advance(1, 'month')), 'NDVI',
      function (c) { return c.median(); });
  });
  return applyMask(ee.ImageCollection.fromImages(monthly).mean()
    .rename('NDVI_' + year));
}


/* ---------------------------------------------------------------------------
 * 3. SENTINEL-1 VH (border noise correction + terrain flattening)
 * ------------------------------------------------------------------------- */

var wrapper = require('users/adugnagirma/gee_s1_ard:wrapper');
var DEM = ee.Image('USGS/SRTMGL1_003');

function s1Collection(year) {
  return wrapper.s1_preproc({
    START_DATE: year + '-' + S1_START_MD,
    STOP_DATE: year + '-' + S1_STOP_MD,
    POLARIZATION: 'VVVH',
    ORBIT: S1_ORBIT,
    GEOMETRY: STUDY_BOX,
    APPLY_BORDER_NOISE_CORRECTION: true,
    APPLY_TERRAIN_FLATTENING: true,
    TERRAIN_FLATTENING_MODEL: 'VOLUME',
    DEM: DEM,
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER: 0,
    APPLY_SPECKLE_FILTERING: false,
    FORMAT: 'DB',
    CLIP_TO_ROI: true,
    SAVE_ASSETS: false
  })[1].select('VH');
}

// Mean computed in linear power, then converted back to dB.
function vhSummer(col, year) {
  var lin = col.map(function (img) {
    return img.multiply(Math.LN10 / 10).exp().rename('VH')
      .copyProperties(img, ['system:time_start']);
  });
  var meanLin = safeReduce(lin, 'VH', function (c) { return c.mean(); });
  return applyMask(meanLin.log10().multiply(10).rename('VH_' + year));
}


/* ---------------------------------------------------------------------------
 * 4. ANNUAL STACKS AND YEAR-OVER-YEAR DIFFERENCES
 * ------------------------------------------------------------------------- */

var years = [];
var ndviList = [];
var vhList = [];
var s1Counts = [];

for (var y = START_YEAR; y <= END_YEAR; y++) {
  var s1 = s1Collection(y);
  years.push(y);
  ndviList.push(ndviSummer(y));
  vhList.push(vhSummer(s1, y));
  s1Counts.push(ee.Feature(null, {year: y, n_S1_scenes: s1.size()}));
}

var ndviStack = ee.Image.cat(ndviList);
var vhStack = ee.Image.cat(vhList);

var dNdviStack = ee.Image.cat(ndviList.slice(1).map(function (img, i) {
  return img.subtract(ndviList[i]).rename('dNDVI_' + years[i + 1]);
}));
var dVhStack = ee.Image.cat(vhList.slice(1).map(function (img, i) {
  return img.subtract(vhList[i]).rename('dVH_' + years[i + 1]);
}));

// Number of S1 scenes per summer.
print('Sentinel-1 scenes per summer:', ee.FeatureCollection(s1Counts));


/* ---------------------------------------------------------------------------
 * 5. MEANS OVER THE SELECTED POLYGONS
 * ------------------------------------------------------------------------- */

var means = ee.Dictionary(ndviStack.addBands(vhStack).reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: STUDY_BOX,
  scale: SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e10,
  tileScale: 4
}));

function safeDiff(a, b) {
  return ee.Algorithms.If(
    ee.Algorithms.IsEqual(a, null), null,
    ee.Algorithms.If(ee.Algorithms.IsEqual(b, null), null,
      ee.Number(a).subtract(ee.Number(b))));
}

var levelFC = ee.FeatureCollection(years.map(function (yr) {
  return ee.Feature(null, {
    year: yr,
    NDVI_mean: means.get('NDVI_' + yr),
    VH_mean_dB: means.get('VH_' + yr)
  });
}));

var diffFC = ee.FeatureCollection(years.slice(1).map(function (yr) {
  return ee.Feature(null, {
    year: yr,
    dNDVI_mean: safeDiff(means.get('NDVI_' + yr), means.get('NDVI_' + (yr - 1))),
    dVH_mean_dB: safeDiff(means.get('VH_' + yr), means.get('VH_' + (yr - 1)))
  });
}));

print('Annual means:', levelFC);
print('Year-over-year differences:', diffFC);


/* ---------------------------------------------------------------------------
 * 6. CHARTS - one time series per municipality (nomcom)
 * ------------------------------------------------------------------------- */

// Mean NDVI and VH of each municipality, for every year.
var perComune = ndviStack.addBands(vhStack).reduceRegions({
  collection: COMUNI,
  reducer: ee.Reducer.mean(),
  scale: SCALE,
  crs: EXPORT_CRS,
  tileScale: 4
});

// Long table: one row per municipality and year (year = numeric x axis).
var longTable = perComune.map(function (f) {
  return ee.FeatureCollection(years.map(function (yr) {
    return ee.Feature(null, {
      comune: f.get('comune'),
      year: yr,
      NDVI: f.get('NDVI_' + yr),
      VH: f.get('VH_' + yr)
    });
  }));
}).flatten();

function timeSeriesChart(yProp, title, vTitle) {
  return ui.Chart.feature.groups({
    features: longTable,
    xProperty: 'year',
    yProperty: yProp,
    seriesProperty: 'comune'
  })
  .setChartType('LineChart')
  .setOptions({
    title: title,
    vAxis: {title: vTitle},
    hAxis: {title: 'Year', format: '####', gridlines: {count: years.length}},
    pointSize: 4, lineWidth: 2, interpolateNulls: false,
    legend: {position: 'right'},
    chartArea: {width: '65%', height: '70%'},
    width: 900, height: 450
  });
}

print(timeSeriesChart('NDVI',
  'Summer NDVI (Jul-Aug) time series by municipality - Vaia areas (' + FALLEN_FIELD + ' > ' + MIN_FALLEN + ')',
  'NDVI'));
print(timeSeriesChart('VH',
  'Summer Sentinel-1 VH (Jul-Aug) time series by municipality - Vaia areas (' + FALLEN_FIELD + ' > ' + MIN_FALLEN + ')',
  'VH'));


/* ---------------------------------------------------------------------------
 * 7. MAP
 * ------------------------------------------------------------------------- */

var vis = {
  ndvi: {min: 0, max: 1, palette: ['d73027', 'fee08b', '1a9850']},
  vh: {min: -22, max: -10, palette: ['000000', 'ffffff']},
  dNdvi: {min: -0.4, max: 0.4, palette: ['d73027', 'ffffff', '1a9850']},
  dVh: {min: -4, max: 4, palette: ['2c7bb6', 'ffffff', 'd7191c']}
};

Map.addLayer(ndviStack.select('NDVI_2018'), vis.ndvi, 'NDVI 2018 (pre-storm)', false);
Map.addLayer(ndviStack.select('NDVI_' + END_YEAR), vis.ndvi, 'NDVI ' + END_YEAR, false);
Map.addLayer(vhStack.select('VH_2018'), vis.vh, 'VH 2018 (pre-storm)', false);
Map.addLayer(vhStack.select('VH_' + END_YEAR), vis.vh, 'VH ' + END_YEAR, false);
Map.addLayer(dNdviStack.select('dNDVI_2019'), vis.dNdvi, 'dNDVI 2019-2018');
Map.addLayer(dVhStack.select('dVH_2019'), vis.dVh, 'dVH 2019-2018 (dB)');
Map.addLayer(ee.Image().paint(DAMAGED, 0, 2), {palette: ['d7191c']}, 'Selected polygons');
Map.addLayer(ee.Image().paint(ee.FeatureCollection([ee.Feature(AOI)]), 0, 2),
  {palette: ['ffff00']}, 'AOI');


/* ---------------------------------------------------------------------------
 * 8. EXPORTS (Google Drive)
 * ------------------------------------------------------------------------- */

function exportImage(image, name) {
  Export.image.toDrive({
    image: image.toFloat(), description: name, folder: EXPORT_FOLDER,
    fileNamePrefix: name, region: AOI, scale: SCALE, crs: EXPORT_CRS,
    maxPixels: 1e13
  });
}


// Tables behind the charts (CSV).
Export.table.toDrive({
  collection: longTable,
  description: 'NDVI_VH_by_municipality_and_year',
  folder: EXPORT_FOLDER, fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: levelFC, description: 'Annual_means_NDVI_VH',
  folder: EXPORT_FOLDER, fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: diffFC, description: 'Year_over_year_differences_NDVI_VH',
  folder: EXPORT_FOLDER, fileFormat: 'CSV'
});