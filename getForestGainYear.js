var startYear=1982;
var endYear=2020;
var years=["1982","1983","1984","1985","1986","1987","1988","1989","1990","1991","1992","1993",
            "1994","1995","1996","1997","1998","1999","2000","2001","2002","2003","2004",
            "2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015",
            "2016","2017","2018","2019","2020"]
var forestNow=ee.Image(0).expression('gedi>0?1:(gedi==0?0:(ig==2?1:0))',{
  'gedi':gedi,
  'ig':image_global
})
var hansenGain=hansen.select('gain')
hansenGain=hansenGain.updateMask(hansenGain.gt(0));
var hansen2000=hansen.select('treecover2000')
var forestUnchanged=hansen2000 
var unforestUnchanged=ee.Image(1)

//----------functions-----------//
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
};
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};
//Function used to calculate NDVI & NBR
var getL8ND = function(img) {
  var b4=img.select('B4')
  b4=b4.multiply(0.9372).add(0.0123)
  var b5=img.select('B5')
  b5=b5.multiply(0.8339).add(0.0448)
  var b7=img.select('B7')
  b7=b7.multiply(0.9165).add(0.0116)
  var etm=b5.addBands(b7)
  var etm0=b5.addBands(b4)
  var ndvi = etm.normalizedDifference().rename('NDVI')
  var nbr = etm.normalizedDifference().rename('NBR')
  return ndvi.addBands(nbr).set('system:time_start', img.get('system:time_start'));
};
var getL457ND = function(img) {
  var ndvi = img.normalizedDifference(['B4', 'B3']).rename('NDVI')
  var nbr = img.normalizedDifference(['B4', 'B7']).rename('NBR');
  return ndvi.addBands(nbr).set('system:time_start', img.get('system:time_start'))
};

//----------load landsat---------//
for(var year=startYear;year<=2020;year++)
{
  var l4_col = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
                  .filterDate(year+'-06-01', year+'-09-01')
                  .map(cloudMaskL457).map(getL457ND);
  var l5_col = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                  .filterDate(year+'-06-01', year+'-09-01')
                  .map(cloudMaskL457).map(getL457ND);
  var l7_col = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                .filterDate(year+'-06-01', year+'-09-01')
                .map(cloudMaskL457).map(getL457ND);
  var l8_col = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                .filterDate(year+'-06-01', year+'-09-01')
                .map(maskL8sr).map(getL8ND);
  var lan_col=l4_col.merge(l5_col).merge(l7_col).merge(l8_col);
  var ndviMin=lan_col.select('NDVI').min();
  var ndviMax=lan_col.select('NDVI').max();
  var nbrYear=lan_col.select('NBR').max();
  var nbrOrg=nbrYear;
  nbrYear=nbrYear.expression('nbr>=-1&&nbr<=1?nbr:-2',{'nbr':nbrYear})
            .set('system:time_start', ee.Date(year+'-02-01').millis());
  forestUnchanged=forestUnchanged.expression('ndvi>=-1&&ndvi<=1?(ndvi>0.75&&f>10?f:0):f',{
    'ndvi':ndviMin,
    'f':forestUnchanged
  })
  unforestUnchanged=unforestUnchanged.expression('f>10?0:(ndvi<0.45&&uf==1?uf:0)',{
    'ndvi':ndviMax,
    'uf':unforestUnchanged,
    'f':hansen2000
  })
  if(year==startYear){
    var nbrYears=ee.ImageCollection([nbrYear]);
    var nbrOrgs=ee.ImageCollection([nbrOrg]);
  }
  else{
    nbrYears=nbrYears.merge(nbrYear);
    nbrOrgs=nbrOrgs.merge(nbrOrg);
  }
}

//----------gap filling
var increase = ee.Date('1992-01-01').millis().subtract(ee.Date('1991-01-01').millis());
var nbrMosaic=nbrOrgs.sort('system:time_start').mosaic()
var nbr2020=nbrYears.filter(ee.Filter.eq('system:time_start',ee.Date('2020-02-01').millis())).first();
nbr2020=nbr2020.expression('n==-2&&m>=-1&&m<=1?m:n',{
  'n':nbr2020,
  'm':nbrMosaic
})
var accumulate = function(image, list) {
  var time=ee.Number(image.get('system:time_start'));
  var yearM1=nbrYears.filterDate(ee.Date(time.subtract(increase)).advance(-1, 'month'),ee.Date(time.subtract(increase)).advance(1, 'month')).first();
  var yearA1=ee.Image(ee.List(list).get(-1));
  var imageCheck=image.expression('nbr!=-2?nbr:(ym1!=-2&&ya1!=-2?(ym1+ya1)/2:(ym1!=-2?ym1:ya1))',{
    'nbr':image,
    'ym1':yearM1,
    'ya1':yearA1
  }).set('system:time_start',time);
  return ee.List(list).add(imageCheck);
}
var first = ee.List([nbr2020]);
var checkCol=nbrYears.filterDate('1983-01-01','2019-12-31').sort('system:time_start',false);
var check = ee.ImageCollection(ee.List(checkCol.iterate(accumulate, first)));
var startYearImg=nbrYears.filter(ee.Filter.eq('system:time_start',ee.Date('1982-02-01').millis())).first();
var startYearCheckImg=startYearImg.expression('nbr!=-2?nbr:pre',{
  'nbr':startYearImg,
  'pre':check.filter(ee.Filter.eq('system:time_start',ee.Date('1983-02-01').millis())).first()
}).set('system:time_start',ee.Date('1982-02-01').millis());
check=check.merge(startYearCheckImg);
nbrYears=check.sort('system:time_start');
nbrYears=nbrYears.map(function(img){
  return img.updateMask(forestUnchanged.eq(0)
            .and(unforestUnchanged.eq(0)))
})
print(nbrYears,'nbrYears')
Map.addLayer(nbrYears,{},'nbrYears')

//----------quality ImageCollection
for(var year=1982;year<=2020;year++)
{
  var l4_col = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR')
                  .filterDate(year+'-06-01', year+'-09-01')
                  .map(cloudMaskL457).map(getL457ND);
  var l5_col = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
                  .filterDate(year+'-06-01', year+'-09-01')
                  .map(cloudMaskL457).map(getL457ND);
  var l7_col = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
                .filterDate(year+'-06-01', year+'-09-01')
                .map(cloudMaskL457).map(getL457ND);
  var l8_col = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
                .filterDate(year+'-06-01', year+'-09-01')
                .map(maskL8sr).map(getL8ND);
  var lan_col=l4_col.merge(l5_col).merge(l7_col).merge(l8_col)
  var qaUnForest=lan_col.select('NBR').count()
  qaUnForest=qaUnForest.expression('q>0?1:0',{
    'q':qaUnForest
  }).rename('qa').set('system:time_start', ee.Date(year+'-02-01').millis())
  if(year==1982){var qaUnForestCol=ee.ImageCollection([qaUnForest])}
  else{qaUnForestCol=qaUnForestCol.merge(qaUnForest)}
}
var qaYears=qaUnForestCol.map(function(img){
  var time=img.get('system:time_start')
  var y=ee.Date(time).get('year')
  var sum=qaUnForestCol.filterDate('1982-01-01',ee.Date(time).advance(2, 'month')).sum()
  return sum.set('system:time_start',time)
})
Map.addLayer(qaYears,{},'qaYears')

//-----------get basic forest gain year for GFC-----------//
var increase = ee.Date('2001-01-01').millis().subtract(ee.Date('2000-01-01').millis());
var accumulate = function(img, list) {
  var time=ee.Number(img.get('system:time_start'))
  var yearM1=nbrYears.filterDate(ee.Date(time.subtract(increase)).advance(-1, 'month'),ee.Date(time.add(increase)).advance(1, 'month')).first()
  var sub=img.subtract(yearM1).set('system:time_start',time)
  return ee.List(list).add(sub);
}
var first = ee.List([
  nbrYears.filter(ee.Filter.eq('system:time_start',ee.Date('2001-02-01').millis())).first().subtract(
    nbrYears.filter(ee.Filter.eq('system:time_start',ee.Date('2000-02-01').millis())).first())
    .set('system:time_start',ee.Date('2001-02-01').millis())
]);
var subInp=nbrYears.filterDate('2002-01-01','2012-12-31')
var subExp = ee.ImageCollection(ee.List(subInp.iterate(accumulate, first)))
var subMax=subExp.max()
var gainYear=ee.Image(0).updateMask(hansenGain.gt(0))
print(subExp,'subExp')
for(var year=2001;year<=2012;year++){
  var sub=subExp.filter(ee.Filter.eq('system:time_start',ee.Date(year+'-02-01').millis())).first();
  gainYear=gainYear.expression('sub>=max?y:g',{
    'sub':sub,
    'max':subMax.subtract(0.01),
    'y':ee.Image(year),
    'g':gainYear
  })
}
var hansenGainYear=gainYear.updateMask(hansenGain.gt(0))

//-----------LandTrendr------------//
var landTrendr=ee.Algorithms.TemporalSegmentation.LandTrendr({
                timeSeries:nbrYears,
                maxSegments:6,
                spikeThreshold:0.9,
                vertexCountOvershoot:3,
                preventOneYearRecovery:true,
                recoveryThreshold:0.25,
                pvalThreshold:0.05,
                bestModelProportion:0.75
              })
var lt=landTrendr.select('LandTrendr').clip(study_area)
var vertexMask = lt.arraySlice(0, 3, 4)
var vertices = lt.arrayMask(vertexMask)
var ltFits=lt.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([years])
// construct segment start and end point years and index values
var left = vertices.arraySlice(1, 0, -1);    // slice out the vertices as the start of segments
var right = vertices.arraySlice(1, 1, null); // slice out the vertices as the end of segments
var startYear = left.arraySlice(0, 0, 1);    // get year dimension of LT data from the segment start vertices
var startVal = left.arraySlice(0, 2, 3);     // get spectral index dimension of LT data from the segment start vertices
var endYear = right.arraySlice(0, 0, 1);     // get year dimension of LT data from the segment end vertices 
var endVal = right.arraySlice(0, 2, 3);      // get spectral index dimension of LT data from the segment end vertices
var dur = endYear.subtract(startYear);       // subtract the segment start year from the segment end year to calculate the duration of segments 
var mag = endVal.subtract(startVal);         // substract the segment start index value from the segment end index value to calculate the delta of segments
var distImg = ee.Image.cat([startYear,endYear, mag, dur, startVal,endVal]).toArray(0);
var distImgSorted = distImg.arraySort(mag.multiply(-1));
var vertexSum=vertexMask.arrayReduce(ee.Reducer.sum(), [1]).arrayProject([1]).arrayFlatten([['sum']])
distImgSorted=distImgSorted.updateMask(vertexSum.gt(2)).updateMask(forestNow.eq(1))
//--------------first
var tempDistImg = distImgSorted.arraySlice(1, 0, 1).unmask(ee.Image(ee.Array([[0],[0],[0],[0],[0],[0],[0]])));
var firstDistImg = ee.Image.cat(tempDistImg.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['preYear']]), 
                                tempDistImg.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['gainYear']]),
                                tempDistImg.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['mag']]),
                                tempDistImg.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['dur']]),
                                tempDistImg.arraySlice(0,4,5).arrayProject([1]).arrayFlatten([['preval']]),
                                tempDistImg.arraySlice(0,5,6).arrayProject([1]).arrayFlatten([['gainval']]);
firstDistImg=firstDistImg.updateMask(firstDistImg.select('mag').gt(0.1))
var gainThrd=ee.Image(0).clip(study_area)
var preThrd=ee.Image(0).clip(study_area)
var gainQa=ee.Image(0).clip(study_area)
var preQa=ee.Image(0).clip(study_area)
var prepreQa=ee.Image(0).clip(study_area)
var firstGainYear=ee.Image(0).clip(study_area)
for(var year=1982;year<=2020;year++){
  var gy=firstDistImg.select('gainYear')
  var preYear=firstDistImg.select('preYear')
  gainThrd=gainThrd.expression('gainYear==year?newThrd:thrd',{
    'gainYear':gy,
    'year':ee.Image(year).clip(study_area),
    'newThrd':ee.Image(thrd[year-1982]).clip(study_area),
    'thrd':gainThrd
  })
  preThrd=preThrd.expression('preYear==year?newThrd:thrd',{
    'preYear':preYear,
    'year':ee.Image(year).clip(study_area),
    'newThrd':ee.Image(thrd[year-1982]).clip(study_area),
    'thrd':preThrd
  })
  firstGainYear=firstGainYear.expression('fg==0&&y-py>1&&y<=gy&&yv>yt?y:fg',{
    'fg':firstGainYear,
    'y':ee.Image(year).clip(study_area),
    'py':preYear,
    'gy':gy,
    'yv':ltFits.select(year.toString()),
    'yt':ee.Image(thrd[year-1982]).clip(study_area)
  })
  gainQa=gainQa.expression('gainYear==year?newQa:qa',{
    'gainYear':firstGainYear,
    'year':ee.Image(year).clip(study_area),
    'newQa':qaYears.filterDate(year+'-01-01',year+'-12-31').first().clip(study_area),
    'qa':gainQa
  })
  preQa=preQa.expression('preYear==year?newQa:qa',{
    'preYear':preYear,
    'year':ee.Image(year).clip(study_area),
    'newQa':qaYears.filterDate(year+'-01-01',year+'-12-31').first().clip(study_area),
    'qa':preQa
  })
  if(year>1982){
    prepreQa=prepreQa.expression('preYear==year?newQa:qa',{
      'preYear':preYear,
      'year':ee.Image(year).clip(study_area),
      'newQa':qaYears.filterDate((year-1).toString()+'-01-01',(year-1).toString()+'-12-31').first().clip(study_area),
      'qa':prepreQa
    })
  }
}
var qa=gainQa.expression('pq>0&&pq-pp>0?(gq-pq)/(dur+1):0',{
  'pq':preQa,
  'pp':prepreQa,
  'gq':gainQa,
  'dur':firstGainYear.subtract(firstDistImg.select('preYear'))
})
var threshold = firstDistImg.select(['preval']).lte(preThrd)
                  .and(qa.gt(0.5))
var firstgainYear=firstGainYear.updateMask(threshold)
firstgainYear=ee.Image(0).expression('gy<=2020&&gy>1982?gy:0',{
  'gy':firstgainYear
}).clip(study_area)
var gainYear=firstgainYear

//----------secondToSixth
for (var i=2;i<=6;i++)
{
  var tempDistImg = distImgSorted.updateMask(vertexSum.gt(i)).arraySlice(1, i-1, i).unmask(ee.Image(ee.Array([[0],[0],[0],[0],[0],[0],[0]])));
  var distImg = ee.Image.cat(tempDistImg.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['preYear']]), 
                                tempDistImg.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['gainYear']]),
                                tempDistImg.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['mag']]),
                                tempDistImg.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['dur']]),
                                tempDistImg.arraySlice(0,4,5).arrayProject([1]).arrayFlatten([['preval']]),
                                tempDistImg.arraySlice(0,5,6).arrayProject([1]).arrayFlatten([['gainval']]);
  distImg=distImg.updateMask(distImg.select('mag').gt(0.1))
  var gainThrd=ee.Image(0).clip(study_area)
  var preThrd=ee.Image(0).clip(study_area)
  var gainQa=ee.Image(0).clip(study_area)
  var preQa=ee.Image(0).clip(study_area)
  var prepreQa=ee.Image(0).clip(study_area)
  var firstGainYear=ee.Image(0).clip(study_area)
  for(var year=1982;year<=2020;year++){
    var gy=distImg.select('gainYear')
    var preYear=distImg.select('preYear')
    gainThrd=gainThrd.expression('gainYear==year?newThrd:thrd',{
      'gainYear':gy,
      'year':ee.Image(year).clip(study_area),
      'newThrd':ee.Image(thrd[year-1982]).clip(study_area),
      'thrd':gainThrd
    })
    preThrd=preThrd.expression('preYear==year?newThrd:thrd',{
      'preYear':preYear,
      'year':ee.Image(year).clip(study_area),
      'newThrd':ee.Image(thrd[year-1982]).clip(study_area),
      'thrd':preThrd
    })
    firstGainYear=firstGainYear.expression('fg==0&&y-py>1&&y<=gy&&yv>yt?y:fg',{
      'fg':firstGainYear,
      'y':ee.Image(year).clip(study_area),
      'py':preYear,
      'gy':gy,
      'yv':ltFits.select(year.toString()),
      'yt':ee.Image(thrd[year-1982]).clip(study_area)
    })
    gainQa=gainQa.expression('gainYear==year?newQa:qa',{
      'gainYear':firstGainYear,
      'year':ee.Image(year).clip(study_area),
      'newQa':qaYears.filterDate(year+'-01-01',year+'-12-31').first().clip(study_area),
      'qa':gainQa
    })
    preQa=preQa.expression('preYear==year?newQa:qa',{
      'preYear':preYear,
      'year':ee.Image(year).clip(study_area),
      'newQa':qaYears.filterDate(year+'-01-01',year+'-12-31').first().clip(study_area),
      'qa':preQa
    })
    if(year>1982){
      prepreQa=prepreQa.expression('preYear==year?newQa:qa',{
        'preYear':preYear,
        'year':ee.Image(year).clip(study_area),
        'newQa':qaYears.filterDate((year-1).toString()+'-01-01',(year-1).toString()+'-12-31').first().clip(study_area),
        'qa':prepreQa
      })
    }
  }
  var qa=gainQa.expression('pq>0&&pq-pp>0?(gq-pq)/(dur+1):0',{
    'pq':preQa,
    'pp':prepreQa,
    'gq':gainQa,
    'dur':firstGainYear.subtract(distImg.select('preYear'))
  })
  var threshold = distImg.select(['preval']).lte(preThrd)
                    .and(qa.gt(0.5))
  var GainYear=firstGainYear.updateMask(threshold)
  GainYear=ee.Image(0).expression('gy<=2020&&gy>1982?gy:0',{
    'gy':GainYear
  }).clip(study_area)
  gainYear=gainYear.expression('g>0?g:(gy>0?gy:g)',{
      'g':gainYear,
      'gy':GainYear
    })
}

//----------mosaic----------//
var gainYearFinal=ee.Image(0).expression('hgy>0&&(gy==0||gy>2012||gy<2001)?hgy:gy',{
  'gy':gainYear,
  'hgy':hansenGainYear
}).clip(study_area)
gainYearFinal=gainYearFinal.reduceNeighborhood(ee.Reducer.mode(), ee.Kernel.square(30,'meters'),'kernel',false)
gainYear=gainYear.reduceNeighborhood(ee.Reducer.mode(), ee.Kernel.square(30,'meters'),'kernel',false)