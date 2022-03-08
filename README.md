# Annual global forest gain maps from 1984 to 2020
JS code and valid samples for mapping global forest gain years in Google Earth Engine.

In this study, a new framework was developed to map temporal and spatial changes of global forest gain using long-term Landsat observations on the GEE platform. First, stable forest regions during 1982-2020 were delineated based on the Global Forest Change (GFC) product. Stable forest samples were then collected using a random sampling strategy, and annual NBR thresholds for forest gain detection were generated based on these samples. Second, the NBR time series and LandTrendr were used to detect the forest gain years at the pixel level. Incorrect detection of forest gain arising from NBR fluctuations were reduced using Global Forest Canopy Height (GFCH), Finer Resolution Observation and Monitoring of Global Land Cover (FROM-GLC) and thresholds obtained in the first step.


## View online maps using GEE experimental app
https://duzhenrong.users.earthengine.app/view/globalforestgainyear
