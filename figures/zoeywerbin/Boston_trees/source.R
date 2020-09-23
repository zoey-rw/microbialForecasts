library(rgeos)
library(sf)

# if changed to "TRUE", you must have a copy of all other associated data
# email zrwerbin@bu.edu if you would like a copy
redownloadData = FALSE
if(redownloadData==TRUE){
  source("prepMapData.R")
  source("prepTreeData.R")
}

# read in prepped data files
tract_data <- readRDS("data/tract_data.rds")
tracts <- readRDS("data/tracts.rds")
app_df <- readRDS("data/treeData.rds")
parcel <- readRDS("data/parcel.rds")

# remove airport from dataset
tract_data <- tract_data[tract_data$TRACTCE10 != "981300",]

# Mean temperature
var1 <- tract_data$MEAN
pal1 <- colorNumeric("YlOrRd", domain=var1)

# Heat Vulnerability Index
var2 <- tract_data$HVI_cat
pal2 <- colorNumeric("YlOrRd", domain=var2)

# Current canopy, %
var3 <- tract_data$current_percent
pal3 <- colorNumeric("YlGn", domain=var3)

# Potential canopy, %
var4 <- tract_data$potential_percent
pal4 <- colorNumeric("YlGn", domain=var4)

# create ! icons
exclamation.point <-  makeIcon(iconUrl = "www/exclamation.png", iconWidth = 20, iconHeight = 30)


### a bunch of text for pop-ups ###
modal_sizing_text <- "<p>When planting a tree, you should take into account the expected canopy spread at the tree's maturity, as well as its height and the extent of its roots. Trees can upheave sidewalks or building foundations if planted too closely.</p> <p> <a href='https://www.arborday.org/trees/righttreeandplace/size.cfm'> Tree sizing guide from ArborDay.org </a></p>"

modal_allergen_text <- "<p>Allergic reactions (also known as 'hay fever') are often more frequent and more severe among sufferers of chronic asthma. If planting in a census tract with high chronic asthma prevalence, consider planting trees that are less likely to trigger allergic reactions. Click on a census tract in this application's interactive map to see the chronic asthma prevalence in your neighborhood. </p> <p>Individuals can have reactions to various types of plants, but allergies to airborne pollen from pine trees (rather than insect-pollinated trees) are common. Firs, cedars, and many fruit trees produce very low amounts of pollen, while elm, oak, zelkova and liquidambar produce higher amounts of pollen.</p> <p>Other resources for tree allergies:</p> <p> <a href='https://weather.com/forecast/allergy/l/USMA0046:1:US'> Allergy forecast from Weather.com </a></p>"

modal_site_text <- "Street trees have different requirements than those planted in parks or spacious yards. Street trees can create significant hazards when not maintained regularly, and they must be more tolerant to pollutants from cars (such as ozone). If power lines are present, smaller ornamental trees should be planted rather than larger shade trees."

modal_light_text <- "<a href='https://www.thespruce.com/what-is-full-sun-partial-shade-1402372'> Here are some guidelines from theSpruce.com. </a>"

modal_data_text <- "<p>Mean summer temperature is based on averages of mid-morning Land Surface Temperature (LST) from 2000-2013, calculated using remotely-sensed data. Data details and limitations (including a comparison with temperatures measured via sensors in Boston) are available from
<a href='https://doi.org/10.1175/JAMC-D-16-0325.1'> Wang (2017).</a></p>

<p> Census tract demographics used to calculate Heat Vulnerability Index (HVI), as well as prevalence of chronic asthma (used to recommend low-allergen trees on the 'Step 2' tab), come from 
<a href='https://www.census.gov/programs-surveys/acs/data/summary-file.2010.html'> the 2010 American Community Survey.</a>
 Full methods are in the publication linked below.</p>
<p>'Current' and 'potential' tree canopy data are derived from maps produced by 
<a href='https://www.uvm.edu/rsenr/sal/'> the University of Vermont Spatial Analysis Lab.</a> 'Potential' tree canopy refers to areas that already have grass cover. It does not include impervious surfaces, like concrete. </p>
<p>'Priority' regions are those that have high HVI or high mean summer temperatures, as well as high 'potential' for tree canopy.</p>
<p>Tree characteristics (on the 'Step 3' tab) come from the fact sheets linked for each tree.</p>
<p>Land ownership comes from the <a href='https://data.boston.gov/dataset/parcels-2016-data-full'>City of Boston's Parcels 2016 Data.</a> Property codes are categorized as follows:<br> City-owned property: 902 <br> State-owned property: 901, 910-929<br>Federal-owned property: 900<br>Boston Housing Authority, Boston Redevelopment Authority, <br>or other government-associated or public land: 903, 908, 965, 973, 978, 984, 986</p>
<br>More details in the preprint available here: [link]"

modal_alternatives_text <- "Green roofs and walls provide substantial cooling benefits, and can serve as educational resources or food gardens as well. The resources below can provide information and services related to installing and maintaining green roofs and walls in Boston:
<p> <a href='http://www.recovergreenroofs.com/'> Recover Green Roofs</a> design and install rooftop gardens and farms. They will also do maintenance for gardens.</p>
<p> <a href='https://greencitygrowers.com/'>Green City Growers</a> helps design and maintain rooftop gardens.</p>
<p> <a href='http://www.highergroundrooftopfarm.com/'>Higher Ground Rooftop Farm</a> helps maintain rooftop farms.</p>
<p> <a href='https://bostoncityscapes.com/'>Boston CityScapes</a> installs and maintains gardens and green walls.</p>
<p> <a href='https://urbanfarminginstitute.org/'>The Urban Farming Institute</a> plans events and workshops related to urban farming</p>
<p> <a href='https://www.freightfarms.com/ '>Freight Farms</a> installs farms inside of old freight shipping containers with a controlled environment</p>
<p> <a href='https://www.greenroofs.com/projects/'>Greenroofs.com</a> View examples of green roofs around the world!</p>"