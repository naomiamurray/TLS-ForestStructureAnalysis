
#This script can be used to process "small" .las or .laz files -- small enough to be read into R in their entirety 
#The largest file I processed this way was approx 1.5GB in size. 
#Anything bigger needs to be handled as a virtual LAScatalog object, or it will crash the program and local machine. 


#lidR and raster for general .las processing and standard metrics (like mean and max height)
library(lidR)
library(raster)
#RCSF for a cloth simulation filter to classify ground points and ground-normalize the point cloud
library(RCSF)
#ForestGapR for identifying canopy gaps 
library(ForestGapR)
#ggplot2 and viridis for customizing a visual canopy height model
library(ggplot2)
library(viridis)
#leafR for calculating plant area density and associated metrics (like foliage height diversity)
library(leafR)

#First, read in a las file
las <- lidR::readLAS("pre_01.las", select = "xyz",filter = "-drop_z_below 0 -drop_z_above 50")

#identify ground points 
las <- classify_ground(las, algorithm = csf())

#if you want to see how the algorithm identified the ground points, you can check with a plot
#plot(las, color = "Classification")

#Now that the ground has been identified, you can ground-normalize the point cloud
las <- lidR::normalize_height(las, algorithm = knnidw(rmax = 1))
#And clean up any anomalous, below-zero points 
las <- filter_poi(las, Z >= 0)

#Take some general metrics on the point cloud
#This will return max, mean, and standard deviation of height, % of height above the mean, along with a laundry list of other metrics
cloud_metrics(las, func = .stdmetrics)

#Trying some more custom metrics - first write a user-defined function 
myMetrics <- function(z, rn){
  first = rn == 1L
  zfirst = z[first]
  
  metrics = list(
    iqr =  IQR(z), # inter-quartile range
    entropy = entropy(z, by = 1, zmax = max(z)),
    vci = VCI(z, zmax = max(z)), # vertical complexity index
    fhd = (entropy(z, zmax = max(z)) * log(max(z))),  #foliar height diversity
    mean = mean(z),
    max = max(z),
    std = sd(zfirst))
  
  return(metrics)}

#Use your new function to generate your metrics on the point cloud
cloud_metrics(las, func = ~myMetrics(z = Z, rn = ReturnNumber), res = 100)


#The next set of metrics utilize a canopy height model, which is easy to generate with the lidR package
chm <- lidR::grid_canopy(las, res = 1, algorithm = p2r())

#If there are "empty" grids in the model, use a triangulated interpolation function to estimate missing values
#I rarely had to do this because TLS produces a very dense point cloud
#chm <- lidR::grid_canopy(las, res = 1, algorithm = p2r(na.fill = tin()))
#plot(chm)

#Rumple is calculated using this canopy height model 
rumple_index(chm)

#Use the ForestGapR package to identify canopy gaps
#This produces a binary raster layer of the gap areas
gaps <- getForestGaps(chm)

#And this will return stats on the gaps 
stats <- GapStats(gap_layer = gaps, chm_layer = chm)
#I specifically wanted the mean gap area and gap height
mean(stats$gap_area)
mean(stats$chm_mean)

#To take some final metrics on the canopy height model, and generate custom plots, turn the raster layer into a dataframe
chm <- rasterToPoints(chm, spatial = TRUE)
chm <- as.data.frame(chm)

#Top rugosity is the standard deviation of top of canopy height
sd(chm$Z)
#And, it's easy to take the mean top of canopy height as well
mean(chm$Z)

#Finally, here is the visualized canopy height model -- using ggplot2, which allows for complete control over the plot style
ggplot(chm, aes(x = x, y = y, fill = Z))+
  geom_raster()+
  scale_fill_viridis(limits = c(0, 50))+
  theme_classic() + xlab(NULL) + ylab(NULL) + 
  guides(fill=guide_legend(title="Height (m)"))


#Now for plant area density plots, using leafR
#The original leafR function wants a path to a normalized .las file, rather than the file itself
#So, I wrote the previously ground-normalized point cloud into a new file 
new.laz <- writeLAS(las, "C://Users/murra/Documents/Gabon_forest/new.las")

#which can then be called to estimate plant area density 
lad.vox <- lad.voxels("new.las", grain.size = 1, k = 1)
lad_profile <- lad.profile(lad.vox, relative = TRUE)
lad_profile_lai <- lad.profile(lad.vox)

#leafR also has a function to calculate FHD, I used the relative LAD for this, and for canopy rugosity
leafR::FHD(lad_profile)
#And I used the standard deviation of leaf area density as a measure of canopy rugosity. 
sd(lad_profile$lad)
#LAI has to be calculated with the absolute LAD
leafR::lai(lad_profile_lai)

#Finally, here is the plant area density plot. 
plot(lad_profile$height ~ lad_profile$lad, type = "l",
     ylab = "Canopy height (m)", xlab = "LAD", ylim = c(0, 50), xlim = c(0, 20))
