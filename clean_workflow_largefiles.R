

#This script can be used to process very large .las or .laz files 
#The files I handled with this script were typically 10-12GB each 

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
#dplyr to amend and extend leafR functions to handle LAScatalog objects
library(dplyr)


#Read in a very large file as a LAScatalog object
#This is a representation of the file, but does not load the whole thing into R
ctg <- readLAScatalog("D:/post_03.las")
#plotting the object just allows you to visualize the boundaries of the point cloud
plot(ctg)

#These parameters specify how R should divide the large file into small chunks that can be read in and processed individually
#I set the chunk size to 10x10m. The LAScatalog engine will spit a warning about such small chunks being "irrelevant," but it will not terminate the function
opt_chunk_size(ctg) <- 10
#This allows you to set a custom buffer size. The default is 30m
opt_chunk_buffer(ctg) <- 5
#Drop anomalous false returns 
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 50"
#Specify a place for R to save the processed outputs, and condense this as a .laz (rather than .las) file
opt_output_files(ctg) <- "D:/ground/{XCENTER}_{YCENTER}"
opt_laz_compression(ctg) <- TRUE

#identify ground points 
las <- classify_ground(ctg, algorithm = csf())

#Now you should have a new folder of .laz files, each a processed 10x10m chunk of the full point cloud
#We can continue to process this folder as a new LAScatalog object

#The next step is to ground-normalize the point cloud
#Read in the "ground" folder
ctg <- readLAScatalog("D:/ground/")

#Specify "chunking" parameters and new output folder
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 50"
opt_laz_compression(ctg) <- TRUE
opt_output_files(ctg) <- "D:/ground_normal/{XCENTER}_{YCENTER}"

#normalize height
las <- lidR::normalize_height(ctg, algorithm = knnidw(rmax = 1))


#Again, using this new folder, we can start taking some structure metrics on the point cloud
ctg <- readLAScatalog("D:/ground_normal/")
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 50"
#No need to specify an output folder because this ouput from the "metrics" function is a small dataframe and can be handled in R

metrics <- pixel_metrics(ctg, func = ~myMetrics(z = Z, rn = ReturnNumber), 10)
metrics <- as.data.frame(metrics)

#The LAScatalog engine will add the calculated values for each chunk to a dataframe 
#To see the full-plot metrics, I typically averaged these chunk values
max(metrics$max)
mean(metrics$mean)
mean(metrics$std)
mean(metrics$entropy)
mean(metrics$vci)
mean(metrics$fhd)


#Now, we can make a canopy height model through the same process
#Again, we don't need to specify an output folder because the raster layer is small enough to be handled in R

ctg <- readLAScatalog("D:/post_01_normalized/")
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 50"

#Create canopy height model 
chm <- lidR::grid_canopy(ctg, res = 1, algorithm = p2r(na.fill = tin()))

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
chm2 <- rasterToPoints(chm, spatial = TRUE)
chm2 <- as.data.frame(chm)

#Top rugosity is the standard deviation of top of canopy height
sd(chm2$Z)
#And, it's easy to take the mean top of canopy height as well
mean(chm2$Z)

#Finally, here is the visualized canopy height model -- using ggplot2, which allows for complete control over the plot style
ggplot(chm, aes(x = x, y = y, fill = Z))+
  geom_raster()+
  scale_fill_viridis(limits = c(0, 50))+
  theme_classic() + xlab(NULL) + ylab(NULL) + 
  guides(fill=guide_legend(title="Height (m)"))


#Finally for plant area density, using leafR
#The leafR package is still in development and is not able to recognize or read LAScatalog objects
#I updated the specific functions I used to calculate plant area density and FHD to handle LAScatalogs

#First, I imported the leafR functions that process a typical .las file 

#Points by ZSlice is a function underlying plant area density calculations
pointsByZSlice = function(Z, maxZ){
  heightSlices = as.integer(Z) # Round down
  zSlice = data.table::data.table(Z=Z, heightSlices=heightSlices) # Create a data.table (Z, slices))
  sliceCount = stats::aggregate(list(V1=Z), list(heightSlices=heightSlices), length) # Count number of returns by slice
  
  colRange = 0:maxZ
  addToList = setdiff(colRange, sliceCount$heightSlices)
  n = length(addToList)
  if (n > 0) {
    bindDt = data.frame(heightSlices = addToList, V1=integer(n))
    sliceCount = rbind(sliceCount, bindDt)
    # Order by height
    sliceCount = sliceCount[order(sliceCount$heightSlices),]
    }
  
  colNames = as.character(sliceCount$heightSlices)
  colNames[1] = "ground_0_1m"
  colNames[-1] = paste0("pulses_", colNames[-1], "_", sliceCount$heightSlices[-1]+1, "m")
  metrics = list()
  metrics[colNames] = sliceCount$V1
  
  return(metrics)
}

#Next, I amended the "lad.voxels" function so that it can read in a LAScatalog object

lad.voxels <- function(las, grain.size = 1, k = 1)
{
  if (is(las, "LAScatalog"))  {
    options <- list(automerge = TRUE, need_buffer = TRUE)
    LAD_VOXELS <- catalog_apply(las, lad.voxels, grain.size = 1, k = 1)
    return(LAD_VOXELS)
  }
  else if (is(las, "LAScluster")) {
    bbox <- st_bbox(las)
    ##MAKING THE Boundary BOX into a dataframe thingy
    bbox <- as.numeric(bbox)
    pos <- c('xmin','ymin', 'xmax', 'ymax')
    nobuf <- data.frame(pos, bbox)
    
    las <- readLAS(las)
    if (is.empty(las)) return(NULL) 
    LAD_VOXELS <- lad.voxels(las, grain.size = 1, k = 1)
    
    ##make the LAD_VOXELinto a dataframe
    LAD_VOXELS <- as.data.frame(LAD_VOXELS)
    
    ##Subset out the buffer points 
    LAD_VOXELS <- subset(LAD_VOXELS, coordinates.X > nobuf$bbox[1])
    LAD_VOXELS <- subset(LAD_VOXELS, coordinates.X < nobuf$bbox[3])
    LAD_VOXELS <- subset(LAD_VOXELS, coordinates.Y > nobuf$bbox[2])
    LAD_VOXELS <- subset(LAD_VOXELS, coordinates.Y < nobuf$bbox[4])
    
    #Turn them back into items to use for the LAD profile 
    coordinates <- data.frame(LAD_VOXELS$coordinates.X, LAD_VOXELS$coordinates.Y)
    LAD <- subset(LAD_VOXELS, select = -c(coordinates.X, coordinates.Y))
    
    #Return them as a single list
    LAD_VOXELS <- list(coordinates, LAD)
    names(LAD_VOXELS) <- c("coordinates", "LAD")
    
    return(LAD_VOXELS)
  }
  else if (is(las, "LAS")) {
    #empty list object that will be fueling with binneds data.frames
    LAD_VOXELS = list()
    Z = NA
    
    #load normalized las cloud
    .las = las
    
    .las@data$Z[.las@data$Z < 0] = 0
    
    maxZ = floor(max(.las@data$Z))
    
    func = formula(paste0("~pointsByZSlice(Z, ", maxZ, ")"))
    t.binneds    = lidR::grid_metrics(.las, func, res = grain.size,
                                      start = c(min(.las@data$X), max(.las@data$Y)))
    t.binneds    = data.frame(sp::coordinates(t.binneds), raster::values(t.binneds))
    names(t.binneds)[1:2] = c("X", "Y")
    
    
    #getting the coordinates X and Y
    #t.binneds$X = coordinates(t.binneds)[,1]
    #t.binneds$Y = coordinates(t.binneds)[,2]
    #t.binneds = as.data.frame(t.binneds) #transforming in a data.frame
    
    #clip product by las files limits
    #t.binneds = t.binneds[t.binneds$X < xmax(.las) &
    #                        t.binneds$X > xmin(.las) &
    #                        t.binneds$Y > ymin(.las) &
    #                        t.binneds$Y < ymax(.las),]
    
    
    #select ground returns
    ground.returns = t.binneds[, grep("ground", names(t.binneds))]
    
    #select columns vegetation above 1m:
    if(nrow(t.binneds) != 1){ #this if is necessary when grain size is the whole plot
      pulses.profile.dz1 = t.binneds[, c(grep("pulses", names(t.binneds)))]
    }else{
      pulses.profile.dz1 = data.frame(matrix(as.numeric(as.character(t.binneds[, c(grep("pulses", names(t.binneds)))])), ncol = length(grep("pulses", names(t.binneds)))))
      names(pulses.profile.dz1) = names(t.binneds)[c(grep("pulses", names(t.binneds)))]
    }
    
    #invert data.frames for the sky be first
    pulses.profile.dz1 = pulses.profile.dz1[,length(pulses.profile.dz1):1] #invert columns
    
    #add grounds returns (0-1m)
    pulses.profile.dz1 = cbind(pulses.profile.dz1, ground.returns)
    rm(ground.returns)
    
    ### total matrix and cumsum.matrix:
    total.pulses.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, sum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1))
    cumsum.matrix.dz1 = matrix(apply(pulses.profile.dz1, 1, cumsum), ncol = length(pulses.profile.dz1), nrow = nrow(pulses.profile.dz1), byrow = TRUE)
    
    rm(pulses.profile.dz1)
    
    #Pulses out for each voxel
    pulse.out.dz1 = total.pulses.matrix.dz1 - cumsum.matrix.dz1
    
    #The pulses.out of voxel 1 is the pulses.in of voxel 2 and so on...
    #Therefore, pulse.in is pulse.out without the last line and adding in the
    #first line the total pulses:
    if(nrow(t.binneds) != 1){ #if used when grain size of the whole plot
      pulse.in.dz1 <- cbind(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
    }else{
      pulse.in.dz1 <- c(total.pulses.matrix.dz1[,1], pulse.out.dz1[,-c(ncol(pulse.out.dz1))])
    } #enf if
    
    rm(total.pulses.matrix.dz1, cumsum.matrix.dz1)
    
    # MacArthur-Horn eqquation
    # LAD = ln(S_bottom/S_top)*(1/(dz*K))
    #k value for LAD equation
    dz = 1
    
    LAD.dz1 = log(pulse.in.dz1/pulse.out.dz1) * 1/k * 1/dz
    
    rm(pulse.in.dz1, pulse.out.dz1)
    
    # Remove infinite and NaN values
    #Inf ocorre qndo pulses.out eh zero
    #NaN ocorre qndo pulses.in eh zero
    LAD.dz1[is.infinite(LAD.dz1)] <- NA; LAD.dz1[is.nan(LAD.dz1)] <- NA;
    
    #remove the first 1 meter close to the ground (and the ground too)
    LAD.dz1 = LAD.dz1[, -c(ncol(LAD.dz1))]
    
    #fuel list object
    LAD_VOXELS[["LAD"]] = LAD.dz1
    LAD_VOXELS[["coordinates"]] = t.binneds[,c("X", "Y")]
    
    rm(LAD.dz1, t.binneds)
    
    return(LAD_VOXELS)
  }
  else {
    stop("This type is not supported.")
  }
}


#The function can now read a LAScatalog

ctg <- readLAScatalog("D:/ground_normal/")
opt_chunk_buffer(ctg) <- 5
opt_filter(ctg) <- "-drop_z_below 0 -drop_z_above 50"

lad_vox <- lad.voxels(ctg, grain.size = 1, k = 1)

big_lad_vox <- bind_rows(lad_vox, .id = "column_label")
big_lad_vox <- subset(big_lad_vox, select = -c(1:2))


#This function was directly imported from the leafR package without modification 

lad.profile = function(VOXELS_LAD, relative = FALSE){
  if(relative == TRUE){
    t.lad.profile = apply(VOXELS_LAD$LAD, 2, mean, na.rm = TRUE)
    t.lad.profile = t.lad.profile/sum(t.lad.profile)*100
  }else{
    t.lad.profile = apply(VOXELS_LAD$LAD, 2, mean, na.rm = TRUE)
  }
  
  max_height = ncol(VOXELS_LAD[[1]]) + .5
  t.lad.profile = data.frame(height = seq(1.5, max_height), lad = t.lad.profile[length(t.lad.profile):1])
  return(t.lad.profile)
}


lad_profile_relative <- lad.profile(big_lad_vox, relative = TRUE)
lad_profile <- lad.profile(big_lad_vox)

#And finally, we can calculate our metrics of interest
leafR::FHD(lad_profile_relative)
sd(lad_profile_relative$lad)
leafR::lai(lad_profile)

#And plot the plant area density
plot(lad_profile_relative$height ~ lad_profile_relative$lad, type = "l",
     ylab = "Canopy height (m)", xlab = "LAD", ylim = c(0, 50), xlim = c(0, 20))

