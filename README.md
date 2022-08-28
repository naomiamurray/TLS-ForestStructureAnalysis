# TLS-ForestStructureAnalysis
Code to calculate common forest structure metrics on TLS point clouds

Both scripts can be used to generate the following metrics:

	- Max height
	- Mean height
	- Standard deviation of height
	- Number of canopy gaps
	- Mean gap size
	- Proportion of gap area
	- Entropy
	- Foliage height diversity (FHD)
	- Top rugosity
	- Canopy rugosity
	- Rumple
	- Plant area density (PAD)
	- Plant area index (PAI) 

 The "smallfiles" script is intended to be used with .las files that can be read into the R environment in their entirety (<2 GB). The "largefiles" script can be used with .las files that are too big to be held within R -- it utilizes the LAScatalog engine to read them in chunks.
