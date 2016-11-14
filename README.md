## List of the R-files supplied with the paper: Koshkina et al. "Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection"

POPA-function.r - Fitting PO PA and Integrated models the data. The data required for the file is stored in the file data.rda and include

	 - s.occupancy - raster with background covariates that effect occupancy
	 - s.detection - raster with background covariates that effect detection
	 - pb.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
	 - pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
	 - y.so - matrix of detection/non detection of the PA surveys
	 - so.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
	 - so.detection - matrix with covariates that effect detection in the locations of PA survey sites in 		each survey


sim-data.r  - Simulated data generated. Simulates covariates $x$ and $w$, PO and PA datasets.  All the data is saved in a file data.rda

functions.r  - Utility and likelihood functions required for POPA-function.r

