###################################################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection
## fitting PO PA and Integrated models the data
## 19/08/2016
###################################################################################################
######### Data required for the function ##########################################################
###################################################################################################
## All the data required for this function is stored in a file data.rda
## s.occupancy - raster with background covatiates that effect occupancy
## s.detection - raster with background covariates that effect detection
## po.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
## po.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
##
## y.pa - matrix of detection/non detection of the PA surveys
## PA.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
## PA.detection - matrix with covariates that effect detection in the locations of PA survey sites in each survey
###################################################################################################

require(raster)
require(fields)
require(mvtnorm)
require(matrixStats)

source("functions.r")

#Loading the data ----
load("data.rda")


#Checking whether occupancy and detection rasters have the same resolution -----
if(sum(res(s.occupancy)!=res(s.detection)))
	stop("Occupancy and detection raster layers have different resolution")

if(ncell(s.occupancy)!=ncell(s.detection))
	stop("Occupancy and detection have different number of cells")


# Plotting covariates that drive occupancy and detection in PO
ppi = 300
png('occupancy-covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.occupancy)
dev.off()

png('PO-detection -covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.detection)
dev.off()

# 1. Preparing the data ========================================================
#Preparing Presence-Only data ------------------------------------------------

#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))

# remove all NA values
tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell
s.area=area.back*nrow(X.back) #study area

# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(po.occupancy))), po.occupancy)
W.po=cbind(rep(1, nrow(as.matrix(po.detection))), po.detection)



#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.pa # matrix of presences and absences (when the species was and wasn't present)
J.pa=ncol(y.pa)
X.pa=cbind(rep(1, nrow(as.matrix(PA.occupancy))), PA.occupancy)
W.pa = array(dim=c(nrow(as.matrix(PA.detection)), J.pa, 2))
W.pa[,,1] = 1
W.pa[,,2] = PA.detection# if it changes

# 2. Analising the data ========================================================

#Analyzing Presence-Only data
po.fit=po.ipp(X.po, W.po,X.back, W.back)

# Analyzing presence-absence data
pa.fit=pa.mackenzie(X.pa,W.pa,y.pa)

# Analyzing presence-only data AND presence-absence data
poANDpa.fit=poANDpa.combined(X.po, W.po,X.back, W.back,X.pa,W.pa,y.pa)


