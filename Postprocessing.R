suppressMessages(library(rgeos))
suppressMessages(library(sf))
suppressMessages(library(sp))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(lwgeom))
# gdal_polygonizeR <- function(x, 
#                              outshape   = NULL, 
#                              gdalformat = 'ESRI Shapefile',
#                              pypath     = NULL, 
#                              readpoly   = TRUE, 
#                              quiet      = TRUE, 
#                              overwrite  = FALSE) {
#   
#   if (is.null(pypath)) {
#     pypath <- Sys.which('gdal_polygonize.py')
#   }
#   
#   if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
#   owd <- getwd()
#   on.exit(setwd(owd))
#   setwd(dirname(pypath))
#   if (!is.null(outshape)) {
#     outshape <- sub('\\.shp$', '', outshape)
#     f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep = '.'))
#     if (any(f.exists)) {
#       if (overwrite == FALSE) {
#         stop(sprintf('File already exists: %s',
#                      toString(paste(outshape, c('shp', 'shx', 'dbf'),
#                                     sep = '.')[f.exists])), call.=FALSE)
#       } else (
#         unlink(paste(outshape, c('shp', 'shx', 'dbf')))
#       )
#     }
#   }else outshape <- tempfile()
#   
#   if (methods::is(x, 'Raster')) {
#     
#     raster::writeRaster(x, {f <- tempfile(fileext = '.tif')})
#     rastpath <- normalizePath(f)
#   } else if (is.character(x)) {
#     rastpath <- normalizePath(x)
#   } else stop('x must be a file path (character string), or a Raster object.')
#   
#   system2('python', args = (paste0(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
#                                            pypath, rastpath, gdalformat, outshape), " -fieldname id")))
#   if (isTRUE(readpoly)) {
#     shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
#     return(shp)
#   }
#   return(NULL)
# }

postprocessing <- function(in.path, 
                           out.path, 
                           dk.thresh)
{
  polygon_sf <- suppressWarnings(st_read(in.path, quiet = TRUE))
  # douglas puecker algorithm
  polygon_sp_simplify <- suppressMessages(gSimplify(as(polygon_sf,"Spatial"), tol = dk.thresh, 
                                   topologyPreserve = TRUE))
  polygon_sf$geometry <- suppressMessages(st_as_sf(polygon_sp_simplify)) 
  polygon_sf_simplify <- suppressMessages(st_sf(data.frame(polygon_sf$id, geom = st_as_sf(polygon_sp_simplify))) %>% st_make_valid %>% st_collection_extract("POLYGON"))
  
  polygon_sfc_valid <- suppressMessages(polygon_sf_simplify[polygon_sf_simplify$polygon_sf.id!=0,])
  
  suppressMessages(st_write(polygon_sfc_valid, out.path, delete_layer = TRUE))
} 

# labelmap <- raster(file.path(paste0("/home/su/Documents/Jupyter/Segmetation_ARS_results/"
#                                     ,"Dollo_Ethiopia_watershed_2500_001.tif")))



arg <- commandArgs(TRUE)

if(length(arg) < 3) stop("At least 3 arguments needed", call. = FALSE)

in.path <- arg[1]
out.path <- arg[2]
dk.thresh <- as.numeric(arg[3])

# in.path <- '/media/su/DataFactory/MappingAfrica/segmentation/Comparison_site1/results/AOI4_70_95_wat_1200_01.shp'
# out.path <- '/media/su/DataFactory/MappingAfrica/segmentation/Comparison_site1/results/dk_000005.shp'
# area.thresh <- 5000
# dk.thresh <- 0.00005

if(file.exists(out.path)){
  file.remove(out.path)
}
  
suppressMessages(postprocessing(in.path = in.path, out.path = out.path, 
               		dk.thresh = dk.thresh))

