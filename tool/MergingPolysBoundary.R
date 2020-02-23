# This script is used to download segmentation geojson, and merge polygons at 
# the boundary. The output is a combined per-aoi geojson. 
# This script is applied after the segmentation production is finished
# Author: Su Ye

library("aws.s3")
library("sf")
library("dplyr")

bucket <- 'activemapper'
targeted_aoi <- '3' # change it as needed
aoi_file <- 'ghana_tiles_merged.geojson'
segmentation_prefix <- 'segmentation_V2'
tile_diam <- 0.005 / 200  * 2000
download.needed <- FALSE
# change folder name as your setting
dst.folder <- '/Users/coloury/Dropbox/MappingAfrica/segmentation/polygon_merging'
src.folder <- file.path(dst.folder, paste0("aoi", targeted_aoi))
if(!dir.exists(src.folder)){
  dir.create(src.folder)
}

# read aoi definition geojson
aoi_extent <- st_read(file.path(dst.folder, aoi_file))
aoi_extent_selected <- aoi_extent %>% filter(production_aoi == targeted_aoi)
# step 1: download
if (download.needed == True){
  for(i in aoi_extent_selected$tile){
    s3_url <- paste0('s3://', bucket, '/', segmentation_prefix, '/tile', i, '_seg.geojson')
    download.cmd <- paste0('aws s3 cp ', s3_url,' ', file.path(dst.folder, paste0("aoi", targeted_aoi)))
    system(download.cmd)
  }
}


# step 2: making original collection for aoi
lf <- list.files(src.folder)

if(!file.exists(file.path(dst.folder, paste0("aoi", targeted_aoi, "_unmerged.geojson")))){
  poly.list <- lapply(1:length(lf), function(x){
    polys <- st_read(file.path(file.path(dst.folder, paste0("aoi", targeted_aoi)),  
                               lf[x]))
    polys[, 'tile'] <- substring(lf[x], 5, 10)
    polys <- st_as_sf(polys)
    polys
  })

  # all polygons for aoi were collected as sf 
  poly.collect <- st_as_sf(do.call(rbind, poly.list)) %>% mutate(id=row_number())

  # check which tiles polygons are missing from downloading. The missing tiles might 
  # be two reasons: 1) failure to process due to unknown bug; 2) no field in tile
  # the operator needs to compare the missing tile list with the log file produced for each
  # production, and found those tiles that fails to be recorded in the logging system, 
  # these tiles should be due to some unknown bugs that needs to be double checked
  outputted.tiles <- unique(poly.collect$tile)
  original.tiles <- unique(aoi_extent_selected$tile)
  # original.tiles[!(original.tiles %in% outputted.tiles)]
  print(paste0("The missing tiles for aoi ",  targeted_aoi, " are: ", 
               paste(original.tiles[!(original.tiles %in% outputted.tiles)], 
                                                collapse = ' ')))
  # save unmerged
  dsn.dir <- file.path(dst.folder, paste0("aoi", targeted_aoi, "_unmerged.geojson"))
  st_write(poly.collect, delete_dsn = TRUE, dsn.dir)
}else{
  poly.collect <- st_read(file.path(dst.folder, paste0("aoi", targeted_aoi, "_unmerged.geojson")))
}

# st_write(poly.collect, '/Users/coloury/Dropbox/MappingAfrica/segmentation/polygon_merging/aoi8_collect.geojson')
aoi_bbox <- st_bbox(aoi_extent_selected)

nrow <- (as.numeric(aoi_bbox['ymax']) - as.numeric(aoi_bbox['ymin']))/ tile_diam
ncol <- (as.numeric(aoi_bbox['xmax']) - as.numeric(aoi_bbox['xmin']))/ tile_diam

intersect_line <- function(x1, x2, y1, y2, poly.collect){
  p1 <- rbind(c(x1,y1), c(x1, y2), c(x2, y2), c(x2, y1), c(x1,y1))
  poly_line <- st_sfc(st_polygon(list(p1)))
  st_crs(poly_line) <- 4326
  polys_intersect <- poly.collect[st_intersects(poly.collect, poly_line, 
                                                sparse = FALSE),]
  
  # sort by boundary. This step is important because it determines that longer 
  # boundary will be processed as the first
  polys_boundary <- st_intersection(poly.collect, poly_line)
  polys_intersect <- polys_intersect[order(st_area(polys_boundary), 
                                           decreasing = TRUE),]
  
  intersect.list <- st_intersects(polys_intersect, polys_intersect)
  
  # rule of thumb: a poly can be merged only once
  alreadymerge_id <- c() # global variables
  merge.polys <- lapply(1:nrow(intersect.list), function(x){
    if (!(x %in% alreadymerge_id)){
      # first select those polys that are not merged ever
      candidates <- intersect.list[x][[1]][!(intersect.list[x][[1]] %in% alreadymerge_id)]
      if(length(candidates) > 1){
        # second select that polys are from two different tiles
        candidates2 <- unlist(lapply(candidates, function(t){
          if(polys_intersect[x,]$tile != polys_intersect[t,]$tile){
            t
          }
        }))
        
        if(length(candidates2) > 0){
          # for multiple intersected polygon, choose the longest intersection line
          if(length(candidates2) > 1){
            # to.merge <- sample(candidates2, 1)
            primary.poly <- st_buffer(polys_intersect[x, ], 0.005 / 200)
            s <- lapply(1:length(candidates2), function(i){
              st_intersection(primary.poly, polys_intersect[candidates2[i], ])})
            s_area <- st_area(do.call(rbind, s))
            to.merge <- candidates2[which(s_area==max(s_area))]
          }else{
            to.merge <- candidates2
          }
          alreadymerge_id <<- append(alreadymerge_id, to.merge)
          alreadymerge_id <<- append(alreadymerge_id, x)
          # print(alreadymerge_id)
          # print(polys_intersect)
          tmp <- suppressWarnings(st_union(polys_intersect[x, ], polys_intersect[to.merge, ]))
          tmp %>% select(c('id', 'tile')) 
        }else{
          # has intersected polygons but no paired with two tile ids
          polys_intersect[x, ]
        }
      }else{
        # no intersected polygons
        polys_intersect[x, ]
      }
    }
  })
  
  merge.poly_whole <- do.call(rbind, merge.polys)
  
  # delete intsected polys
  # poly.collect <- poly.collect[-polys_intersect$id]
  poly.collect <- poly.collect[!(poly.collect$id %in% polys_intersect$id), ]
  
  # add after-merge polys
  poly.collect <- rbind(poly.collect, merge.poly_whole[, c('id', 'tile')])
  
  return(poly.collect)
}

# step 3: merging boundary per row and per column
for(nncol in 1: ncol){
  # adding a buffer of a pixel
  x1 = as.numeric(aoi_bbox['xmin']) + as.numeric(tile_diam) * nncol - 0.005 / 200
  x2 = as.numeric(aoi_bbox['xmin']) + as.numeric(tile_diam) * nncol + 0.005 / 200
  y1 = as.numeric(aoi_bbox['ymin']) 
  y2 = as.numeric(aoi_bbox['ymax']) 
  poly.collect <- suppressMessages(intersect_line(x1=x1, x2=x2, y1=y1, y2=y2, 
                                                  poly.collect=poly.collect))
}

for(nnrow in 1: nrow){
  # adding a buffer of a pixel
  x1 = as.numeric(aoi_bbox['xmin']) 
  x2 = as.numeric(aoi_bbox['xmax']) 
  y1 = as.numeric(aoi_bbox['ymin']) + as.numeric(tile_diam) * nnrow - 0.005 / 200
  y2 = as.numeric(aoi_bbox['ymin']) + as.numeric(tile_diam) * nnrow + 0.005 / 200
  poly.collect <- suppressMessages(intersect_line(x1=x1, x2=x2, y1=y1, y2=y2, 
                                                  poly.collect=poly.collect))
}

dsn.dir <- file.path(dst.folder, paste0("aoi", targeted_aoi, "_boundarymerge.geojson"))
st_write(poly.collect, delete_dsn = TRUE, dsn.dir)


# # upload
# s3_url <- paste0('s3://', bucket, '/segmentation_V2')
# upload.cmd <- paste0('aws s3 cp ',file.path(dst.folder, paste0("aoi", targeted_aoi)), ' ', s3_url,
#                        ' --recursive')

# system(upload.cmd)

