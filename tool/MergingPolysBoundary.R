# this script is used to merge polygons at intersecting lines across different tiles. 

library("aws.s3")
library("sf")
library("dplyr")

# the users need to double check the below four parameter
targeted_aoi <- '1' # your aoi id
aoi_file <- '/Users/coloury/Dropbox/MappingAfrica/segmentation_congo/congo_tiles_run1_v2.geojson' # tile geojson
dst.folder <- '/Users/coloury/Dropbox/MappingAfrica/segmentation_congo/polygon_merging' # the folder to store temp and final result
segmentation_prefix <- 'segmentation_congo' # the folder name for segmentation result in s3

# the below three can be changed in some special conditions
dsn.dir <- file.path(dst.folder, paste0("aoi", targeted_aoi, "_boundarymerge.geojson")) 
bucket <- 'activemapper'
tile_diam <- 0.005 / 200  * 2000

# this function is used to merge polygons only at intersecting lines
intersect_line <- function(x1, x2, y1, y2, poly.collect){
  p1 <- rbind(c(x1,y1), c(x1, y2), c(x2, y2), c(x2, y1), c(x1,y1))
  poly_line <- st_sfc(st_polygon(list(p1)))
  st_crs(poly_line) <- 4326
  polys_intersect <- poly.collect[st_intersects(poly.collect, poly_line, 
                                                sparse = FALSE),]
  
  # sort by boundary. This step is important because it determines that longer 
  # boundary will be processed as the first
  polys_boundary <- suppressMessages(st_intersection(poly.collect, poly_line))
  polys_intersect <- polys_intersect[order(st_area(polys_boundary), 
                                           decreasing = TRUE),]
  
  intersect.list <- suppressMessages(st_intersects(polys_intersect, polys_intersect))
  
  # rule of thumb: a poly is allowed to be merged for only once
  alreadymerge_id <- c() # global variables
  
  if(nrow(intersect.list) > 0){
    merge.polys <- lapply(1:nrow(intersect.list), function(x){
      if (!(x %in% alreadymerge_id)){
        # first, select those polys that are not merged ever
        candidates <- intersect.list[x][[1]][!(intersect.list[x][[1]] %in% alreadymerge_id)]
        if(length(candidates) > 1){
          # second, select that polys are from two different tiles
          candidates2 <- unlist(lapply(candidates, function(t){
            if(polys_intersect[x,]$tile != polys_intersect[t,]$tile){
              t
            }
          }))
          
          if(length(candidates2) > 0){
            # for multiple intersected polygon, choose the one which has the 
            # longest intersection line
            if(length(candidates2) > 1){
              # to.merge <- sample(candidates2, 1)
              # created a polyline centered at  
              primary.poly <- suppressMessages(st_buffer(polys_intersect[x, ], 
                                                         0.005 / 200))
              s <- lapply(1:length(candidates2), function(i){
                suppressMessages(st_intersection(primary.poly, polys_intersect[candidates2[i], ]))})
              s_area <- st_area(do.call(rbind, s))
              to.merge <- candidates2[which(s_area==max(s_area))]
            }else{
              to.merge <- candidates2
            }
            # print(paste0("tomerge", to.merge))
            alreadymerge_id <<- append(alreadymerge_id, to.merge)
            alreadymerge_id <<- append(alreadymerge_id, x)
            # print(alreadymerge_id)
            # print(polys_intersect)
            tmp <- suppressMessages(st_union(polys_intersect[x, ], polys_intersect[to.merge, ]))
            tmp %>% select(c('id', 'tile')) 
          }else{
            polys_intersect[x, ]
          }
        }else{
          polys_intersect[x, ]
        }
      }
    })
    
    merge.poly_whole <- do.call(rbind, merge.polys)
    
    # delete intsected polys
    poly.collect.new <- poly.collect[-polys_intersect$id, ]
    
    # add after-merge polys
    poly.collect <- rbind(poly.collect.new, merge.poly_whole[, c('id', 'tile')])
    
    return(poly.collect)
  }else{
    return(poly.collect)
  }
  
}

# # step 1: download
src.folder <- file.path(dst.folder, paste0("aoi", targeted_aoi, "_new"))
if(!dir.exists(src.folder)){
  dir.create(src.folder)
}
s3_url <- paste0('s3://', bucket, '/', segmentation_prefix)
download.cmd <- paste0('aws s3 cp ', s3_url,' ', src.folder,
                       ' --recursive')
system(download.cmd)

# # step 2: making original collection for aoi
lf <- list.files(src.folder)

poly.list <- lapply(1:length(lf), function(x){
  polys <- st_read(file.path(src.folder,
                             lf[x]))
  polys[, 'tile'] <- substring(lf[x], 5, 10)
  polys <- st_as_sf(polys)
  polys
})

poly.collect <- st_as_sf(do.call(rbind, poly.list)) %>% mutate(id=row_number())

dsn.dir <- file.path(dst.folder, paste0("aoi", targeted_aoi, "_unmerged.geojson"))
st_write(poly.collect, delete_dsn = TRUE, dsn.dir)

poly.collect <- st_read(dsn.dir) %>% mutate(id=row_number())

# read aoi definition geojson
aoi.boundary <- st_read(aoi_file)
aoi_bbox <- st_bbox(aoi.boundary[aoi.boundary$aoi==targeted_aoi,])
#aoi_bbox <- st_bbox(aoi.boundary[aoi.boundary$aois==targeted_aoi,])

nrow <- (as.numeric(aoi_bbox['ymax']) - as.numeric(aoi_bbox['ymin']))/ tile_diam - 1
ncol <- (as.numeric(aoi_bbox['xmax']) - as.numeric(aoi_bbox['xmin']))/ tile_diam - 1

# step 3: merging boundary per row and per column
for(nncol in 1: ncol){
  print(nncol)
  # adding a buffer of a pixel
  x1 = as.numeric(aoi_bbox['xmin']) + as.numeric(tile_diam) * nncol - 0.005 / 200
  x2 = as.numeric(aoi_bbox['xmin']) + as.numeric(tile_diam) * nncol + 0.005 / 200
  y1 = as.numeric(aoi_bbox['ymin']) 
  y2 = as.numeric(aoi_bbox['ymax']) 
  poly.collect <- intersect_line(x1=x1, x2=x2, y1=y1, y2=y2, poly.collect=poly.collect)
  poly.collect <- poly.collect %>% mutate(id=row_number())
}


for(nnrow in 1: nrow){
  print(nnrow)
  # adding a buffer of a pixel
  x1 = as.numeric(aoi_bbox['xmin']) 
  x2 = as.numeric(aoi_bbox['xmax']) 
  y1 = as.numeric(aoi_bbox['ymin']) + as.numeric(tile_diam) * nnrow - 0.005 / 200
  y2 = as.numeric(aoi_bbox['ymin']) + as.numeric(tile_diam) * nnrow + 0.005 / 200
  poly.collect <- intersect_line(x1=x1, x2=x2, y1=y1, y2=y2, poly.collect=poly.collect)
  poly.collect <- poly.collect %>% mutate(id=row_number())
}

# st_write(merge.poly_whole, delete_dsn = TRUE, '/Users/coloury/Dropbox/MappingAfrica/segmentation/polygon_merging/aoi8_intersect_union.geojson')
# st_write(polys_intersect, delete_dsn = TRUE, '/Users/coloury/Dropbox/MappingAfrica/segmentation/polygon_merging/aoi8_intersect_origin.geojson')

st_write(poly.collect, delete_dsn = TRUE, dsn.dir)


# # upload
# s3_url <- paste0('s3://', bucket, '/segmentation_V2')
# upload.cmd <- paste0('aws s3 cp ',file.path(dst.folder, paste0("aoi", targeted_aoi)), ' ', s3_url,
#                        ' --recursive')

# system(upload.cmd)

