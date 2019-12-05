suppressMessages(library(yaml))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(units))
suppressMessages(library(sf))
suppressMessages(library(DBI))

res_perpixel <- 2.5 # approximate resolution in GCS of 3-m resolution

point_to_gridpoly <- function(xy, w, OldCRSobj, NewCRSobj) {
  dw <- list("x" = c(-w, w, w, -w, -w), "y" = c(w, w, -w, -w, w))
  pols <- do.call(rbind, lapply(1:nrow(xy), function(i) {  # i <- 1
    xs <- unlist(sapply(dw$x, function(x) unname(xy[i, "x"] + x)))
    ys <- unlist(sapply(dw$y, function(x) unname(xy[i, "y"] + x)))
    p1 <- list(t(sapply(1:length(xs), function(i) c(xs[i], ys[i]))))
    ## sf, sfc, sfg are three basic classes used to repreesent simple featurtes in sf package 
    ## see https://cran.r-project.org/web/packages/sf/vignettes/sf1.html
    ## create a geometry (sfg) from a list of points, e.g., point, polygon, multipolygon
    pol <- st_polygon(p1)
    ## create a sfc, which is a list of sfg
    poldf <- st_sfc(pol)
    ## create a sf, a table which contains feature atributes and feature geometries (sfc)
    ## .() is actually just an alias to ‘list()’. It returns a data table, whereas not using ‘.()’ only returns a vector
    polsf <- st_sf(xy[i, .(name)], geom = poldf)
    st_crs(polsf) <- OldCRSobj # first set GCS
    polsf <- st_transform(polsf, crs = NewCRSobj) # then transform into PRS
    polsf
  }))
  return (pols)
}

preprocessing <- function(
  aoi_id
)
{
  # tile.aoi <- '/home/ubuntu/source/ghana_tiles_merged.geojson'
  # yaml.dir <- '/home/ubuntu/source/segmenter_config.yaml'
  # The below is for local run
  workingfolder <- '/home/su/Dropbox/MappingAfrica/segmentation/source/'
  # workingfolder <- '/home/ubuntu/source/'
  tile.aoi <- paste0(workingfolder, 'ghana_tiles_merged.geojson')
  yaml.dir <- paste0(workingfolder, 'segmenter_config.yaml')
  # decide cluster id
  if(aoi_id < 7){
    cluster_id <- 1
  }else if (aoi_id == 7 || aoi_id == 8 || aoi_id == 9|| aoi_id == 12 || aoi_id == 15){
    cluster_id <- 2
  }else{
    cluster_id <- 3
  }
    
  static.cluster.dir <- paste0(workingfolder, 'incoming_names_static_cluster', cluster_id,'.csv')
  qsite <- TRUE
  params <- read_yaml(yaml.dir)
  
  uname <- "postgis"
  dbname <- 'Africa'
  password <- params$mapper$db_password
  # print(password)
  # print(aoi_id)
  host <- paste0('labeller', aoi_id, '.crowdmapper.org')
  # host <- 'ec2-3-229-217-195.compute-1.amazonaws.com'
  con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), host = host,
                        dbname = dbname,
                        user = uname,
                        password = password)
  
  aoi_merge <- suppressMessages(st_read(tile.aoi))
  
  diam <- 0.005 / 2 ## new master grid diameter
  prjsrid <- 102022
  prjstr <- as.character(tbl(con, "spatial_ref_sys") %>% 
                           filter(srid == prjsrid) %>% 
                           dplyr::select(proj4text) %>% collect())
  
  gcsstr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  # step 1: query grid name in database
  label.sql <- paste0("select name from incoming_names where processed = 'TRUE'", 
                      "and label = 'TRUE'")
  label.dynamic <- suppressWarnings(DBI::dbGetQuery(con, label.sql))$name
  
  # step 2: query grid name in static csv within aoi
  aoi <- aoi_merge  %>% filter(production_aoi==aoi_id)
  aoi.union <- st_union(aoi)
  
  static_label <- as.character(read.csv(static.cluster.dir)$name)
  
  label_within <- lapply(1:length(static_label), function(x){
    kmlid <- static_label[x]
    xy_tabs <- data.table(tbl(con, "master_grid") %>% 
                            filter(name == kmlid) %>% 
                            dplyr::select(x, y, name) %>% collect())
    #dplyr::select(x, y, name) %>% collect())
    
    grid.poly <- point_to_gridpoly(xy = xy_tabs, w = diam, NewCRSobj = gcsstr, 
                                   OldCRSobj = gcsstr)
    
    if(nrow(suppressWarnings(st_intersection(grid.poly, aoi.union))) > 0)
    {
      kmlid
    }
  })
  
  label.static <- as.character(
    do.call(rbind, label_within[!is.na(label_within)])
  )
  
  label.merge <- c(label.dynamic, label.static)
  
  # query all worker ids
  workerid.sql <- paste0("select worker_id from worker_data")  
  all_workerid <- (DBI::dbGetQuery(con, workerid.sql))$worker_id
  
  worker.score <- lapply(1:length(all_workerid), function(x){
    # read all scored assignments including 'Approved' and 'Rejected'
    # for calculating history field and no field likelihood of worker i
    
    ############################# assignment history #########################
    histassignid.sql <- paste0("select assignment_id from", 
                               " assignment_data where worker_id = '", 
                               all_workerid[x], "'",
                               " and (status = 'Approved' OR status = ",
                               "'Rejected') and score IS NOT NULL")
    historyassignmentid <- suppressWarnings(
      DBI::dbGetQuery(con, histassignid.sql)$assignment_id)
    
    if (is.null(historyassignmentid)){
      userhistories <- NA
    } else {
      # query all valid likelihood and score from history assignments
      userhistories <- lapply(c(1:length(historyassignmentid)), function(x){
        
        # query likelihood and scores
        likelihood.sql <- paste0("select new_score, field_skill, nofield_skill",
                                 " from accuracy_data  where assignment_id='", 
                                 historyassignmentid[x], "'") 
        measurements <- suppressWarnings(DBI::dbGetQuery(con, 
                                                         likelihood.sql))
        
        # field_skill and nofield_skill are alias of max.field.lklh 
        # and max.nofield.lklh 
        if (nrow(measurements) != 0) {
          c('ml.field' = as.numeric(measurements$field_skill), 
            'ml.nofield' = as.numeric(measurements$nofield_skill), 
            'score.hist' = as.numeric(measurements$new_score))
        } else {
          NA
        }
      })
    }
    
    ##########################################################################
    # combine regular assignment hisotry for the user and delete NA values
    if(length(userhistories[!is.na(userhistories)]) != 0) {
      userhistories <- data.frame(
        do.call(rbind, userhistories)
      )
    } else {
      userhistories <- data.frame(
        'ml.field' = NA, 'ml.nofield' = NA, 'score.hist' = NA)
    }
    
    # calculating mean max likelihood and score from qual and non-qual history
    score.hist <- mean(userhistories$score.hist, 
                       na.rm = TRUE)
    if(is.na(score.hist) == TRUE)
    {
      score.hist = 0
    }
    c('workid' = all_workerid[x], 'score' = score.hist)
  })
  
  worker.score.dt <- as.data.table(suppressWarnings(do.call(rbind, worker.score)))
  
  # only use dynamic cells that are within each instance (not cluster instance)
  bayes.polys.collection <- lapply(1:length(label.dynamic), function(ii){
    # print(ii)
    assignment.sql <- paste0("select assignment_id from assignment_data",
                             " INNER JOIN hit_data USING (hit_id)",
                             " where name ='", label.dynamic[ii], 
                             "' and status = 'Approved'", 
                             " order by assignment_id")  
    assignmentid <- (DBI::dbGetQuery(con, assignment.sql))$assignment_id
    
    assignmentid <- unlist(assignmentid)
    
    # query polygons for each assignment
    bayes.polys <- lapply(1:length(assignmentid), function(x) {
      
      # workerid
      workerid.sql <- paste0("select worker_id from assignment_data", 
                             " where assignment_id ='", assignmentid[x], 
                             "' order by assignment_id")
      workerid <- (DBI::dbGetQuery(con, workerid.sql))$worker_id
      
      # read all scored assignments including 'Approved' and 'Rejected'
      # for calculating history field and no field likelihood of worker i
      
      if(params$mapper$mapping_category1 == 'field') {
        # read'field-group' polygons
        user.sql <- paste0("select name, geom_clean ",
                           "FROM user_maps INNER JOIN categories ", 
                           "USING (category) where assignment_id='",  
                           assignmentid[x], "' ", "AND categ_group ='field'")
        user.polys <- suppressWarnings(
          DBI::dbGetQuery(con, gsub(", geom_clean", "", user.sql))
        )
        
        
        user.hasfields <- ifelse(nrow(user.polys) > 0, "Y", "N")
        
        # if user maps have field polygons
        if(user.hasfields == "Y") {
          user.polys <- suppressWarnings(st_read(con, query = user.sql))
          
          # select only polygons
          user.polys <- user.polys %>% filter(st_is(. , "POLYGON"))
          
          # union user polygons
          user.poly <- suppressWarnings(suppressMessages(
            st_buffer(user.polys, 0)))
          
          # if for F sites, we need to first intersection user maps by grid 
          # to retain those within-grid parts for calculation
          
          if(length(user.poly) == 0){
            geometry.user <- st_polygon()
          } else {
            geometry.user <-  user.poly$geom_clean
          }
        }
        else {
          # if users do not map field, set geometry as empty polygon
          geometry.user <- st_polygon()
        }  
      } 
      
      # we give 0.5 as posterior probability to unsure, meaning that the user
      # thinks it has only 50% to be a field
      # bayes.poly will consist two sf rows, the first is that the surely-labeled
      # fields, and the second is that unsure fields
      bayes.poly <- st_sf('prior' = worker.score.dt[worker.score.dt$workid == workerid]$score, 
                          geometry = st_sfc(geometry.user))
      
      # set crs
      st_crs(bayes.poly) <- gcsstr
      bayes.poly
      
    }) # end bayes.polys <- lapply(assignmentid, function(x)
    
    bayes.polys <- suppressWarnings(do.call(rbind, bayes.polys))
    
    bayes.polys <- bayes.polys %>% filter(prior==max(bayes.polys$prior))
    bayes.polys
  })
  
  bayes.polys.collection.bind <- suppressWarnings(do.call(rbind, bayes.polys.collection))
  
  # select those non-zero polygons for area estimation
  valid.polys <- st_area(bayes.polys.collection.bind)[as.numeric(st_area(bayes.polys.collection.bind)) > 0]

  min.thres <- sort(valid.polys)[round(length(valid.polys) * 0.01)]
  
  max.thres <- sort(valid.polys)[round(length(valid.polys) * 0.99)]
  
  params$segmenter$mmu <- min.thres / (res_perpixel * res_perpixel)
  
  params$segmenter$max_field_size <- max.thres / (res_perpixel * res_perpixel)
  
  # write new yaml
  write_yaml(params, yaml.dir)
  
}

options(warn=-1)

arg <- commandArgs(TRUE)

if(length(arg) < 1) stop("At least 1 arguments needed", call. = FALSE)

aoi_id <- arg[1]

suppressMessages(preprocessing(aoi_id = aoi_id))

# hist(as.numeric(st_area(bayes.polys.collection.bind)), 
#      breaks = seq(0, as.numeric(max(st_area(bayes.polys.collection.bind))), length.out = 200), 
#      xlab  = 'polygon area (m^2)', ylab = 'frequency')
