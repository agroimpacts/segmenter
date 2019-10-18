# This script is to implement segmentation based on probability images and original composite imagies
# author: Su Ye

import boto3
import pandas as pd
import geopandas as gpd
import gdal
import yaml
from yaml import Loader
from datetime import datetime
from shapely.geometry import mapping
import logging
import os
import numpy as np
from osgeo import gdal_array
from skimage.filters import sobel
from sklearn import preprocessing
from skimage.segmentation import slic, watershed
from osgeo import ogr
import cv2 as cv
from skimage.future import graph
import subprocess
from osgeo import osr
from sklearn import preprocessing
from fixed_thread_pool_executor import FixedThreadPoolExecutor
import time
import click
import multiprocessing
from pytz import timezone
from botocore.exceptions import ClientError

def run_cmd(cmd, logger):
    """
    using os to run a command line
    arg:
        cmd: a command line
        logger: handler of logging file
    """
    try:
        os.system(cmd)
    except OSError as e:
        logger.error("Runing command line '{}' fails: {}".format(cmd, e))
        raise


def is_valid_image(path):
    """
    check if the image is valid or not
    arg:
        path: the path of image to check
    """
    ds = gdal.Open(path)
    if ds is None:
        return False

    rasterArray = np.array(ds.GetRasterBand(1).ReadAsArray())
    unique_val = np.unique(rasterArray)
    if len(unique_val) == 1:
        del ds
        return False
    else:
        del ds
        return True


def parse_catalog_from_s3(bucket, prefix, catalog_name):
    """
    read bucket, prefix from yaml.
    arg:
        bucket: Name of the S3 bucket.
        prefix: prefix for yaml file
        catalog_name: name of catalog file
    return:
        'catalog' pandas object
    """
    s3 = boto3.client('s3')
    obj = s3.get_object(Bucket=bucket, Key='{}/{}'.format(prefix, catalog_name))
    catalog = pd.read_csv(obj['Body'], sep=" ")
    return catalog

def parse_yaml_from_s3(bucket, prefix):
    """
    read bucket, prefix from yaml.
    arg:
        bucket: Name of the S3 bucket.
        prefix: the name for yaml file
    return:
        yaml object
    """
    s3 = boto3.resource('s3')
    obj = s3.Bucket(bucket).Object(prefix).get()['Body'].read()
    return yaml.load(obj)

def get_colrow_geojson(foc_gpd_tile, left_corner_x, left_corner_y, per_tile_width):
    """
    get col and row based on inputted gpd object
    arg:
        foc_gpd_tile: the geopandas object as inputted
        left_corner_x: the x coordinate of left corner for the coarse layout (see learner for definition of layer out)
        left_corner_y: the y coordinate of left corner for the coarse layout
        per_tile_width: the width of tile (in GCS)
    return:
        (col, row)
    """
    extent_geojson = mapping(foc_gpd_tile['geometry'])
    center_x = (extent_geojson['bbox'][0] + extent_geojson['bbox'][2]) / 2
    center_y = (extent_geojson['bbox'][1] + extent_geojson['bbox'][3]) / 2
    row = int(round((left_corner_y - center_y +  per_tile_width / 2) / per_tile_width)) - 1
    col = int(round((center_x - left_corner_x + per_tile_width / 2) / per_tile_width)) - 1
    return col, row


def weight_mean_color(graph, src, dst, n):
    """Callback to handle merging nodes by recomputing mean color.

    The method expects that the mean color of `dst` is already computed.

    Parameters
    ----------
    graph : RAG
        The graph under consideration.
    src, dst : int
        The vertices in `graph` to be merged.
    n : int
        A neighbor of `src` or `dst` or both.

    Returns
    -------
    data : dict
        A dictionary with the `"weight"` attribute set as the absolute
        difference of the mean color between node `dst` and `n`.
    """

    diff = graph.node[dst]['mean color'] - graph.node[n]['mean color']
    diff = np.linalg.norm(diff)
    return {'weight': diff}


def merge_mean_color(graph, src, dst):
    """Callback called before merging two nodes of a mean color distance graph.

    This method computes the mean color of `dst`.

    Parameters
    ----------
    graph : RAG
        The graph under consideration.
    src, dst : int
        The vertices in `graph` to be merged.
    """
    graph.node[dst]['total color'] += graph.node[src]['total color']
    graph.node[dst]['pixel count'] += graph.node[src]['pixel count']
    graph.node[dst]['mean color'] = (graph.node[dst]['total color'] /
                                     graph.node[dst]['pixel count'])


def gdal_save_file_tif_3bands(out_path, r, g, b, gdal_type, trans, proj, rows, cols):
    """
    save file
    Parameters
    ----------
    out_path : full outputted path
    r : the first band (numpy array)
    g : the second band (numpy array)
    b : the third band (numpy array)
    gdal_type: gdal type
    trans: transform coefficients
    proj: projection
    Returns
    -------
    TRUE OR FALSE
    """
    outdriver = gdal.GetDriverByName("GTiff")
    outdata = outdriver.Create(out_path, rows, cols, 3, gdal_type)
    if outdata is None:
        return False
    outdata.GetRasterBand(1).WriteArray(r)
    outdata.FlushCache()
    outdata.GetRasterBand(2).WriteArray(g)
    outdata.FlushCache()
    outdata.GetRasterBand(3).WriteArray(b)
    outdata.FlushCache()
    outdata.SetGeoTransform(trans)
    outdata.FlushCache()
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None
    return True


def gdal_save_file_tif_1bands(out_path, array, gdal_type, trans, proj, rows, cols):
    """
    save file
    Parameters
    ----------
    out_path : full outputted path
    array : numpy array to be saved
    gdal_type: gdal type
    trans: transform coefficients
    proj: projection
    Returns
    -------
    TRUE OR FALSE
    """
    outdriver = gdal.GetDriverByName("GTiff")
    outdata = outdriver.Create(out_path, rows, cols, 3, gdal_type)
    if outdata is None:
        return False
    outdata.GetRasterBand(1).WriteArray(array)
    outdata.FlushCache()
    outdata.SetGeoTransform(trans)
    outdata.FlushCache()
    outdata.SetProjection(proj)
    outdata.FlushCache()
    outdata = None
    return True


def segmentation_season(tile_id, season, uri_composite_gdal, uri_prob_gdal, working_dir, mmu, prob_threshold, buf,
                        logger, verbose):
    """
    segmentation for signle season
    tile_id: tile id
    season: season
    uri_composite_gdal: s3 url for composite image
    uri_prob_gdal: s3 url for probability image
    working_dir: outputted folder
    mmu: minimum mapping unit
    prob_threshold: threshold for probability
    buf: buffer of composite image (the default is 11 pixels)
    logger: logger
    verbose: True or False, if outputted intermediate results
    return: True or False
    """

    proj = ''

    # a quick way to read image as numpy array
    array_composite = gdal_array.LoadFile(uri_composite_gdal)
    array_prob = gdal_array.LoadFile(uri_prob_gdal)

    # get channels, cols, rows
    [nchannels, cols, rows] = array_composite.shape

    # temporal change for dealing with Ron's Rf probability image
    # array_prob = array_prob[buf:cols - buf, buf:rows - buf]

    # STEP 1
    # meanshift algorithm to smooth image and filter out noise
    B1, b2, b3, b4 = array_composite

    # scale to int8 for opencv processing
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 255))

    # get min and max
    b2_scale = min_max_scaler.fit_transform(b2.astype(np.float64).reshape(cols * rows, 1)).astype(np.uint8)
    b3_scale = min_max_scaler.fit_transform(b3.astype(np.float64).reshape(cols * rows, 1)).astype(np.uint8)
    b4_scale = min_max_scaler.fit_transform(b4.astype(np.float64).reshape(cols * rows, 1)).astype(np.uint8)

    # keep records for min and max for original bands
    b2_min = np.min(b2)
    b2_max = np.max(b2)
    b3_min = np.min(b3)
    b3_max = np.max(b3)
    b4_min = np.min(b4)
    b4_max = np.max(b4)

    mat_norm = cv.merge([b2_scale.reshape(cols, rows), b3_scale.reshape(cols, rows), b4_scale.reshape(cols, rows)])
    
    # meanshift algorithm
    dst_norm = cv.pyrMeanShiftFiltering(mat_norm, 5, 5, termcrit=(cv.TERM_CRITERIA_EPS | cv.TERM_CRITERIA_COUNT, 1, 5))

    # rescale to int16 (i.e., original band data range)
    b2_m, b3_m, b4_m = cv.split(dst_norm)
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(b2_min, b2_max))
    b2_filter = min_max_scaler.fit_transform(b2_m.reshape(cols * rows, 1).astype(np.float32))
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(b3_min, b3_max))
    b3_filter = min_max_scaler.fit_transform(b3_m.reshape(cols * rows, 1).astype(np.float32))
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(b4_min, b4_max))
    b4_filter = min_max_scaler.fit_transform(b4_m.reshape(cols * rows, 1).astype(np.float32))

    # mat_filter_rescale = cv.merge([b2_filter.reshape(cols, rows),
    #                    b3_filter.reshape(cols, rows),
    #                    b4_filter.reshape(cols, rows)])
    # r, g, b = cv.split(mat_filter_rescale)


    # save the results for checking
    # read metadata info
    if verbose is True:
        metadata = gdal.Open(uri_composite_gdal)
        trans = metadata.GetGeoTransform()
        # proj = metadata.GetProjection()
        out_path = os.path.join(working_dir, 'tile{}_{}_meanshift.tif'.format(tile_id, season))
        gdal_save_file_tif_3bands(out_path, b2_filter.reshape(cols, rows),
                                  b3_filter.reshape(cols, rows),
                                  b4_filter.reshape(cols, rows),
                                  gdal.GDT_Float32, trans, proj, rows, cols)

    # STEP 2
    # sober filtering + watershed
    r_sobel = sobel(b2_filter.reshape(cols, rows))
    g_sobel = sobel(b3_filter.reshape(cols, rows))
    b_sobel = sobel(b4_filter.reshape(cols, rows))
    gradient = r_sobel + g_sobel + b_sobel
    gradient_subset = gradient[buf:cols - buf, buf:rows - buf]
    # peak_gradient = peak_local_max(-gradient, num_peaks=2400, min_distance=10, indices=False)
    # markers = ndi.label(peak_gradient)[0]

    # For a 2D image, a connectivity of 1 corresponds to immediate neighbors up, down, left, and right,
    # while a connectivity of 2 also includes diagonal neighbors.
    segments_watershed = watershed(gradient_subset, markers=2400, connectivity=2, compactness=0).astype(np.int16)

    # read metadata info
    metadata = gdal.Open(uri_prob_gdal)
    trans = metadata.GetGeoTransform()
    # proj = metadata.GetProjection()

    if verbose is True:
        out_path = os.path.join(working_dir, 'tile{}_{}_meanshift_sober.tif'.format(tile_id, season))
        gdal_save_file_tif_1bands(out_path, gradient_subset, gdal.GDT_Float32, trans, proj, rows - 2 * buf,
                                  cols - 2 * buf)
        out_path = os.path.join(working_dir, 'tile{}_{}_meanshift_watershed.tif'.format(tile_id, season))
        gdal_save_file_tif_1bands(out_path, segments_watershed, gdal.GDT_Int16, trans, proj, rows - 2 * buf,
                                  cols - 2 * buf)

    # segments_slic = slic(arr_filter_rescale, n_segments=2400,compactness=0.1)
    # out_path = os.path.join(working_path, 'split_tile1_os_slic.tif')
    # gdal_save_file_tif_1bands(out_path, segments_slic, gdal.GDT_Int16, trans, proj, rows, cols)

    # STEP 3
    # overlapped with prob image, and selected those polygons that have average probability over 0.5
    # segments_watershed_merge_sieve = gdal_array.LoadFile(os.path.join(working_dir, tile_id + '_watershed_merge_sieve.tif'))
    # this way is much faster than overlapping with vectorized polygons
    for i in range(1, np.max(segments_watershed)):
        condition = np.equal(segments_watershed, i)
        prob_condition = np.extract(condition, array_prob)
        if np.mean(prob_condition) < prob_threshold:
            segments_watershed[condition] = 0

    if verbose is True:
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap.tif'.format(tile_id, season))
        gdal_save_file_tif_1bands(out_path, segments_watershed, gdal.GDT_Int16, trans, proj, rows - 2 * buf,
                                  cols - 2 * buf)

    # STEP 4
    # connected small polygons using graph theory + sieving filter
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    array_original_rescale = np.dstack((min_max_scaler.fit_transform(b2.reshape(cols * rows, 1).astype(np.float32))
                                        .reshape(cols, rows),
                                        min_max_scaler.fit_transform(b3.reshape(cols * rows, 1).astype(np.float32))
                                        .reshape(cols, rows),
                                        min_max_scaler.fit_transform(b4.reshape(cols * rows, 1).astype(np.float32))
                                        .reshape(cols, rows)))

    b2_filter_withingrid = array_original_rescale[buf:cols - buf, buf:rows - buf, 0]
    field_b2_std = np.std(b2_filter_withingrid[array_prob > prob_threshold])
    # field_b2_std = np.std(b2_filter_withingrid[arr_prob > prob_threshold])
    # nofield_b2_mean = np.mean(b2_filter_withingrid[array_prob <= prob_threshold])
    # nofield_b2_std = np.std(b2_filter_withingrid[arr_prob <= prob_threshold])

    b3_filter_withingrid = array_original_rescale[buf:cols - buf, buf:rows - buf, 1]
    field_b3_std = np.std(b3_filter_withingrid[array_prob > prob_threshold])
    # nofield_b3_mean = np.mean(b3_filter_withingrid[array_prob <= threshold])

    b4_filter_withingrid = array_original_rescale[buf:cols - buf, buf:rows - buf, 2]
    field_b4_std = np.std(b4_filter_withingrid[array_prob > prob_threshold])
    # nofield_b4_mean = np.mean(b4_filter_withingrid[array_prob <= prob_threshold])

    # calculate norm to adaptively decide merging threshold
    # diff_norm = np.linalg.norm([nofield_b2_mean - field_b2_mean, nofield_b3_mean - field_b3_mean,
    #                             nofield_b4_mean - field_b4_mean])
    diff_norm = np.linalg.norm([field_b2_std, field_b3_std, field_b4_std])

    array_original_subset = array_original_rescale[buf:cols - buf, buf:rows - buf]

    # assign -9999 so background pixels won't be merged
    array_original_subset[segments_watershed == 0] = [-9999, -9999, -9999]

    g = graph.rag_mean_color(array_original_subset, segments_watershed,
                             mode='distance')

    segments_watershed_merge = graph.merge_hierarchical(segments_watershed, g, thresh=diff_norm / 2, rag_copy=False,
                                                        in_place_merge=True,
                                                        merge_func=merge_mean_color,
                                                        weight_func=weight_mean_color)

    # note: after merging, it sometimes reset polygons id, which cause the id of no field to be not zero any more,
    # the below is the function to fix this issue
    for i in range(np.max(segments_watershed_merge)):
        condition = np.equal(segments_watershed_merge, i)
        id_condition = np.extract(condition, segments_watershed)
        counts = np.bincount(id_condition)
        if np.argmax(counts) == 0:
            nofield_id = i
            break

    if nofield_id is not 0:
        # for nofield id is not 0, switch 0 and nofield id
        segments_watershed_merge[segments_watershed_merge == 0] = 9999
        segments_watershed_merge[segments_watershed_merge == nofield_id] = 0
        segments_watershed_merge[segments_watershed_merge == 9999] = nofield_id

    # required to save it for sieve filter
    out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_merge.tif'.format(tile_id, season))
    gdal_save_file_tif_1bands(out_path, segments_watershed_merge, gdal.GDT_Int16, trans, proj, rows - 2 * buf,
                              cols - 2 * buf)

    # sieving filter
    cmd = 'gdal_sieve.py -q -st {} -8 {} {}'.format(mmu, os.path.join(working_dir,
                                                                      'tile{}_{}_watershed_overlap_merge.tif'.format(
                                                                          tile_id, season)),
                                                    os.path.join(working_dir,
                                                                 'tile{}_{}_watershed_overlap_merge_sieve.tif'.format(
                                                                     tile_id, season)))
    run_cmd(cmd, logger)


    # STEP 5
    # polyganize
    out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_merge_sieve.tif'.format(tile_id, season))

    # check if it is valid image before polygonization
    if not is_valid_image(out_path):
        logger.error("Segmentation failed for {} at {} season".format(tile_id, season))

    src_ds = gdal.Open(out_path)

    srcband = src_ds.GetRasterBand(1)

    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)
    dst_layername = os.path.join(working_dir, 'tile{}_{}_watershed_merge_sieve_overlap'.format(tile_id, season))
    drv = ogr.GetDriverByName("GeoJSON")
    dst_ds = drv.CreateDataSource(dst_layername + ".geojson")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)

    fieldName = "id"
    fd = ogr.FieldDefn(fieldName, ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("id")

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
    dst_ds = None  # it guaranttee shapefile is successfully created
    src_ds = None

    # STEP 6
    # post-processing
    infile = dst_layername + ".geojson"
    outfile = os.path.join(working_dir, 'tile{}_{}_seg.geojson'.format(tile_id, season))
    command = '/usr/bin/Rscript'
    path_script = "Postprocessing.R"
    args = [infile, outfile, '0.00005']
    if os.path.isfile(path_script) is False:
        logger.error("Fail to find Postprocessing.R")

    # check_output will run the command and store to result
    cmd = [command, path_script] + args
    x = subprocess.call(cmd, stdout=open(os.devnull, 'wb'))

    # remove temporal file
    if verbose is False:
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_merge.tif'.format(tile_id, season))
        os.remove(out_path)
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_merge_sieve.tif'.format(tile_id, season))
        os.remove(out_path)
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_merge_sieve_overlap.geojson'.format(tile_id, season))
        os.remove(out_path)


def segmentation_execution_doubleseasons(s3_bucket, planet_directory, prob_directory, tile_id, dry_lower_ordinal,
                                         dry_upper_ordinal, wet_lower_ordinal, wet_upper_ordinal, tile_col, tile_row,
                                         working_dir, mmu, prob_threshold, buf, output_s3_prefix, logger, verbose):
    """
    s3_bucket: s3 bucket name
    planet_directory: s3 directory for planet composite
    prob_directory: s3 directory for probability image
    tile_id: tile id
    dry_lower_ordinal: lower bound (ordinal date) for dry season
    dry_upper_ordinal: upper bound (ordinal date) for dry season
    wet_lower_ordinal: lower bound (ordinal date) for wet season
    wet_upper_ordinal: lower bound (ordinal date) for wet season
    tile_col: column of tiles in our coarse layout
    tile_row: row of tiles in our coarse layout
    working_dir: outputted directory
    mmu: minimum mapping units
    prob_threshold:  probability threshold
    buf: buffer of composite image
    output_s3_prefix: outputed prefix for segmentation result in s3
    logger: logger object
    verbose: True or False, if outputted intermediate results
    return:
    True or False
    """

    uri_composite_gdal_os = "/vsis3/{}/{}/{}/tile{}_{}_{}.tif".format(s3_bucket, planet_directory, 'OS', tile_id,
                                                                       dry_lower_ordinal, dry_upper_ordinal)
    uri_composite_gdal_gs = "/vsis3/{}/{}/{}/tile{}_{}_{}.tif".format(s3_bucket, planet_directory, 'GS', tile_id,
                                                                       wet_lower_ordinal, wet_upper_ordinal)
    uri_prob_gdal = "/vsis3/{}/{}/image_c{}_r{}.tif".format(s3_bucket, prob_directory, str(tile_col), str(tile_row))

    # segmentation for off-season
    segmentation_season(tile_id, 'OS', uri_composite_gdal_os, uri_prob_gdal, working_dir, mmu, prob_threshold, buf,
                        logger, verbose)

    # segmentation for growing-season
    segmentation_season(tile_id, 'GS', uri_composite_gdal_gs, uri_prob_gdal, working_dir, mmu, prob_threshold, buf,
                        logger, verbose)

    ############################################################
    #             upload compositing image to s3             #
    ############################################################
    s3 = boto3.client('s3')
    out_path_dry = os.path.join(working_dir, 'tile{}_OS_seg.geojson'.format(tile_id))
    try:
        s3.upload_file(out_path_dry, s3_bucket, '{}/OS/tile{}_{}_{}_seg.geojson'.format(output_s3_prefix, tile_id,
                                                                                 dry_lower_ordinal, dry_upper_ordinal))
    except ClientError as e:
        logger.error("S3 uploading fails for tile{}_{}_{}_seg.geojson: {}".format(tile_id, dry_lower_ordinal, dry_upper_ordinal, e))
        raise

    out_path_wet = os.path.join(working_dir, 'tile{}_GS_seg.geojson'.format(tile_id))

    try:
        s3.upload_file(out_path_wet, s3_bucket, '{}/GS/tile{}_{}_{}_seg.geojson'.format(output_s3_prefix,  tile_id, wet_lower_ordinal,
                                                                                wet_upper_ordinal))
    except ClientError as e:
        logger.error("S3 uploading fails for tile{}_{}_{}_seg.geojson: {}".format(tile_id, wet_lower_ordinal, wet_upper_ordinal, e))
        raise

    if os.path.exists(out_path_dry):
        os.remove(out_path_dry)
    else:
        logger.error("Segmentation fails: couldn't find {}").format(out_path_dry)

    if os.path.exists(out_path_wet):
        os.remove(out_path_wet)
    else:
        logger.error("Segmentation fails: couldn't find {}").format(out_path_wet)

@click.command()
@click.option('--config_filename', default='segmenter_config.yaml', help='The name of the config to use.')
@click.option('--tile_id', default=None, help='only used for debug mode, user-defined tile_id')
@click.option('--csv_pth', default=None, help='csv path for providing a specified tile list')
@click.option('--aoi', default=None, help='specify production AOI id in ghana_tiles.geojson')
@click.option('--s3_bucket', default='activemapper', help='s3 bucket name')
@click.option('--threads_number', default= 4, help='output folder prefix')
@click.option('--verbose', default= False , help='output folder prefix')


def main(config_filename, tile_id, csv_pth, aoi, s3_bucket, threads_number, verbose):

    buf = 11 # buffer of composite image
    left_corner_x = -17.541
    left_corner_y = 37.54
    per_tile_width = 0.005 * 10  # 0.005 degree is the width of cells, 1 tile has 10*10 cells
    working_dir = '/tmp'
    verbose = False

    log_path = '%s/log/segmenter_%s.log' % (os.environ['HOME'], str(aoi))
    logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
    logger = logging.getLogger(__name__)

    # read yaml from local
    with open("/home/ubuntu/source/segmenter_config.yaml", 'r') as yaml_obj:
        params = yaml.safe_load(yaml_obj)['segmenter']

    # read yaml from s3
    # params = parse_yaml_from_s3(s3_bucket, config_filename)['mapper']

    if params is None:
        logger.error("Failed to open yaml file")

    prefix = params['planet_prefix']
    planet_directory = params['planet_directory']
    prob_directory = params['prob_directory']
    tiles_geojson_path = params['tile_geojson_path']
    mmu = params['mmu']
    prob_threshold = params['prob_threshold']
    output_s3_prefix = params['output']
    dry_lower_ordinal = params['dry_lower_ordinal']  # 2018/12/01
    dry_upper_ordinal = params['dry_upper_ordinal']  # 2019/02/28
    wet_lower_ordinal = params['wet_lower_ordinal']  # 2018/05/01
    wet_upper_ordinal = params['wet_upper_ordinal']  # 2018/09/30


    uri_tile = "s3://{}/{}/{}".format(s3_bucket, prefix, tiles_geojson_path)
    gpd_tile = gpd.read_file(uri_tile)
    if gpd_tile is None:
        logger.error("reading geojson tile '{}' failed".format(uri_tile))

    if tile_id is None: # for aoi or list mode, parallel should be applied
        # determine thread number to be used
        if threads_number == 'default':
            threads_number = multiprocessing.cpu_count() * 2
        else:
            threads_number = int(threads_number)

        segmentation_composition_executor = FixedThreadPoolExecutor(size=threads_number)
        # mode 1: tile-csv based
        if aoi is None:
            log_path = '%s/log/segmentation.log' % os.environ['HOME']
            logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
            logger = logging.getLogger(__name__)

            # time zone
            tz = timezone('US/Eastern')
            logger.info(
                "Progress: starting a segmentation task ({})".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

            # read a geopandas object for tile geojson
            uri_tile = "s3://{}/{}/{}".format(s3_bucket, prefix, tiles_geojson_path)
            gpd_tile = gpd.read_file(uri_tile)
            if gpd_tile is None:
                logger.error("reading geojson tile '{}' failed".format(uri_tile))

            if csv_pth is None:
                logger.error("Please provide tile_id, csv_path or aoi_id ({}))"
                             .format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
                return

            # tile-list processing
            alltiles = pd.read_csv(csv_pth)['tile']
        # mode 2: aoi based (production-based)
        else:
            # define log path
            log_path = '%s/log/segmentation_aoi%s.log' % (os.environ['HOME'], aoi)
            logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
            logger = logging.getLogger(__name__)

            # time zone
            tz = timezone('US/Eastern')
            logger.info("Progress: starting a segmentation task ({})".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

            # read a geopandas object for tile geojson
            uri_tile = "s3://{}/{}/{}".format(s3_bucket, prefix, tiles_geojson_path)
            gpd_tile = gpd.read_file(uri_tile)
            if gpd_tile is None:
                logger.error("reading geojson tile '{}' failed".format(uri_tile))

            # here we used merged aoi
            alltiles = gpd_tile.loc[gpd_tile['merged_aoi'] == float(aoi)]['tile']
            # alltiles = gpd_tile.loc[gpd_tile['aoi'] == float(aoi)]['tile']

        failure_count = 0
        success_count = 0
        # looping over each tile
        for i in range(len(alltiles)):

            # retrieve all tile info for focused tile_id
            tile_id = int(alltiles.iloc[i])

            # calculate global column and row for this tile based on center coordinates of tiles
            foc_gpd_tile = gpd_tile[gpd_tile['tile'] == int(tile_id)]
            (tile_col, tile_row) = get_colrow_geojson(foc_gpd_tile, left_corner_x, left_corner_y, per_tile_width,
                                                      logger)

            # implement segmentation for both seasons
            segmentation_composition_executor.submit(segmentation_execution_doubleseasons, s3_bucket, planet_directory,
                                                     prob_directory, tile_id, dry_lower_ordinal, dry_upper_ordinal,
                                                     wet_lower_ordinal, wet_upper_ordinal, tile_col, tile_row,
                                                     working_dir, mmu, prob_threshold, buf, output_s3_prefix, logger, verbose)

        # await all tile finished
        segmentation_composition_executor.drain()

        # await threadpool to stop
        segmentation_composition_executor.close()

        logger.info("Progress: finished segmentation task for aoi {}; the total tile number to be processed is {}; "
                    "the success_count is {}; the failure_count is {} ({})"
                    .format(aoi, len(alltiles), success_count, failure_count, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

    else:
        # define log path
        log_path = '%s/log/segmentation_tile%s.log' % (os.environ['HOME'], str(tile_id))
        logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
        logger = logging.getLogger(__name__)

        # time zone
        tz = timezone('US/Eastern')
        logger.info("Progress: starting a segmentation task ({})".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


        # read a geopandas object for tile geojson
        uri_tile = "s3://{}/{}/{}".format(s3_bucket, prefix, tiles_geojson_path)
        gpd_tile = gpd.read_file(uri_tile)
        if gpd_tile is None:
            logger.error("reading geojson tile '{}' failed". format(uri_tile))

        # calculate global column and row for this tile based on center coordinates of tiles
        foc_gpd_tile = gpd_tile[gpd_tile['tile'] == int(tile_id)]
        (tile_col, tile_row) = get_colrow_geojson(foc_gpd_tile, left_corner_x, left_corner_y, per_tile_width)

        try:
            segmentation_execution_doubleseasons(s3_bucket, planet_directory, prob_directory, tile_id, dry_lower_ordinal,
                                             dry_upper_ordinal, wet_lower_ordinal, wet_upper_ordinal, tile_col, tile_row,
                                             working_dir, mmu, prob_threshold, buf, output_s3_prefix, logger, verbose)
        except (OSError, ClientError, subprocess.CalledProcessError) as e:
            logger.error("Segmentation failed for tile_id {}, and the total finished tiles is {} ({}))"
                         .format(tile_id,  1, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        else:
            logger.info("Progress: finished segmentation for tile_id {}, and the total finished tiles is {} ({}))"
                        .format(tile_id,  1, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


if __name__ == '__main__':
    main()
