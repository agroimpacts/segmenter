# This script is to implement segmentation based on probability images and original composite images
# The function 'segmentation_season' is the main function for segmentation steps
# The whole algorithm needs two inputs: probability image (which you can obtain from any classifer) and original RS
# images, and the output will be a vector map
# The processing is consisted of 6 steps: 1) mean-shift smoothing; 2) watershed segmentation (needs
# over-segment); 3) remain polygons that the average probability is over a threshold; 4) sieving filter + hierarchical
# merging; 5) polygonization; 6) check validity of polygons and simplify boundary
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
from skimage.future.graph import RAG
import numpy as np
import math
import heapq
from skimage.measure import perimeter
import sys


def _revalidate_node_edges(rag, node, heap_list):
    """Handles validation and invalidation of edges incident to a node.
    This function invalidates all existing edges incident on `node` and inserts
    new items in `heap_list` updated with the valid weights.
    rag : RAG
        The Region Adjacency Graph.
    node : int
        The id of the node whose incident edges are to be validated/invalidated
        .
    heap_list : list
        The list containing the existing heap of edges.
    """
    # networkx updates data dictionary if edge exists
    # this would mean we have to reposition these edges in
    # heap if their weight is updated.
    # instead we invalidate them

    for nbr in rag.neighbors(node):
        data = rag[node][nbr]
        try:
            # invalidate edges incident on `dst`, they have new weights
            data['heap item'][3] = False
            _invalidate_edge(rag, node, nbr)
        except KeyError:
            # will handle the case where the edge did not exist in the existing
            # graph
            pass

        wt = data['weight']
        heap_item = [wt, node, nbr, True]
        data['heap item'] = heap_item
        heapq.heappush(heap_list, heap_item)


def _rename_node(graph, node_id, copy_id):
    """ Rename `node_id` in `graph` to `copy_id`. """

    graph._add_node_silent(copy_id)
    graph.nodes[copy_id].update(graph.nodes[node_id])

    for nbr in graph.neighbors(node_id):
        wt = graph[node_id][nbr]['weight']
        graph.add_edge(nbr, copy_id, {'weight': wt})

    graph.remove_node(node_id)


def _invalidate_edge(graph, n1, n2):
    """ Invalidates the edge (n1, n2) in the heap. """
    graph[n1][n2]['heap item'][3] = False


def merge_hierarchical_customized(labels, rag, thresh, rag_copy, in_place_merge, _maximum_fieldsize, _shape_threshold,
                                  _include_texture, merge_func, weight_func, weight_func_texture):
    """Perform hierarchical merging of a RAG.
    Greedily merges the most similar pair of nodes until no edges lower than
    `thresh` remain.
    Parameters
    ----------
    labels : ndarray
        The array of labels.
    rag : RAG
        The Region Adjacency Graph.
    thresh : float
        Regions connected by an edge with weight smaller than `thresh` are
        merged.
    rag_copy : bool
        If set, the RAG copied before modifying.
    in_place_merge : bool
        If set, the nodes are merged in place. Otherwise, a new node is
        created for each merge..
    _maximum_fieldsize: int
        The maximum field size. When reached this threshold, the merging will stop
    _shape_threshold: double
        The maximum shape index for regularity measurement for polygons
    _include_texture: bool
        Indicate if texture info is included for merging
    merge_func : callable
        This function is called before merging two nodes. For the RAG `graph`
        while merging `src` and `dst`, it is called as follows
        ``merge_func(graph, src, dst)``.
    weight_func : callable
        The function to compute the new weights of the nodes adjacent to the
        merged node. This is directly supplied as the argument `weight_func`
        to `merge_nodes`.
    weight_func_texture:
        The function to compute the weights by incoporating both color and texture (variance)
    Returns
    -------
    out : ndarray
        The new labeled array.
    """
    if rag_copy:
        rag = rag.copy()

    edge_heap = []
    for n1, n2, data in rag.edges(data=True):
        # Push a valid edge in the heap
        wt = data['weight']
        heap_item = [wt, n1, n2, True]
        heapq.heappush(edge_heap, heap_item)

        # Reference to the heap item in the graph
        data['heap item'] = heap_item

    while len(edge_heap) > 0 and edge_heap[0][0] < thresh:
        _, n1, n2, valid = heapq.heappop(edge_heap)
        # Ensure popped edge is valid, if not, the edge is discarded
        if valid:
            # test 1: size test
            be_mergy = True
            if rag.nodes[n1]['pixel count'] + rag.nodes[n2]['pixel count'] > _maximum_fieldsize:
                be_mergy = False
            else:
                if rag.nodes[n1]['pixel count'] is not 0 and rag.nodes[n2]['pixel count'] is not 0:
                    shape_idx_1 = 0.25 * perimeter(np.isin(labels, rag.nodes(data=True)[n1]['labels']), neighbourhood=4) \
                                  / np.sqrt(rag.nodes[n1]['pixel count'])
                    if shape_idx_1 < _shape_threshold:
                        shape_idx_2 = 0.25 * perimeter(np.isin(labels, rag.nodes(data=True)[n2]['labels']), neighbourhood=4) \
                                      / np.sqrt(rag.nodes[n2]['pixel count'])
                        if shape_idx_2 < _shape_threshold:
                            be_mergy =False

            if be_mergy is True:
                # Invalidate all neigbors of `src` before its deleted

                for nbr in rag.neighbors(n1):
                    _invalidate_edge(rag, n1, nbr)

                for nbr in rag.neighbors(n2):
                    _invalidate_edge(rag, n2, nbr)

                if not in_place_merge:
                    next_id = rag.next_id()
                    _rename_node(rag, n2, next_id) # next_id return the id for the new mode to be inserted
                    src, dst = n1, next_id
                else:
                    src, dst = n1, n2

                merge_func(rag, src, dst)

                if _include_texture is True:
                    new_id = rag.merge_nodes(src, dst, weight_func_texture)
                else:
                    new_id = rag.merge_nodes(src, dst, weight_func)

                _revalidate_node_edges(rag, new_id, edge_heap)

    label_map = np.arange(labels.max() + 1)
    for ix, (n, d) in enumerate(rag.nodes(data=True)):
        for label in d['labels']:
            label_map[label] = ix

    return label_map[labels]


def rag_mean_variance_color(image, labels, _include_texture, connectivity=2, mode='distance',
                   sigma=255.0):
    """Compute the Region Adjacency Graph using mean colors.
    Given an image and its initial segmentation, this method constructs the
    corresponding Region Adjacency Graph (RAG). Each node in the RAG
    represents a set of pixels within `image` with the same label in `labels`.
    The weight between two adjacent regions represents how similar or
    dissimilar two regions are depending on the `mode` parameter.
    Parameters
    ----------
    image : ndarray, shape(M, N, [..., P,] 3)
        Input image.
    labels : ndarray, shape(M, N, [..., P])
        The labelled image. This should have one dimension less than
        `image`. If `image` has dimensions `(M, N, 3)` `labels` should have
        dimensions `(M, N)`.
    _include_texture: bool
        Indicate if texture info is included for weight calculation
    connectivity : int, optional
        Pixels with a squared distance less than `connectivity` from each other
        are considered adjacent. It can range from 1 to `labels.ndim`. Its
        behavior is the same as `connectivity` parameter in
        ``scipy.ndimage.generate_binary_structure``.
    mode : {'distance', 'similarity'}, optional
        The strategy to assign edge weights.
            'distance' : The weight between two adjacent regions is the
            :math:`|c_1 - c_2|`, where :math:`c_1` and :math:`c_2` are the mean
            colors of the two regions. It represents the Euclidean distance in
            their average color.
            'similarity' : The weight between two adjacent is
            :math:`e^{-d^2/sigma}` where :math:`d=|c_1 - c_2|`, where
            :math:`c_1` and :math:`c_2` are the mean colors of the two regions.
            It represents how similar two regions are.
    sigma : float, optional
        Used for computation when `mode` is "similarity". It governs how
        close to each other two colors should be, for their corresponding edge
        weight to be significant. A very large value of `sigma` could make
        any two colors behave as though they were similar.
    Returns
    -------
    out : RAG
        The region adjacency graph.
    References
    ----------
    .. [1] Alain Tremeau and Philippe Colantoni
           "Regions Adjacency Graph Applied To Color Image Segmentation"
           :DOI:`10.1109/83.841950`
    """
    graph = RAG(labels, connectivity=connectivity)

    for n in graph:
        graph.nodes[n].update({'labels': [n],
                               'pixel count': 0,
                               'total color': np.array([0, 0, 0, 0], dtype=np.double),
                               'total color square': np.array([0, 0, 0, 0], dtype=np.double),
                               'variance color': np.array([0, 0, 0, 0], dtype=np.double),
                               'shape index': 0})

    for index in np.ndindex(labels.shape):
        current = labels[index]
        graph.nodes[current]['pixel count'] += 1
        graph.nodes[current]['total color'] += image[index]
        graph.nodes[current]['total color square'] += image[index] * image[index]

    for n in graph:
        graph.nodes[n]['mean color'] = (graph.nodes[n]['total color'] /
                                        graph.nodes[n]['pixel count'])
        graph.nodes[n]['variance color'] = (graph.nodes[n]['total color square'] / graph.nodes[n]['pixel count']) - \
                                           (graph.nodes[n]['mean color'] * graph.nodes[n]['mean color'])

    for x, y, d in graph.edges(data=True):
        if _include_texture is True:
            diff = np.absolute(graph.nodes[x]['mean color'] - graph.nodes[y]['mean color']) * 0.8 + \
                   np.sqrt(np.absolute(graph.nodes[x]['variance color'] - graph.nodes[y]['variance color'])) * 0.2
        else:
            diff = graph.node[x]['mean color'] - graph.node[y]['mean color']
        diff = np.linalg.norm(diff)
        if mode == 'similarity':
            d['weight'] = math.e ** (-(diff ** 2) / sigma)
        elif mode == 'distance':
            d['weight'] = diff
        else:
            raise ValueError("The mode '%s' is not recognised" % mode)

    return graph


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
    diff = graph.nodes[dst]['mean color'] - graph.nodes[n]['mean color']
    diff = np.linalg.norm(diff)
    return {'weight': diff}


def weight_mean_color_texture(graph, src, dst, n):
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

    diff = np.absolute(graph.nodes[dst]['mean color'] - graph.nodes[n]['mean color']) * 0.8 + \
            np.sqrt(np.absolute(graph.nodes[dst]['variance color'] - graph.nodes[n]['variance color'])) * 0.2
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
    graph.nodes[dst]['total color'] += graph.nodes[src]['total color']
    graph.nodes[dst]['pixel count'] += graph.nodes[src]['pixel count']
    graph.nodes[dst]['total color square'] += graph.nodes[src]['total color square']
    graph.nodes[dst]['mean color'] = (graph.nodes[dst]['total color'] /
                                     graph.nodes[dst]['pixel count'])
    graph.nodes[dst]['variance color'] = (graph.nodes[dst]['total color square'] / graph.nodes[dst]['pixel count']) - \
                                        (graph.nodes[dst]['mean color'] * graph.nodes[dst]['mean color'])


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


def segmentation_season(tile_id, season, uri_composite_gdal, uri_prob_gdal, working_dir, mmu, maximum_field_size,
                        prob_threshold, buf, logger, verbose):
    """
    segmentation for signle season
    tile_id: tile id
    season: season
    uri_composite_gdal: s3 url for composite image
    uri_prob_gdal: s3 url for probability image
    working_dir: outputted folder
    mmu: minimum mapping unit (unit: pixel)
    maximum_field_size: maximum field size (unit: pixel)
    maximum_field_size:
    prob_threshold: threshold for probability
    buf: buffer of composite image (the default is 11 pixels)
    logger: logger
    verbose: True or False, if outputted intermediate results
    return: True or False
    """

    proj = ''
    shape_threshold = 1.15
    include_texture = False

    # a quick way to read image as numpy array
    array_composite = gdal_array.LoadFile(uri_composite_gdal)
    array_prob = gdal_array.LoadFile(uri_prob_gdal)

    # get channels, cols, rows
    [nchannels, cols, rows] = array_composite.shape

    # temporal change for dealing with Ron's Rf probability image
    # array_prob = array_prob[buf:cols - buf, buf:rows - buf]

    # STEP 1
    # meanshift algorithm to smooth image and filter out noise
    b1, b2, b3, b4 = array_composite

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
    segments_watershed = watershed(gradient_subset, markers=3600, connectivity=1, compactness=0).astype(np.int16)

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

    # required to save image, as needed for gdal sieving function
    out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap.tif'.format(tile_id, season))
    gdal_save_file_tif_1bands(out_path, segments_watershed, gdal.GDT_Int16, trans, proj, rows - 2 * buf,
                              cols - 2 * buf)


    # STEP 4
    # connected small polygons using sieving filter + hierachical merging

    # first apply sieving filter, as sieving substantially reduce the node number in graph, thereby improve efficiency
    cmd = 'gdal_sieve.py -q -st {} -8 {} {}'.format(int(mmu), os.path.join(working_dir, 'tile{}_{}_watershed_overlap.tif'
                                                                      .format(tile_id, season)),
                                                         os.path.join(working_dir, 'tile{}_{}_watershed_overlap_sieve.tif'
                                                                      .format(tile_id, season)))
    os.system(cmd)
    # reopen sieved image
    segments_watershed_sieve = gdal_array.LoadFile(os.path.join(working_dir, 'tile{}_{}_watershed_overlap_sieve.tif'
                                                                .format(tile_id, season)))

    # then hierachical merging
    # min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
    # array_original_stack = np.dstack((b1, b2, b3, b4))
    array_original_subset = np.dstack((b1, b2, b3, b4))[buf:cols - buf, buf:rows - buf]

    # calculate adaptive threshold based stand deviation for 'field' pixels
    std_band = np.std(array_original_subset[segments_watershed_sieve > 0], axis=0)
    vec_std = np.linalg.norm(std_band)

    # assign -9999 so background pixels won't be merged
    array_original_subset[segments_watershed_sieve == 0] = [-9999, -9999, -9999, -9999]

    g = rag_mean_variance_color(array_original_subset, segments_watershed_sieve, connectivity=1,
                                mode='distance', _include_texture=include_texture)

    segments_watershed_merge = merge_hierarchical_customized(segments_watershed_sieve, g, thresh=vec_std/3,
                                                             rag_copy=False,
                                                             in_place_merge=False,
                                                             _maximum_fieldsize=maximum_field_size,
                                                             _shape_threshold=shape_threshold,
                                                             _include_texture=include_texture,
                                                             merge_func=merge_mean_color,
                                                             weight_func=weight_mean_color,
                                                             weight_func_texture=weight_mean_color_texture)

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
    
    # required to save the merged polygon cuz polygonization needs
    out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_sieve_merge.tif'.format(tile_id, season))
    gdal_save_file_tif_1bands(out_path, segments_watershed_merge, gdal.GDT_Int16, trans, proj, rows - 2 * buf,
                              cols - 2 * buf)

    # STEP 5
    # polyganize

    # check if it is valid image before polygonization
    if not is_valid_image(out_path):
        logger.error("Segmentation failed for {} at {} season".format(tile_id, season))

    src_ds = gdal.Open(out_path)

    srcband = src_ds.GetRasterBand(1)

    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj)
    dst_layername = os.path.join(working_dir, 'tile{}_{}_watershed_sieve_merge_overlap'.format(tile_id, season))
    drv = ogr.GetDriverByName("GeoJSON")
    dst_ds = drv.CreateDataSource(dst_layername + ".geojson")
    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)

    fieldName = "id"
    fd = ogr.FieldDefn(fieldName, ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("id")

    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
    
    # release memory
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
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap.tif'.format(tile_id, season))
        os.remove(out_path)
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_sieve.tif'.format(tile_id, season))
        os.remove(out_path)
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_overlap_sieve_merge.tif'.format(tile_id, season))
        os.remove(out_path)
        out_path = os.path.join(working_dir, 'tile{}_{}_watershed_sieve_merge_overlap.geojson'.format(tile_id, season))
        os.remove(out_path)


def segmentation_execution_doubleseasons(s3_bucket, planet_directory, prob_directory, tile_id, dry_lower_ordinal,
                                         dry_upper_ordinal, wet_lower_ordinal, wet_upper_ordinal, tile_col, tile_row,
                                         working_dir, mmu, maximum_field_size, prob_threshold,
                                         buf, output_s3_prefix, logger, verbose):
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
    mmu: minimum mapping units (unit: pixel)
    maximum_field_size: maximum field size (unit: pixel)
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
    uri_prob_gdal = "/vsis3/{}/{}/image_c{}_r{}_8_run0_iteration4.tif".format(s3_bucket, prob_directory, str(tile_col), str(tile_row))
    #uri_prob_gdal = "/vsis3/{}/{}/image_c{}_r{}.tif".format(s3_bucket, prob_directory, str(tile_col), str(tile_row))
    # segmentation for off-season
    segmentation_season(tile_id, 'OS', uri_composite_gdal_os, uri_prob_gdal, working_dir, mmu, maximum_field_size,
                        prob_threshold, buf, logger, verbose)

    # segmentation for growing-season
    segmentation_season(tile_id, 'GS', uri_composite_gdal_gs, uri_prob_gdal, working_dir, mmu, maximum_field_size,
                        prob_threshold, buf, logger, verbose)

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
@click.option('--config_filename', type=str, default='segmenter_config.yaml', help='The name of the config to use.')
@click.option('--tile_id', type=int, default=None, help='only used for debug mode, user-defined tile_id')
@click.option('--csv_pth', type=str, default=None, help='csv path for providing a specified tile list')
@click.option('--aoi', type=int, default=None, help='specify production AOI id in ghana_tiles.geojson')
@click.option('--s3_bucket', type=str, default='activemapper', help='s3 bucket name')
@click.option('--threads_number', type=int,default= 4, help='output folder prefix')
@click.option('--be_minmax_analysis', is_flag=True, help='if extract min and max from worker labels')
@click.option('--verbose', is_flag=True, help='if output folder prefix')
def main(config_filename, tile_id, csv_pth, aoi, s3_bucket, threads_number, be_minmax_analysis, verbose):
    # some constants are defined here
    buf = 11 # buffer of composite image
    left_corner_x = -17.541
    left_corner_y = 37.54
    per_tile_width = 0.005 * 10  # 0.005 degree is the width of cells, 1 tile has 10*10 cells
    working_dir = '/tmp'
    
    log_path = '%s/log/segmenter_%s.log' % (os.environ['HOME'], str(aoi))
    logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
    logger = logging.getLogger(__name__)

    # print(be_minmax_analysis)
    if be_minmax_analysis is True:
        command = '/usr/bin/Rscript'
        path_script = "Preprocessing.R"
        args = ['1']
        if os.path.isfile(path_script) is False:
            logger.error("Fail to find Preprocessing.R")
            sys.exit()

        # check_output will run the command and store to result
        cmd = [command, path_script] + args
        print(cmd)
        try:
            x = subprocess.call(cmd, stdout=open(os.devnull, 'wb'))
        except:
            logger.error("MinMax Analysis failed!")
            sys.exit()

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
    maximum_field_size = params['max_field_size']
    if mmu <= 0 or maximum_field_size <= 0:
        logger.error("Min and max polygon size were not set correctly!")
        sys.exit()
    else:
        print("min_max analysis finished!")

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

    # for aoi or list mode, parallel should be applied
    if tile_id is None:
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
            alltiles = gpd_tile.loc[gpd_tile['production_aoi'] == float(aoi)]['tile']
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
                                                     working_dir, mmu, maximum_field_size, prob_threshold, buf,
                                                     output_s3_prefix, logger, verbose)

        # await all tile finished
        segmentation_composition_executor.drain()

        # await threadpool to stop
        segmentation_composition_executor.close()

        logger.info("Progress: finished segmentation task for aoi {}; the total tile number to be processed is {}; "
                    "the success_count is {}; the failure_count is {} ({})"
                    .format(aoi, len(alltiles), success_count, failure_count, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
    # single tile processing
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
                                                 dry_upper_ordinal, wet_lower_ordinal, wet_upper_ordinal, tile_col,
                                                 tile_row, working_dir, mmu, maximum_field_size, prob_threshold, buf,
                                                 output_s3_prefix, logger, verbose)
        except (OSError, ClientError, subprocess.CalledProcessError) as e:
            logger.error("Segmentation failed for tile_id {}, and the total finished tiles is {} ({}))"
                         .format(tile_id,  1, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        else:
            logger.info("Progress: finished segmentation for tile_id {}, and the total finished tiles is {} ({}))"
                        .format(tile_id,  1, datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))


if __name__ == '__main__':
    main()
