# -*- coding: utf-8 -*-

import mapnik

try:
    import json
except ImportError:
    import simplejson as json


# returns true if pixels are not nearly identical
def compare_pixels(pixel1, pixel2, alpha=True, pixel_threshold=0):
    if pixel1 == pixel2:
        return False
    r_diff = abs((pixel1 & 0xff) - (pixel2 & 0xff))
    g_diff = abs(((pixel1 >> 8) & 0xff) - ((pixel2 >> 8) & 0xff))
    b_diff = abs(((pixel1 >> 16) & 0xff)- ((pixel2 >> 16) & 0xff))
    if alpha:
        a_diff = abs(((pixel1 >> 24) & 0xff) - ((pixel2 >> 24) & 0xff))
        if(r_diff > pixel_threshold or
           g_diff > pixel_threshold or
           b_diff > pixel_threshold or
           a_diff > pixel_threshold):
            return True
    else:
        if(r_diff > pixel_threshold or
           g_diff > pixel_threshold or
           b_diff > pixel_threshold):
            return True
    return False

# compare two images and return number of different pixels
def compare(actual, expected, alpha=True):
    im1 = mapnik.Image.open(actual)
    im2 = mapnik.Image.open(expected)
    pixels = im1.width() * im1.height()
    delta_pixels = (im2.width() * im2.height()) - pixels
    #diff = 0
    if delta_pixels != 0:
        return delta_pixels
    #for x in range(0,im1.width(),2):
    #    for y in range(0,im1.height(),2):
    #        if compare_pixels(im1.get_pixel(x,y),im2.get_pixel(x,y),alpha=alpha):
    #            diff += 1
    #return diff
    return im1.compare(im2, 0, alpha)

def compare_grids(actual, expected, threshold=0, alpha=True):
    global errors
    global passed
    im1 = json.loads(open(actual).read())
    im2 = json.loads(open(expected).read())
    # TODO - real diffing
    if not im1['data'] == im2['data']:
        return 99999999
    if not im1['keys'] == im2['keys']:
        return 99999999
    grid1 = im1['grid']
    grid2 = im2['grid']
    # dimensions must be exact
    width1 = len(grid1[0])
    width2 = len(grid2[0])
    if not width1 == width2:
        return 99999999
    height1 = len(grid1)
    height2 = len(grid2)
    if not height1 == height2:
        return 99999999
    diff = 0;
    for y in range(0,height1-1):
        row1 = grid1[y]
        row2 = grid2[y]
        width = min(len(row1),len(row2))
        for w in range(0,width):
            if row1[w] != row2[w]:
                diff += 1
    return diff
