"""Tests the luminance of PNG files against time-of-flight DAT files"""

import cv2
import numpy


def compare(first, second):
    if len(first) != len(second):
        print 'Comparison error: lists are different sizes'
        return
    diff = [first[x] - second[x] for x in range(0, len(first))]
    print '=== Comparing difference between two data sets ==='
    print 'Average Difference: ' + str(numpy.mean(diff))
    print 'Standard Deviation: ' + str(numpy.std(diff))
    print 'Average Percent Error: ' + str(abs(round(numpy.mean(diff) / numpy.mean(first) * 100, 2))) + "%"


def load(filename):
    extension = str(filename).rsplit(".", 1)[1]
    if extension == 'exr':
        return _loadexr(filename)
    if extension == 'png':
        return _loadpng(filename)
    if extension == 'dat':
        return _loaddat(filename)


def _loadexr(filename):
    img = cv2.imread(filename, cv2.CV_LOAD_IMAGE_UNCHANGED)
    if img is None:
        print 'Error loading image'
        return []
    rgb = [x for y in img for x in y]
    return [0.212671 * pixel[2] + 0.715160 * pixel[1] + 0.072169 * pixel[0] for pixel in rgb]


def _loadpng(filename):
    img = cv2.imread(filename, cv2.CV_LOAD_IMAGE_COLOR)
    if img is None:
        print 'Error loading image'
        return []
    rgb = [x for y in img for x in y]
    return [0.212671 * pixel[0] + 0.715160 * pixel[1] + 0.072169 * pixel[2] for pixel in rgb]


def _loaddat(filename):
    with open(filename, 'r') as f:
        data = f.read(50)
    n = len([x for x, y in enumerate(data.split(" ", 6)) if y == "#"])
    if n == 2:
        return _loadgroundtruth(filename)
    if n == 1:
        return _loadhistogram(filename)
    else:
        print 'Unrecognized DAT file format'
        return []


def _loadgroundtruth(filename):
    with open(filename, 'r') as f:
        data = f.read()
    return [float(y) for x, y in enumerate(data.split(" ")) if (x + 1) % 5 == 0]


def _loadhistogram(filename):
    with open(filename, 'r') as f:
        data = f.read()
    pixels = [x.strip() for x in data.split("#")[1:]]
    sums = [sum([float(y) for x, y in enumerate(pixel.split(" ")[2:]) if (x + 1) % 2 == 0]) for pixel in pixels]
    return sums
