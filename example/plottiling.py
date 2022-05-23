#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

import sys, argparse


def parse_argument():

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    parser = argparse.ArgumentParser()

    parser.add_argument('--tiling_plot', nargs=1, metavar="file", help='filename for the tiling plot')
    parser.add_argument('--tiling', nargs=1, metavar="file", help='file for the coordinates of the beams')
    parser.add_argument('--boresight', nargs=2, metavar=('RA', 'DEC'), help='boresight position in RA and DEC in h:m:s.s d:m:s.s')
    parser.add_argument('--beamshape', nargs=3, metavar=('x', 'y', 'deg'), help='semi-axis and orientation angle of the beam')
    parser.add_argument('--beam_size_scaling', nargs=1, metavar="scaling", help='scaling factor for the size of the beam')
    parser.add_argument("--flip", action="store_true", help='flip the orientation of the beam')
    parser.add_argument("--inner", nargs=1, metavar="number", help='highlight the most [number] inner beams')
    parser.add_argument("--extra_source", nargs=1, metavar="file", help='extra point sources to plot')

    args = parser.parse_args()

    return args

args = parse_argument()

coord_file = args.tiling[0]
if args.beam_size_scaling is not None:
    scaling = float(args.beam_size_scaling[0])
else:
    scaling = 1.0
coords = np.genfromtxt(coord_file, dtype=None, encoding='utf-8')
equatorialCoordinates = SkyCoord(coords[:,1].astype(str), coords[:,2].astype(str),
        frame='fk5', unit=(u.hourangle, u.deg))
indice = coords[:,0].astype(str)
equatorialCoordinates = np.array([equatorialCoordinates.ra.astype(float), equatorialCoordinates.dec.astype(float)]).T
axis1, axis2, angle = (float(args.beamshape[0])*scaling, float(args.beamshape[1])*scaling,
        180-float(args.beamshape[2]) if args.flip else float(args.beamshape[2]))
equatorialBoresight = SkyCoord(args.boresight[0], args.boresight[1],
        frame='fk5', unit=(u.hourangle, u.deg))
boresight = (equatorialBoresight.ra.deg , equatorialBoresight.dec.deg)
fileName = args.tiling_plot[0]

if args.inner is not None:
    inner = int(args.inner[0])
else:
    inner = 0

index = True

step = 1/10000000000.
wcs_properties = wcs.WCS(naxis=2)
wcs_properties.wcs.crpix = [0, 0]
wcs_properties.wcs.cdelt = [-step, step]
wcs_properties.wcs.crval = boresight
wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

center = boresight
resolution = step
inner_idx = []

fig = plt.figure(figsize=(1600./96, 1600./96), dpi=96)
axis = fig.add_subplot(111,aspect='equal', projection=wcs_properties)


scaled_pixel_coordinats = wcs_properties.wcs_world2pix(equatorialCoordinates, 0)
beam_coordinate = np.array(scaled_pixel_coordinats)
if inner > 0:
    index_sort = np.argsort(np.sum(np.square(beam_coordinate), axis=1))
    beam_coordinate = beam_coordinate.take(index_sort, axis=0)
    indice = indice.take(index_sort, axis=0)

for idx in range(len(beam_coordinate)):
    coord = beam_coordinate[idx]
    if index == True:
        num = indice[idx].split('cfbf')[-1]
        axis.text(coord[0], coord[1], int(num), size=6, ha='center', va='center')
    ellipse = Ellipse(xy=coord,
            width=2.*axis1/resolution,height=2.*axis2/resolution, angle=angle)
    ellipse.fill = False
    if inner > 0 and idx < inner:
        ellipse.fill = True
        ellipse.edgecolor = 'auto'
        inner_idx.append(int(num))
        if idx == 0:
            ellipse.facecolor = '#4169E1'
        else:
            ellipse.facecolor = '#0F52BA'
    axis.add_artist(ellipse)
margin = 1.1 * max(np.sqrt(np.sum(np.square(beam_coordinate), axis=1)))
axis.set_xlim(center[0]-margin, center[0]+margin)
axis.set_ylim(center[1]-margin, center[1]+margin)

if args.extra_source is not None:
    extra_coords = np.genfromtxt(args.extra_source[0], dtype=None)
    if len(extra_coords.shape) == 1:
        extra_coords = extra_coords.reshape(1, -1)
    extra_equatorial_coordinates = SkyCoord(extra_coords[:,1].astype(str),
            extra_coords[:,2].astype(str),
        frame='fk5', unit=(u.hourangle, u.deg))
    extra_equatorial_coordinates = np.array([extra_equatorial_coordinates.ra.astype(float),
            extra_equatorial_coordinates.dec.astype(float)]).T
    scaled_extra_pixel_coordinats = np.array(wcs_properties.wcs_world2pix(
                extra_equatorial_coordinates, 0))
    axis.scatter(scaled_extra_pixel_coordinats[:,0], scaled_extra_pixel_coordinats[:,1],s=40)

ra = axis.coords[0]
dec = axis.coords[1]
ra.set_ticklabel(size=20)
dec.set_ticklabel(size=20, rotation="vertical")
dec.set_ticks_position('l')
ra.set_ticks_position('b')
ra.set_axislabel("RA", size=20)
dec.set_axislabel("DEC", size=20)

plt.savefig(fileName, dpi=96)

if len(inner_idx) != 0:
    print(np.sort(inner_idx))
