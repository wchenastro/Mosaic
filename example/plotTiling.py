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

index = True

step = 1/10000000000.
wcs_properties = wcs.WCS(naxis=2)
wcs_properties.wcs.crpix = [0, 0]
wcs_properties.wcs.cdelt = [-step, step]
wcs_properties.wcs.crval = boresight
wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

center = boresight
resolution = step

fig = plt.figure(figsize=(1600./96, 1600./96), dpi=96)
axis = fig.add_subplot(111,aspect='equal', projection=wcs_properties)


scaled_pixel_coordinats = wcs_properties.wcs_world2pix(equatorialCoordinates, 0)
beam_coordinate = np.array(scaled_pixel_coordinats)

for idx in range(len(beam_coordinate)):
    coord = beam_coordinate[idx]
    ellipse = Ellipse(xy=coord,
            width=2.*axis1/resolution,height=2.*axis2/resolution, angle=angle)
    ellipse.fill = False
    axis.add_artist(ellipse)
    if index == True:
        num = indice[idx].split('cfbf')[-1]
        axis.text(coord[0], coord[1], int(num), size=6, ha='center', va='center')
margin = 1.1 * max(np.sqrt(np.sum(np.square(beam_coordinate), axis=1)))
axis.set_xlim(center[0]-margin, center[0]+margin)
axis.set_ylim(center[1]-margin, center[1]+margin)

ra = axis.coords[0]
dec = axis.coords[1]
ra.set_ticklabel(size=20)
dec.set_ticklabel(size=20, rotation="vertical")
dec.set_ticks_position('l')
ra.set_ticks_position('b')
ra.set_axislabel("RA", size=20)
dec.set_axislabel("DEC", size=20)

plt.savefig(fileName, dpi=96)
