# Mosaic: Multibeamformed observation simulation and interferometry characterization

A software package consists of an interferometric pattern simulator and characterizer, an optimized tiling generator and a beamforming weights calculator. This document only describes the package from version 1.0.0. __Try Mosaic in your web browser [here](https://wchenastro.github.io/mosaic_web).__

## Dependent

For python 3.8.5+

- numpy
- scipy
- matplotlib
- astropy
- nvector
- katpoint

## Installation

Try Mosaic in your web browser [here](https://wchenastro.github.io/mosaic_web) without installing anything.

To use offline, install the package via pip:
```
pip3 install https://github.com/wchenastro/Mosaic/archive/refs/heads/master.zip
```

## Usage

There is a helper script `example/maketiling.py` to demonstrate the interface of the package

### Simulate the interferometric pattern and generate a fits file and a plot

```
python3 ./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 \
--datetime 2020.05.02 06:02:13.663903 --verbose --subarray 000, 001, 002:0.7, 003:0.5+0.1j \
--size 900 --resolution 40 --psf_plot psf.png --psf_fit psf.fits --weight
```

`--ants`: the file containing the antenna specification.

`--freq`: the frequency at which the interferometric pattern is simulated (Hz).

`--source`: the equatorial coordinates of the source in `hh:mm:ss.s dd:mm:ss.s` format.

`--datetime`: the date and time of the observation in UTC and in `yyyy.mm.dd hh:mm:ss.s` format.

 `--subarray`: a list of index for selection of antennas in the file specified by `--ants`. Optional scale or complex weight can be attached after each antenna index separated by a colon.

`--resolution`: the scale of one single pixel in the pattern in arc seconds, default is None which means it is determined by the code.

`--size`: the total number of pixels in the simulation, default is 400 which corresponds to a pattern of 20x20 in dimension.

`--field`: the side length of the field of view in sexagesimal unit, such as `60s` or `20m` or `1d`, default is None. `--resolution`, `--size` and `--field` can not be used all together, it can be combination of any two.

`--psf_plot`: filename of the plot of the pattern, the file format can be anything that matplotlib supports, such as "jpeg, pdf".

`--psf_fits`: filename of the fits file of the pattern

`--weight`: a switch for individual weight for each antenna, the weight values in `--subarray` will not be effective without this argument.

`--verbose`: print logs containing the input parameter and result, the input parameter listed in the log should reproduce the same result.

Example output:
![psf](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/psf.png)

### Generate a tiling in specified overlap ratio and overlay some point sources on top of it

```
python3 ./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 \
--datetime 2020.05.02 06:02:13.663903 --beamnum 400 --verbose --overlap 0.7 \
--subarray 000, 017, 036, 038, 041, 043, 044 --tiling_method variable_size \
--tiling_shape circle --tiling_plot tiling.png --overlay_source overlay_sources
```

`--beamnum`: the requesting beam number in the tiling, the actual number in the generated tiling is less than or equal to this number. The default is 400.

`--tiling_method`: the method to use for generating the tiling, possible values are

- "`variable_size`": given an overlap ratio, the code decide the size of the tiling
- "`variable_ovelap`": given a size of the tiling, the code decide the overlap between the beams.

`--overlap`: The beams in the tiling overlap with each other in their power levels equal to this ratio, only available in the "`variable_size`" method. The default is 0.5.

`--tiling_shape`: the shape of the tiling boundaries, possible values are: "circle", "hexagon", "ellipse", "polygon", "annulus". The "`variable_size`" method only supports the first two shapes.

`--tiling_plot`: the filename for the plot of the tiling, the file format can be anything that matplotlib supports, such as "jpeg, pdf".

`--overlay_source`: the file containing the point sources to overlay, one per line,  in `identification RA DEC` format. for example: "C 00:23:50.3546 -72:04:31.5048"

Example output:
![circular_tiling](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/circular_tiling.png)

### Generate an elliptical shape tiling,  let the code decide a suitable overlap and output the coordinates

```
python3 ./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 \
--datetime 2020.05.02 06:02:13.663903 --beamnum 400 --verbose --subarray 000, 001, 002, 003 \
--tiling_method variable_overlap --tiling_shape ellipse --tiling_parameter 0.07 0.05 45 \
--tiling_plot tiling.png --tiling_coordinate coordinate.csv
```

`--tiling_coordinate`: the filename for the equatorial coordinates in degrees.

`--tiling_parameter_coordinate_type`: the coordinate type of the parameter,  default is image coordinate.

`--tiling_parameter`: the parameter of the tiling, for example:

- "`--tiling_shape circle --tiling_parameter 0.05`": a circular shape tiling with radius of 0.05 degree
-  "`--tiling_shape hexagon --tiling_parameter 0.07 45`": a hexagonal shape tiling with its circumradius and orientation in degrees
-  "`--tiling_shape ellipse --tiling_parameter 0.07 0.05 45`": an elliptical shape tiling with its two semi-axis and orientation in degrees
-  "`--tiling_shape polygon --tiling_parameter 6.1522476, -72.0506681, 5.9448280, -72.0557907, 5.8695621, -72.0879815, 6.0670744, -72.1139826`": a polygonal shape tiling with its vertices in "RA1, DEC1, RA2, DEC2, RA3, DEC3" format.

Example output:
![elliptical_tiling](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/elliptical_tiling.png)


### Generate a polygon shape tiling using a boundary region file and generate a region file for all the beams

```
python3 ./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 \
--datetime 2020.05.02 06:02:13.663903 --beamnum 400 --verbose --overlap 0.7 \
--subarray 000, 001, 002, 003 --tiling_method variable_overlap --tiling_shape polygon \
--tiling_parameter_file polygon.reg --tiling_region tiling.reg
```

`--tiling_parameter_file`: the filename of the polygon boundary region file from ds9

`--tiling_region`: the filename for the region file of the generated tiling which can be imported into ds9

Example:

|                                       Create a region in ds9                                        |                                        Create a tiling within the region                                        |
| :-------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------: |
| ![ds9 region](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/ds9_region.png) | ![ds9 region with tiling](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/ds9_region_tiling.png) |

### Generate an annulus shape tiling with different shape of boundaries.

```
python3 ./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 --datetime 2020.05.02 06:02:13.663903 --beamnum 400 --verbose --subarray 000, 001, 002, 003 --tiling_method variable_overlap --tiling_shape annulus --tiling_parameter polygon 9.00,-72.5,8.5,-71.2,3,-71.5,2,-73:ellipse 0.4 0.6 100 --tiling_plot tiling.png --tiling_coordinate coordinate.csv
```

`--tiling_parameter polygon 9.00,-72.5,8.5,-71.2,3,-71.5,2,-73:ellipse 0.4 0.6 100`: an annulus shape with a polygon as the outer boundary and an ellipse as the inner boundary. The parameters of outer and inner boundaries are separated with an "`:`", each set of parameters starts with the name of the shape. Currently, only boundaries of polygon and ellipse shape are supported. If the vectors of the polygon are provided by a file, then it can be specified as `--tiling_parameter polygon:ellipse 0.4 0.6 100 --tiling_parameter_file polygon.reg`

Example:
![annulus_tiling](https://gist.githubusercontent.com/wchenastro/eb0159359511808ff7d0363db9b32d8b/raw/d57ba74209627de65307f99809d1fa92e34b8c79/annulus_tiling.png)
## License

[MIT](https://github.com/wchenastro/Mosaic/blob/master/LICENSE)

