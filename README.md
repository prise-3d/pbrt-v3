PrISE-3D - pbrt, Version 3
==========================

This repository holds the source code to the version of pbrt for the PrISE-3D ANR project. 
PBRT is described in the third edition of *Physically Based Rendering: From
Theory to Implementation*, by [Matt Pharr](http://pharr.org/matt), [Wenzel
Jakob](http://www.mitsuba-renderer.org/~wenzel/), and Greg Humphreys.  As
before, the code is available under the BSD license.

The [pbrt website](http://pbrt.org) has general information about both the
*Physically Based Rendering* book as well as many other resources for pbrt.
As of October 2018, the full [text of the book](http://www.pbr-book.org) is
now available online, for free.

Please refer to the pbrt-v3 github repository for further information about pbrt.

Clone the project
-----------------

To check out pbrt together with all dependencies, be sure to use the
`--recursive` flag when cloning the repository, i.e.
```bash
$ git clone --recursive https://github.com/mmp/pbrt-v3/
```
If you accidentally already cloned pbrt without this flag (or to update an
pbrt source tree after a new submodule has been added, run the following
command to also fetch the dependencies:
```bash
$ git submodule update --init --recursive
```

### Makefile builds (Linux, other Unixes, and Mac) ###

```
mkdir build
cd build
cmake ..
make -j
```

### Updates

#### Running `pbrt` 

Possibility to generate few images in one run with different random seed. Hence, scene is loading only once time.

```sh
cd build
./pbrt --images 10 --samples 100 --folder output ../path/to/file.pbrt
```

**Parameters:**
- `--samples`: number of samples per pixel per image
- `--images`: number of images with `n` samples expected as output
- `--folder`: specific output folder (by default `temp`)
- `--startindex`: specific image index to start with (if necessary to relaunch pbrt) and to not erase previously generated image
- `--independent`: boolean which specifies if each generated image are independent or image convergence is used (samples merge at each new image generated)
- `--digit`: number of digits used when saving image index (if 5, index in filename is equal to `00001`)
- `monk`: number `k` clusters to use for Median of meanNs

**Output information:**

- Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-{imageIndex}.{ext}`
- Normals Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-normals.rawls`
- Z-Buffer Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-zbuffer.rawls`

#### Information

All updates made into this version are encapsulated by these comments:
```cpp
//////////////////////
// PrISE-3D Updates //
//////////////////////

...
...

//////////////////////////
// END PrISE-3D Updates //
//////////////////////////
```