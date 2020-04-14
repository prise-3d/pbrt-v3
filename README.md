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

Install procedure with libtorch:
```
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.4.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.4.0+cpu.zip

mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=../libtorch ..
make -j
```

### Updates

Possibility to generate few images in one run with different random seed. 

```sh
cd build
./pbrt --images 10 --samples 100 --folder output ../path/to/file.pbrt
```

**Parameters:**
- *--samples*: number of samples per pixel per image
- *--images*: number of images with `n` samples expected as output
- *--folder*: specific output folder (by default `temp`)
- *--startindex*: specific image index to start with (if necessary to relaunch pbrt) and to not erase previously generated image
- *--independent*: boolean which specifies if each generated image are independent or image convergence is used (samples merge at each new image generated)
- *--digit*: number of digits used when saving image index (if 5, index in filename is equal to `00001`)
- *--normals*: boolean which specifies if normals map is saved after generating first image
- *--zbuffer*: boolean which specifies if Z-Buffer is saved after generating first image
- *--model_path*: pytorch model to use (need `.pt` file), follow this [tutorial](https://pytorch.org/tutorials/advanced/cpp_export.html) to extract your model.

Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-{imageIndex}.{ext}`

Normals Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-normals.rawls`

Z-Buffer Images are saved into `{folder}/{image_name}` folder with the following name convention `{image_name}-S{samples}-zbuffer.rawls`

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