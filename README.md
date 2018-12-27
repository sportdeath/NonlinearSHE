# NonlinearSHE

An implementation of the algorithms described in the included [paper](https://github.com/sportdeath/NonlinearSHE/raw/master/paper/paper.pdf).
The algorithms allow for efficient computation of nonlinear functions, like division, boolean comparisons, logarithms and exponentials under homomorphic encryption.
These computations are built solely from these multiplication and addition operations exposed by homomorphic encryption.
The algorithms described and implemented are provably optimal as they meet an asymptotic lower bound on the number of multiplications performed, which is the bottleneck in most if not all homomorphic encryption schemes.
The encryption scheme used is the somewhat homomorphic encryption scheme [YASHE](https://eprint.iacr.org/2013/075.pdf).

## Dependencies
This library relies on the number theoretic library [NTL](http://www.shoup.net/ntl/) for large integer and polynomial math as well as [Boost](http://www.boost.org/) for serialization (storing and retrieving information from disk). Both are available through most pacakge managers (apt-get, pacman, brew). The library is built with [CMake](https://cmake.org/). The unit tests are evaluated with [Google Test](https://github.com/google/googletest), but the framework is downloaded on the fly through CMake and not required as a dependency.

Several examples can be built optionally that depend on [CImg](http://www.cimg.eu), a basic C++ image processing library.

## Building

Install the dependencies, then clone the library with git submodules:
    
    git clone --recurse-submodules https://github.com/sportdeath/NonlinearSHE

The library can be built with CMake as follows:

    cd NonlinearSHE
    mkdir build
    cd build
    cmake ..
    make

To build the library with example programs you can instead add the flag

    cmake -DBUILD_EXAMPLES=ON ..
    make

## Running Tests
The unit tests can be run simply by running the following command in the build directory:

    ctest

For more indepth output, the tests can be run individually

    ./numberTheory_test
    ./functions_test
    ./yashe_test

## Running Examples

A set of examples involving homomorphic image processing are located in the examples folder.
If the library has been made with `-DBUILD_EXAMPLES=ON`, they will be built to `build/examples/FHEImageProcessing`.
Some appropriately sized pictures can be found in the `resources` subfolder.

### Transforming an RGB image into a YCbCr image

The following example encrypts an image as RGB channels and then homomorphically transforms
it into YCbCr channels and displays the result.

    cd build/examples/FHEImageProcessing
    ./RGBtoYCbCr resources/mona.png

For an image with no more than 5376 pixels represented in 24-bit color and a security parameter of 128, this process completes in under 4 minutes:

* 30 seconds to generate the YASHE parameters (requires polynomial factorization)
* 15 seconds to encrypt the image
* 155 seconds to compute the homomorphic transformation
* 30 seconds to decrypt the image

The input file displayed next to the individual channels Y, Cb, and Cr of the output and their combination

![RGBtoYCbCr](https://github.com/sportdeath/NonlinearSHE/raw/master/images/RGBtoYCbCr.png)

### Applying a HSV filter to an RGB image

The following example encrypts an image as 8 bit RGB and then homomorphically transforms it
into HSV (hue, saturation, value). Then it rotates the hue and transforms the image back to RGB and displays the result.

    cd build/examples/FHEImageProcessing
    ./imageTransform resources/marilyn8Bit.png

Since this example only uses 8 bit color, it is significantly faster than the other examples and completes in under a minute:

* 30 seconds to encrypt the image
* 17 seconds to compute the homomorphic transformation
* 9 seconds to decrypt the image

The input image displayed next to the output image:

![ColorTransform](https://github.com/sportdeath/NonlinearSHE/raw/master/images/ColorTransform.png)

### Taking the mean of two images
The following example encrypts two images as 24-bit RGB images, then takes the mean of their pixel values and displays the result.

    cd build/examples/FHEImageProcessing
    ./mean resources/mona.png resources/lena.jpg 

The result is computed in about 3 minutes. The input images displayed next to the output image:

![Mean](https://github.com/sportdeath/NonlinearSHE/raw/master/images/Mean.png)
