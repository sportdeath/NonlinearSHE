# YASHE
An implementation of the somewhat homomorphic encryption scheme [YASHE](https://eprint.iacr.org/2013/075.pdf).

## Dependencies
This library relies on the number theoretic library [NTL](http://www.shoup.net/ntl/) for large integer and polynomial math as well as [Boost](http://www.boost.org/) for serialization (storing and retrieving information from disk). Both are available through most pacakge managers (apt-get, pacman, brew). The library is built with [CMake](https://cmake.org/). The unit tests are evaluated with [Google Test](https://github.com/google/googletest), but the framework is downloaded on the fly through CMake and not required as a dependency.

Several examples can be built optionally that depend on [CImg](http://www.cimg.eu), a basic C++ image processing library.

## Installation
The library can be built with CMake as follows:
```bash
cd (directory of YASHE)
mkdir build
cd build
cmake ..
make
```
To build the library with example programs you can run
```bash
cmake -DBUILD_EXAMPLES=ON ..
```
## Running Tests
The test cases can be run simply by running the following command in the build directory:
```bash
ctest
```
For more indepth output, the tests can be run individually
```bash
./yashe_test
```

## Running Examples
A set of examples involving homomorphic image processing are located in the examples folder.
If the library has been made with `-DBUILD_EXAMPLES=ON`, they will be built to `build/examples/FHEImageProcessing`.
Some appropriately sized pictures can be found in the `resources` subfolder.

The following example encrypts an image as RGB channels and then homomorphically transforms
it into YCbCr channels and displays the result.
```bash
./RGBtoYCbCr resources/mona.png
```
For an image with no more than 5376 pixels represented in 24-bit color and a security parameter of 128, this process completes in about 5 minutes:
* 100 seconds to encrypt the image
* 170 seconds to compute the homomorphic transformation
* 30 seconds to decrypt the image
