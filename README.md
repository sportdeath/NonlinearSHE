# YASHE
An implementation of the somewhat homomorphic encryption scheme [YASHE](https://eprint.iacr.org/2013/075.pdf).

# Dependencies
This library relies on the number theoretic library [NTL](http://www.shoup.net/ntl/) for large integer and polynomial math as well as [Boost](http://www.boost.org/) for serialization (storing and retrieving information from disk). Both are available through most pacakge managers (apt-get, pacman, brew). The library is built with [CMake](https://cmake.org/). The unit tests are evaluated with [Google Test](https://github.com/google/googletest), but the framework is downloaded on the fly through CMake and not required as a dependency.

Several examples can be built optionally that depend on [CImg](http://www.cimg.eu), a basic C++ image processing library.

# Installation
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
cmake -DBUILD_EXAMPLES ..
```
