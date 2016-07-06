# YASHE
An implementation of the somewhat homomorphic encryption scheme [YASHE](https://eprint.iacr.org/2013/075.pdf).

# Dependencies
This library relies on the number theoretic library [NTL](http://www.shoup.net/ntl/) for large integer and polynomial math as well as [Boost](http://www.boost.org/) for serialization. Both are available through most pacakge managers (apt-get, pacman, brew). The library is built with CMake.

# Installation
Most common:

cd <directory of YASHE>
mkdir build
cd build
cmake ../
make
