# ising

Ferromagnetic and antiferromagnetic Ising model simulator on a hypercubic lattice. 

This repository was forked from the Ising model simulator [here](https://github.com/zeehio/ising). 
The primary differences from the parent repository are the following:
- This code supports the antiferromagnetic case.
- This code **only** supports the Wolff MCMC algorithm.


## Building the source code

Note: This will need autotools (`autoconf`, `automake`, and `libtool`).

Run the following:

    ./autogen.sh
    ./configure
    make
    make install

The **ising** program will then be found in the `install/bin` directory.


# Usage

Move to the default install directory:

    cd install/bin

Run the following to get a list of command-line arguments and to see an example:

    ./ising --help
