# Program SPHERES

SPHERES is a program package written in C++ for optimising and analysing cluster structures.

## Compilation

The program was developed using the clang++-6.0 compiler on Ubuntu 18.04 LTS.

### Prerequisits

The program uses the machine learning library [dlib](https://github.com/davisking/dlib) to carry out optimisations. Version 19.2 should be installed on your system.

### Installation

Open the included Makefile and set the `MATLIB` variable to the dlib-19.2 installation directory.

```
make
```

Your are done!




