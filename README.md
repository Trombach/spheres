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

## Run

```
optimize [filename]
```
`[filename]` must contain a list of 3D coordinates. Individual structures must be separated by a single blank line. The file must end in a single blank line.

```
analyze
```
Automatically opens the file `coords` produced by `optimize [filename]`.

```
match [filename1] [filename2]
```
`[filename1]` and `[filename2]` must be of the same format as `[filename]`.

```
global [number-of-atoms] [number-of-steps]
```
Starts a basinhopping run with the given parameters.





