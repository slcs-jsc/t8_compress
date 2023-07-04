# t8_compress

## Features

- Software package to learn how to apply the adaptive mesh library `t8code` for compression of 2-D or 3-D meteorological data fields.

## Installation

- Requires an installation of `t8code` in a separate directory, see https://github.com/DLR-AMR/t8code for details.

## Compilation

- Edit the `Makefile` according to your needs. In particular, set `INCDIR` and `LIBDIR` to point to the `t8code` library.

- Use `make` to compile the code.

## Running a test case

- Use `make check` to run a test case.

- Results can be inspected with ParaView, e.g. `paraview forest_with_data.pvtu`.
