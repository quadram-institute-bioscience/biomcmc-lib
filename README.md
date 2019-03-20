# Biomcmc-lib 

This library has low level functions used by other C programs that work with phylogenetic analyses. 
It borrows many functions from the [guenomu software](https://bitbucket.org/leomrtns/guenomu/) for phylogenomic species tree inference and 
also some ideas from the [genefam-dist library](https://github.com/leomrtns/genefam-dist).

It includes a version of [argtable3](https://www.argtable.org/) downloaded from
[github](https://github.com/argtable/argtable3) on 2019.03.12. Argtable3 is released under a BSD license. 

## Installation
This library usually is not installed directly, but as a submodule of another project. It includes the 'makefile.am' for
autotools.

## License 
Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

biomcmc-lib is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

The files `argtable3.c` and `argtable3.h` are part of the [argtable3](https://www.argtable.org/) library, maintained by
Tom G. Huang at the time of this writting  and are distributed under a BSD license. Please refer to 
https://github.com/argtable/argtable3 for the list of authors.
