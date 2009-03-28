#!/bin/sh
autoreconf --force --install -I config  -I m4
echo "Now run ./configure --with-boost=</path/to/boost> --prefix=$HOME && make install"
echo "Add --disable-shared to the configure line if building on Mac OS X"
