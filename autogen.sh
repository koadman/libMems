#!/bin/sh
(libtoolize || glibtoolize) && \
aclocal -I m4 && \
autoheader && \
automake -a && \
autoconf
echo "Now run ./configure --with-boost=</path/to/boost> --prefix=$HOME && make install"
echo "Add --disable-shared to the configure line if building on Mac OS X"
echo "If boost doesn't create default library links, you may need to add"
echo "--with-boost-filesystem=gcc --with-boost-program-options=gcc"
