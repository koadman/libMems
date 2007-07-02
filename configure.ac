dnl Process this file with autoconf to produce a configure script.
AC_INIT(libMems/Match.h)

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

GENERIC_LIBRARY_NAME=libMems

#release versioning
GENERIC_MAJOR_VERSION=1
GENERIC_MINOR_VERSION=5
GENERIC_MICRO_VERSION=0

#API version (often = GENERIC_MAJOR_VERSION.GENERIC_MINOR_VERSION)
GENERIC_API_VERSION=1.5
AC_SUBST(GENERIC_API_VERSION)

#shared library versioning
GENERIC_LIBRARY_VERSION=1:0:0
#                       | | |
#                +------+ | +---+
#                |        |     |
#             current:revision:age
#                |        |     |
#                |        |     +- increment if interfaces have been added
#                |        |        set to zero if interfaces have been removed
#                                  or changed
#                |        +- increment if source code has changed
#                |           set to zero if current is incremented
#                +- increment if interfaces have been added, removed or changed



dnl --------------------------------
dnl Package name and version number
dnl --------------------------------

AC_SUBST(GENERIC_LIBRARY_VERSION)

PACKAGE=$GENERIC_LIBRARY_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)

GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_VERSION)

VERSION=$GENERIC_VERSION

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define)



AC_PREFIX_DEFAULT(/usr/local)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S 
AM_PROG_LIBTOOL
AC_SYS_LARGEFILE

dnl Checks for header files.
AC_HEADER_STDC

dnl Check what compiler we're using
AM_CONDITIONAL(ICC, test x$CXX = xicc )

dnl Allow debugging compilation
AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

dnl Get location of Boost
AX_PATH_BOOST(1.33)
AC_CXX_NAMESPACES
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS

dnl Get location of libGenome/libClustalW Headers
PKG_CHECK_MODULES(DEPS, libGenome-1.3 >= 1.0.0  libClustalW-1.0 >= 1.8.4  libMUSCLE-3.6 >= 1.0.0)
AC_SUBST(DEPS_CFLAGS)

dnl ensure portability for OS X with these checks
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(sys/aio.h)
AC_CHECK_HEADERS(aio.h)
AC_CHECK_HEADERS(features.h)

dnl certain parts of the library need async io and threads
AC_CHECK_LIB(pthread, pthread_mutex_unlock)
AC_CHECK_LIB(rt, aio_write)
AC_SUBST(DEPS_LIBS)


AM_CONFIG_HEADER(config.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
dnl AC_C_BIGENDIAN
AC_HEADER_TIME

dnl Checks for library functions.
AC_PROG_GCC_TRADITIONAL

dnl SAVE_LIBRARY_VERSION
AC_SUBST(LIBTOOL_VERSION_INFO)

AC_OUTPUT(Makefile libMems/Makefile libMems/dmSML/Makefile doc/Makefile libMems-1.5.pc )
#doc/html/Makefile  doc/man/Makefile  doc/man/man3/Makefile)
