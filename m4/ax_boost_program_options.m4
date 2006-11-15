dnl @synopsis AX_BOOST_PROGRAM_OPTIONS
dnl
dnl This macro checks to see if the Boost.ProgramOptions library is
dnl installed. It also attempts to guess the currect library name using
dnl several attempts. It tries to build the library name using a user
dnl supplied name or suffix and then just the raw library.
dnl
dnl If the library is found, HAVE_BOOST_PROGRAM_OPTIONS is defined and
dnl BOOST_PROGRAM_OPTIONS_LIB is set to the name of the library.
dnl
dnl This macro calls AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB).
dnl
dnl @category InstalledPackages
dnl @author Thomas Porschberg <thomas@randspringer.de>
dnl @version 2005-12-28
dnl @license GPLWithACException

AC_DEFUN([AX_BOOST_PROGRAM_OPTIONS],
[AC_REQUIRE([AC_CXX_NAMESPACES])dnl
AC_CACHE_CHECK(whether the Boost::Program_Options library is available,
ax_cv_boost_program_options,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include <boost/program_options.hpp>]],
                                   [[boost::program_options::options_description generic("Generic options");
                                   return 0;]]),
                   ax_cv_boost_program_options=yes, ax_cv_boost_program_options=no)
 AC_LANG_RESTORE
])
if test "$ax_cv_boost_program_options" = yes; then
  AC_DEFINE(HAVE_BOOST_PROGRAM_OPTIONS,,[define if the Boost::PROGRAM_OPTIONS library is available])
  dnl Now determine the appropriate file names
  AC_ARG_WITH([boost-program-options],AS_HELP_STRING([--with-boost-program-options],
  [specify the boost program-options library or suffix to use]),
  [if test "x$with_boost_program_options" != "xno"; then
    ax_program_options_lib=$with_boost_program_options
    ax_boost_program_options_lib=boost_program_options-$with_boost_program_options
    axpols=$ax_program_options_lib-s
    axbpols=$ax_boost_program_options_lib-s
  fi])
  for ax_lib in $ax_program_options_lib $ax_boost_program_options_lib boost_program_options $axpols $axbpols boost_program_options-s; do
    AC_CHECK_LIB($ax_lib, main, [BOOST_PROGRAM_OPTIONS_LIB=$ax_lib break])
  done
  AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB)
fi
])dnl
