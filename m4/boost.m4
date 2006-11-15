dnl
dnl A macro to detect the boost install path and version information
dnl Taken from the boost mailing list by Aaron Darling
dnl

AC_DEFUN([AX_PATH_BOOST],
[
	AC_ARG_WITH([boost],
			AC_HELP_STRING([--with-boost=DIR],
						   [use boost (default is NO) , specify the root directory for boost library (optional)]), [
		if test "$withval" = "no"; then
			want_boost="no"
		elif test "$withval" = "yes"; then
			want_boost="yes"
			ac_boost_path=""
		else
			want_boost="yes"
			ac_boost_path="$withval"
		fi
	],[want_boost="no"])

	if test "x$want_boost" = "xyes"; then
		boost_lib_version_req=ifelse([$1], ,1.20.0,$1)
		boost_lib_version_req_shorten=`expr $boost_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
		boost_lib_version_req_major=`expr $boost_lib_version_req : '\([[0-9]]*\)'`
		boost_lib_version_req_minor=`expr $boost_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
		boost_lib_version_req_sub_minor=`expr $boost_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
		if test "x$boost_lib_version_req_sub_minor" = "x" ; then
			boost_lib_version_req_sub_minor="0"
	    fi
		WANT_BOOST_VERSION=`expr $boost_lib_version_req_major \* 100000 \+ $boost_lib_version_req_minor \* 100 \+ $boost_lib_version_req_sub_minor`
		AC_MSG_CHECKING(for boostlib >= $boost_lib_version_req)
		succeeded=no

		dnl first we check the system location for boost libraries
		dnl this location ist chosen if boost libraries are installed with the --layout=system option
		dnl or if you install boost with RPM
		if test "$ac_boost_path" != ""; then
			BOOST_LDFLAGS="-L$ac_boost_path/lib"
			BOOST_CPPFLAGS="-I$ac_boost_path/include"
		else
			for ac_boost_path_tmp in /usr /usr/local /opt ; do
				if test -d "$ac_boost_path_tmp/include/boost" && test -r "$ac_boost_path_tmp/include/boost"; then
					BOOST_LDFLAGS="-L$ac_boost_path_tmp/lib"
					BOOST_CPPFLAGS="-I$ac_boost_path_tmp/include"
					break;
				fi
			done
	
		fi
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS
		
		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

	     AC_TRY_COMPILE(
	       [
#include <boost/version.hpp>
],
       [
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
# error Boost version is too old
#endif

		],
	    [
	 AC_MSG_RESULT(yes ($_version))
		 succeeded=yes
		 found_system=yes
	 ifelse([$2], , :, [$2])
       ],
       [
       ])

		dnl if we found no boost with system layout we search for boost libraries
		dnl built and installed without the --layout=system option or for a staged(not installed) version
		
		if test "x$succeeded" != "xyes"; then
			_version=0
			if test "$ac_boost_path" != ""; then
				BOOST_LDFLAGS="-L$ac_boost_path/lib"
				if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
					for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
						_version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
						V_CHECK=`expr $_version_tmp \> $_version`
						if test "$V_CHECK" = "1" ; then
							_version=$_version_tmp
						fi
						VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
						BOOST_CPPFLAGS="-I$ac_boost_path/include/boost-$VERSION_UNDERSCORE"
					done
				fi
			else
				for ac_boost_path in /usr /usr/local /opt ; do
					if test -d "$ac_boost_path" && test -r "$ac_boost_path"; then
						for i in `ls -d $ac_boost_path/include/boost-* 2>/dev/null`; do
							_version_tmp=`echo $i | sed "s#$ac_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
							V_CHECK=`expr $_version_tmp \>$_version`
							if test "$V_CHECK" = "1" ; then
								_version=$_version_tmp
								best_path=$ac_boost_path
							fi
						done
					fi
				done

				VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
				BOOST_CPPFLAGS="-I$best_path/include/boost-$VERSION_UNDERSCORE"
				BOOST_LDFLAGS="-L$best_path/lib"

			    if test "x$BOOST_ROOT" != "x"; then
					if test -d "$BOOST_ROOT" && test -r "$BOOST_ROOT"; then
						version_dir=`expr //$BOOST_ROOT : '.*/\(.*\)'`
						stage_version=`echo $version_dir | sed 's/boost_//' | sed 's/_/./g'`
						stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
						V_CHECK=`expr $stage_version_shorten \>\=$_version`
						if test "$V_CHECK" = "1" ; then
							AC_MSG_NOTICE(We will use a staged boost library from $BOOST_ROOT)
							BOOST_CPPFLAGS="-I$BOOST_ROOT"
							BOOST_LDFLAGS="-L$BOOST_ROOT/stage/lib"
						fi
					fi
			    fi
			fi

			CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
			export CPPFLAGS
			LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
			export LDFLAGS

		     AC_TRY_COMPILE(
		       [
#include <boost/version.hpp>
],
       [
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
# error Boost version is too old
#endif

		],
	    [
	 AC_MSG_RESULT(yes ($_version))
		 succeeded=yes
	 ifelse([$2], , :, [$2])
       ],
       [
	 AC_MSG_RESULT(no ($_version))
	 ifelse([$3], , :, [$3])
       ])
		fi

		if test "$succeeded" != "yes" ; then
			if test "$_version" = "0" ; then
				AC_MSG_ERROR('We could not detect the boost libraries.
			   If you have a staged boost library (still not installed)
			   please specify \$BOOST_ROOT in your environment and do not
			   give a PATH to --with-boost option')
			else
				AC_MSG_ERROR('Your boost libraries seems to old (version $_version).
			   We need at least $boost_lib_version_shorten')
			fi
		else
			AC_MSG_NOTICE(BOOST_CPPFLAGS: $BOOST_CPPFLAGS)
			AC_MSG_NOTICE(BOOST_LDFLAGS: $BOOST_LDFLAGS)
			AC_SUBST(BOOST_CPPFLAGS)
			AC_SUBST(BOOST_LDFLAGS)
			AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])

			AC_CACHE_CHECK(whether the Boost::Filesystem library is available,
						   ax_cv_boost_filesystem,
						[AC_LANG_SAVE
			 AC_LANG_CPLUSPLUS
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include <boost/filesystem/path.hpp>]],
				   [[using namespace boost::filesystem;
				   path my_path( "foo/bar/data.txt" );
				   return 0;]]),
						   ax_cv_boost_filesystem=yes,ax_cv_boost_filesystem=no)
								 AC_LANG_RESTORE
			])
			if test "$ax_cv_boost_filesystem" = "yes"; then
				AC_DEFINE(HAVE_BOOST_FILE,,[define if the Boost::FILESYSTEM library is available])
				if test "x$found_system" = "xyes"; then
					ax_lib="boost_filesystem"
				else
					ax_lib="boost_filesystem-$CC"
				fi
				AC_MSG_NOTICE(ax_lib=$ax_lib)

			    AC_CHECK_LIB($ax_lib, main, [BOOST_FILESYSTEM_LIB=$ax_lib AC_SUBST(BOOST_FILESYSTEM_LIB) break],[link_filesystem="no"])
				if test "x$link_filesystem" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi

			AC_CACHE_CHECK(whether the Boost::Program_Options library is available,
						   ax_cv_boost_program_options,
						[AC_LANG_SAVE
			 AC_LANG_CPLUSPLUS
			 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include <boost/program_options.hpp>]],
				   [[boost::program_options::options_description generic("Generic options");
				   return 0;]]),
		   ax_cv_boost_program_options=yes,ax_cv_boost_program_options=no)
			 AC_LANG_RESTORE
			])
			if test "$ax_cv_boost_program_options" = yes; then
				AC_DEFINE(HAVE_BOOST_PROGRAM_OPTIONS,,[define if the Boost::PROGRAM_OPTIONS library is available])
				if test "x$found_system" = "xyes"; then
					ax_lib="boost_program_options"
				else
					ax_lib="boost_program_options-$CC"
				fi
				AC_MSG_NOTICE(ax_lib=$ax_lib)

			    AC_CHECK_LIB($ax_lib, main, [BOOST_PROGRAM_OPTIONS_LIB=$ax_lib AC_SUBST(BOOST_PROGRAM_OPTIONS_LIB) break], [link_program_options="no"])
				if test "x$link_program_options" = "xno"; then
					AC_MSG_NOTICE(Could not link against $ax_lib !)
				fi
			fi
		fi
	fi
])
 
