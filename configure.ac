# Process this file with reconfigure.sh to update configure script.
# See also config/config.*.


# General Configurations

AC_INIT([lis],[2.1.6],[devel@ssisc.org])
AC_CONFIG_SRCDIR([src/fortran/lisf_matrix.c])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([include/lis_config.h])

AM_MAINTAINER_MODE

AC_ARG_VAR(MPICC,[MPI C compiler])
AC_ARG_VAR(MPIF77,[MPI Fortran compiler])
AC_ARG_VAR(MPIRUN,[command to execute MPI programs])
AC_ARG_VAR(OMPFLAG,[OpenMP flags])
AC_ARG_VAR(MPFLAG,[floating point precision flags])
AC_ARG_VAR(TARGET,[system architecture])


# Build Options

AC_ARG_ENABLE(fortran, AS_HELP_STRING([--enable-fortran],[enable FORTRAN 77 compatible interfaces [[default=no]]]))
AC_ARG_ENABLE(f90, AS_HELP_STRING([--enable-f90],[enable Fortran 90 compatible interface [[default=no]]]))
AC_ARG_ENABLE(saamg, AS_HELP_STRING([--enable-saamg],[enable SA-AMG preconditioner [[default=no]]]))
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi],[build with MPI library [[default=no]]]))
AC_ARG_ENABLE(omp, AS_HELP_STRING([--enable-omp],[build with OpenMP library [[default=no]]]))
AC_ARG_ENABLE(quad, AS_HELP_STRING([--enable-quad],[enable quadruple precision operations [[default=no]]]))
AC_ARG_ENABLE(sse2, AS_HELP_STRING([--enable-sse2],[use Intel Streaming SIMD Extensions [[default=yes]]]),[],enable_sse2=yes)
AC_ARG_ENABLE(fma, AS_HELP_STRING([--enable-fma],[use fused multiply add [[default=no]]]))
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[enable debugging [[default=no]]]))
AC_ARG_ENABLE(gprof, AS_HELP_STRING([--enable-gprof],[enable profiling [[default=no]]]))
AC_ARG_ENABLE(64bit, AS_HELP_STRING([--enable-64bit],[build 64bit binary on AIX and Solaris x86 [[default=no]]]))
AC_ARG_ENABLE(32bit, AS_HELP_STRING([--enable-32bit],[build 32bit binary on Mac OS X [[default=no]]]))
AC_ARG_ENABLE(longlong, AS_HELP_STRING([--enable-longlong],[enable long long int support [[default=no]]]))
AC_ARG_ENABLE(longdouble, AS_HELP_STRING([--enable-longdouble],[enable long double support [[default=no]]]))
AC_ARG_ENABLE(complex, AS_HELP_STRING([--enable-complex],[enable complex scalar support [[default=no]]]))
AC_ARG_ENABLE(test, AS_HELP_STRING([--enable-test],[enable test programs [[default=yes]]]))

AC_MSG_CHECKING([for a library suffix to use])
AC_ARG_WITH(library-suffix,
	[  --with-library-suffix=LIBSUFFIX    tag a suffix to the library names [[default=]]],
	[LIBSUFFIX="${withval}"],
	[LIBSUFFIX=])
AC_SUBST(LIBSUFFIX)


# Configurations for Some Specific Systems

if test "$TARGET" != ""; then
	case "$TARGET" in
# Cray XT3 Cross Compiler
	  cray_xt3_cross) CC="cc"
                    F77="ftn"
                    FC="ftn"
                    MPICC="$CC"
                    MPIF77="$F77"
                    MPIFC="$FC"
                    MPIRUN="yod"
                    CFLAGS="-O3 -B -fastsse -tp k8-64"
                    FFLAGS="-O3 -fastsse -tp k8-64 -Mpreprocess"
                    FCFLAGS="-O3 -fastsse -tp k8-64 -Mpreprocess"
                    FCLDFLAGS="-Mnomain"
                    enable_mpi="yes"
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
                    ax_f77_mangling="lower case, underscore, no extra underscore"
	  ;;
# Fujitsu Compiler for PRIMERGY
	  fujitsu_pg) CC="fcc"
                    F77="frt"
                    FC="frt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex -w"
                    FFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCFLAGS="-Kfast,ocl,preex -Cpp -Am -fs"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    cross_compiling="yes"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-KOMP"
                    AMDEFS="-PG"
	  ;;
# Fujitsu Compiler for PRIMEQUEST
	  fujitsu_pq) CC="fcc"
                    F77="frt"
                    FC="frt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex -w"
                    FFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCFLAGS="-Kfast,ocl,preex -Cpp -Am -fs"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-KOMP"
	  ;;
# Fujitsu Technical Computing Suite (for FX and PRIMERGY)
          fujitsu_tcs) CC="fcc"
                    F77="frt"
                    FC="frt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex -w"
                    FFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-Kopenmp"
          ;;
# Fujitsu Cross Compiler for FX10
          fujitsu_fx10_cross) CC="fccpx"
                    F77="frtpx"
                    FC="frtpx"
                    MPICC="mpifccpx"
                    MPIF77="mpifrtpx"
                    MPIFC="mpifrtpx"
                    CFLAGS="-Kfast,ocl,preex -w"
                    FFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCLDFLAGS="-mlcmain=main"
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
		    enable_sse2="no"
		    target="s64fx"
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-Kopenmp"
          ;;
# Fujitsu Compiler for K 
          K_shared) CC="fcc"
                    F77="frt"
                    FC="frt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex,noeval -w"
                    FFLAGS="-Kfast,ocl,preex,noeval -Cpp -fs -KPIC"
                    FCFLAGS="-Kfast,ocl,preex,noeval -Cpp -fs -KPIC"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    OMPFLAG="-Kopenmp"
          ;;
# Fujitsu Compiler for K 
          K_mpi_shared) CC="fcc"
                    F77="mpifrt"
                    FC="mpifrt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex,noeval -w"
                    FFLAGS="-Kfast,ocl,preex,noeval -Cpp -fs -KPIC"
                    FCFLAGS="-Kfast,ocl,preex,noeval -Cpp -fs -KPIC"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-Kopenmp"
          ;;
# Fujitsu Technical Computing Suite (for PRIMEHPC FX1000/FX700)
          fujitsu_arm) CC="fcc"
                    F77="frt"
                    FC="frt"
                    MPICC="mpifcc"
                    MPIF77="mpifrt"
                    MPIFC="mpifrt"
                    CFLAGS="-Kfast,ocl,preex -w"
                    FFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCFLAGS="-Kfast,ocl,preex -Cpp -fs"
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="fujitsu"
                    ax_cv_f77_compiler_vendor="fujitsu"
                    MPIRUN="mpiexec"
                    MPINP="-n"
                    OMPFLAG="-Kopenmp -Nlibomp"		    
          ;;
# Hitachi Optimizing Compiler for SR16000
          hitachi_sr16k) CC="cc"
                    F77="f90"
                    FC="f90"
                    MPICC="mpicc"
                    MPIF77="mpif90"
                    MPIFC="mpif90"
                    CFLAGS="-Os -64 -nolimit"
                    FFLAGS="-Oss -64 -nolimit -cpp"
                    FCFLAGS="-Oss -64 -nolimit -cpp -c"
                    FCLDFLAGS="-lf90s"
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="hitachi"
                    ax_cv_f77_compiler_vendor="hitachi"
                    MPIRUN="mpirun"
                    MPINP="-np"
                    OMPFLAG="-parallel -omp"
          ;;
# Hitachi Optimizing Compiler for HA8000
          hitachi_ha8k) CC="cc"
                    F77="f90"
                    FC="f90"
                    MPICC="mpicc"
                    MPIF77="mpif90"
                    MPIFC="mpif90"
                    CFLAGS="-Os -64 -nolimit"
                    FFLAGS="-Oss -64 -nolimit -cpp"
                    FCFLAGS="-Oss -64 -nolimit -cpp -c"
                    FCLDFLAGS="-lf90s"
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, underscore, no extra underscore"
                    ax_cv_c_compiler_vendor="hitachi"
                    ax_cv_f77_compiler_vendor="hitachi"
                    MPIRUN="mpirun"
                    MPINP="-np"
                    OMPFLAG="-parallel -omp"
          ;;
# IBM XL Compiler for Blue Gene Solution
	  ibm_bgl_cross) CC="blrts_xlc"
                    F77="blrts_xlf90"
                    FC="blrts_xlf90"
                    MPICC="$CC"
                    MPIF77="$F77"
                    MPIFC="$FC"
                    CFLAGS="-O3 -qarch=440d -qtune=440 -qstrict -I/bgl/BlueLight/ppcfloor/bglsys/include"
                    FFLAGS="-O3 -qarch=440d -qtune=440 -qsuffix=cpp=F -qfixed=72 -w -I/bgl/BlueLight/ppcfloor/bglsys/include"
                    FCFLAGS="-O3 -qarch=440d -qtune=440 -qsuffix=cpp=F90 -w -I/bgl/BlueLight/ppcfloor/bglsys/include"
                    MPILIBS="-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts"
                    enable_mpi="yes"
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
                    ax_cv_f77_compiler_vendor="ibm"
                    ax_cv_c_compiler_vendor="ibm"
                    ax_f77_mangling="lower case, no underscore, no extra underscore"
	  ;;
# NEC Compiler for SX-6i
	  nec_sx6i) CC="c++"
	            F77="f90"
                    FC="f90"
                    RANLIB="true"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
	  ;;
# NEC Cross Compiler for SX-6i
	  nec_sx6i_cross) CC="sxc++"
	            F77="sxf90"
                    FC="sxf90"
                    AR="sxar"
                    RANLIB="true"
                    STRIP="sxstrip"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
	  ;;
# NEC Cross Compiler for Earth Simulator
	  nec_es_cross) CC="esmpic++"
	            F77="esmpif90"
                    FC="esmpif90"
                    MPICC="$CC"
                    MPIF77="$F77"
                    MPIFC="$FC"
                    AR="esar"
                    RANLIB="true"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    LDFLAGS=""
                    enable_mpi="yes"
                    enable_omp="yes"
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
                    OMPFLAG="-Popenmp"
	  ;;
# NEC Compiler for SX-9
          nec_sx9)  CC="sxmpic++"
                    F77="sxmpif90"
                    FC="sxmpif90"
                    MPICC="$CC"
                    MPIF77="$F77"
                    MPIFC="$FC"
                    AR="ar"
                    RANLIB="true"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    LDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
                    OMPFLAG="-Popenmp"
          ;;
# NEC Cross Compiler for SX-9
          nec_sx9_cross) CC="sxmpic++"
                    F77="sxmpif90"
                    FC="sxmpif90"
                    MPICC="$CC"
                    MPIF77="$F77"
                    MPIFC="$FC"
                    AR="sxar"
                    RANLIB="true"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    LDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    cross_compiling="yes"
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
                    OMPFLAG="-Popenmp"
          ;;
# NEC Compiler for SX-Aurora
          nec_sxaurora)  CC="ncc"
                    F77="nfort"
                    FC="nfort"
                    MPICC="mpincc"
                    MPIF77="mpinfort"
                    MPIFC="mpinfort"
                    AR="ar"
                    RANLIB="true"
                    CFLAGS=""
                    FFLAGS=""
                    FCFLAGS=""
                    FCLDFLAGS=""
                    LDFLAGS=""
                    ac_cv_sizeof_void_p=8
                    ax_f77_mangling="lower case, no underscore, extra underscore"
                    ax_vector_machine="yes"
                    OMPFLAG="-fopenmp"
          ;;
# Intel Compiler for MIC
          intel_mic) CC="icc"
                    F77="ifort"
                    FC="ifort"
                    MPICC="mpiicc"
                    MPIF77="mpiifort"
                    MPIFC="mpiifort"
                    CFLAGS="-mmic"
                    FFLAGS="-mmic"
                    FCFLAGS="-mmic"
                    FCLDFLAGS="-mmic"
                    LDFLAGS="-mmic"		    
                    ax_cv_c_compiler_vendor="intel"
                    ax_cv_f77_compiler_vendor="intel"
	            OMPFLAG="-openmp"
          ;;
# Intel Cross Compiler for MIC
          intel_mic_cross) CC="icc"
                    F77="ifort"
                    FC="ifort"
                    MPICC="mpiicc"
                    MPIF77="mpiifort"
                    MPIFC="mpiifort"
                    CFLAGS="-mmic"
                    FFLAGS="-mmic"
                    FCFLAGS="-mmic"
                    FCLDFLAGS="-mmic"
                    LDFLAGS="-mmic"		    
                    cross_compiling="yes"
		    host="x86_64-pc-linux-gnu"
		    host_alias="x86_64-linux-gnu"
		    host_cpu="x86_64"   
		    host_os="linux-gnu"
		    host_vendor="pc"
		    target="k1om-mpss-linux-gnu"
		    target_alias="k1om-mpss-linux"
		    target_cpu="k1om"
		    target_os="linux-gnu"
		    target_vendor="mpss"
                    ax_cv_c_compiler_vendor="intel"
                    ax_cv_f77_compiler_vendor="intel"
		    OMPFLAG="-openmp"
          ;;
	  *)  echo 'TARGET' $TARGET 'cannot be specified.'
              exit 1
          ;;
	esac
fi

# GNU Compiler for Apple Mac OS X
if test "$build_vendor" = "apple"; then
   if test -z "$CC"; then
      CC="gcc"
   fi
   if test -z "$FC"; then   
      FC="gfortran"
   fi
   if test "$enable_32bit" = "yes"; then
      CFLAGS="$CFLAGS -m32"
      FCFLAGS="$FCFLAGS -m32"
   else
      CFLAGS="$CFLAGS -m64"
      FCFLAGS="$FCFLAGS -m64"
   fi
   if test "$enable_omp" = "yes"; then   
      if $CC --version 2>&1 | grep clang > /dev/null ; then
      	 OMPFLAG="-Xclang -fopenmp"
       	 LDFLAGS="$LDFLAGS -lomp"
      fi	 
   fi
fi

enable_dependency_tracking="no"


# Configurations for C and Fortran Compilers

########## Check C Compilers ##########
cflags_save="$CFLAGS"
AC_PROG_CC([icc icx ecc xlc_r xlc fcc pgcc pathcc cc gcc clang])
CFLAGS="$cflags_save"
if test "$TARGET" = ""; then
  if test "$CC" = "fcc"; then
    CFLAGS="-O3 -w"
    FFLAGS="-O3 -Cpp -fs"
    FCFLAGS="-O3 -Cpp -Am -fs"
    ax_f77_mangling="lower case, underscore, no extra underscore"
    MPIRUN="mpiexec"
    MPINP="-n"
  fi
fi

if test "$TARGET" = ""; then
	if test "$enable_mpi" = "yes"; then
		AC_CHECK_PROGS(MPICC, $MPICC mpicc mpifcc hcc mpcc_r mpcc mpxlc cmpicc blrts_xlc esmpic++, $CC)
		save_cc="$CC"
		CC="$MPICC"
		if test x = x"$MPILIBS"; then
			AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])
		fi
		if test x = x"$MPILIBS"; then
			AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
		fi
		if test x = x"$MPILIBS"; then
			AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
		fi
		if test x != x"$MPILIBS"; then
			AC_MSG_CHECKING([for mpi.h])
			AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]], [[]])],[AC_MSG_RESULT(yes)],[MPILIBS="" 
			AC_MSG_RESULT(no)])
	        else
	        	enable_mpi="not found"
		fi
		CC="$save_cc"
		CC=$MPICC
	fi

########## Check C Flags ########## 
	for ventest in intel:__ICC,__ICX,__ECC,__INTEL_COMPILER ibm:__xlc__,__xlC__,__IBMC__,__IBMCPP__ gnu:__GNUC__ sun:__SUNPRO_C,__SUNPRO_CC hp:__HP_cc,__HP_aCC dec:__DECC,__DECCXX,__DECC_VER,__DECCXX_VER borland:__BORLANDC__,__TURBOC__ comeau:__COMO__ cray:_CRAYC kai:__KCC lcc:__LCC__ metrowerks:__MWERKS__ sgi:__sgi,sgi microsoft:_MSC_VER watcom:__WATCOMC__ portland:__PGI pathscale:__PATHCC__,__PATHSCALE__; do
	    vencpp="defined("`echo $ventest | cut -d: -f2 | sed 's/,/) || defined(/g'`")"
	    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
	#if !($vencpp)
	      thisisanerror;
	#endif
	])], [ax_cv_]_AC_LANG_ABBREV[_compiler_vendor=`echo $ventest | cut -d: -f1`; break])
	done
	AC_MSG_CHECKING([for C compiler vendor])
	AC_MSG_RESULT($[ax_cv_]_AC_LANG_ABBREV[_compiler_vendor])
	  case $ax_cv_c_compiler_vendor in
	    dec) CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host $CFLAGS" ;;
	    gnu) CFLAGS="-O3 -fomit-frame-pointer $CFLAGS" ;;
	    intel) CFLAGS="-ansi_alias $CFLAGS" ;;	    
	    portland) CFLAGS="-O3 -B -fastsse $CFLAGS" ;;
	    pathscale) CFLAGS="-O3 $CFLAGS" ;;
	    ibm) CFLAGS="-O3 -qansialias -w -qarch=auto -qtune=auto $CFLAGS"
                 AR="ar -X any"
		 ;;
	    sun) case `(uname -m) 2>/dev/null` in
	           "i86pc") CFLAGS="-xtarget=native -xO5 -dalign $CFLAGS"
                            LDFLAGS="-xtarget=native -xO5"  
			    ;;      
	                 *) CFLAGS="-xtarget=native64 -xO5 -dalign $CFLAGS"
                            LDFLAGS="-xtarget=native64 -xO5" 
                            ;;			                     
          	 esac
                 ;;		 	    
	    hp)  CFLAGS="+Oall +Optrs_ansi +DSnative $CFLAGS" ;;
	  esac

	AC_DEFUN([_AX_PROG_FC_V_OUTPUT],
	[AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

	# Compile and link our simple test program by passing a flag (argument
	# 1 to this macro) to the Fortran compiler in order to get
	# "verbose" output that we can then parse for the Fortran linker
	# flags.
	ac_save_CFLAGS=$CFLAGS
	CFLAGS="$CFLAGS m4_default([$1], [$ac_cv_prog_[]_AC_LANG_ABBREV[]_v])"
	(eval echo $as_me:__oline__: \"$ac_link\") >&AS_MESSAGE_LOG_FD
	ac_[]_AC_LANG_ABBREV[]_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 | grep -v 'Driving:'`
	echo "$ac_[]_AC_LANG_ABBREV[]_v_output" >&AS_MESSAGE_LOG_FD
	CFLAGS=$ac_save_CFLAGS

	rm -f conftest*

	# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
	# /foo, /bar, and /baz are search directories for the Fortran linker.
	# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
	ac_[]_AC_LANG_ABBREV[]_v_output="`echo $ac_[]_AC_LANG_ABBREV[]_v_output |        grep 'LPATH is:' |        sed 's,.*LPATH is\(: *[[^ ]]*\).*,\1,;s,: */, -L/,g'` $ac_[]_AC_LANG_ABBREV[]_v_output"

	case $ac_[]_AC_LANG_ABBREV[]_v_output in
	  # If we are using xlf then replace all the commas with spaces.
	  *xlcentry*)
	    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/,/ /g'` ;;

	  # With Intel ifc, ignore the quoted -mGLOB_options_string stuff (quoted
	  # $LIBS confuse us, and the libraries appear later in the output anyway).
	  *mGLOB_options_string*)
	    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/\"-mGLOB[[^\"]]*\"/ /g'` ;;

	esac
	])

	AC_DEFUN([_AX_PROG_FC_V],
	[AC_CACHE_CHECK([how to get verbose linking output from C],
	                [ac_cv_prog_[]_AC_LANG_ABBREV[]_v],
	[AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
	[ac_cv_prog_[]_AC_LANG_ABBREV[]_v=
	# Try some options frequently used verbose output
	for ac_verb in -v -verbose --verbose -V -\#\#\#; do
	  _AX_PROG_FC_V_OUTPUT($ac_verb)
	  # look for -l* and *.a constructs in the output
	  for ac_arg in $ac_[]_AC_LANG_ABBREV[]_v_output; do
	     case $ac_arg in
	        [[\\/]]*.a | ?:[[\\/]]*.a | -[[lLRu]]*)
	          ac_cv_prog_[]_AC_LANG_ABBREV[]_v=$ac_verb
	          break 2 ;;
	     esac
	  done
	done
	if test -z "$ac_cv_prog_[]_AC_LANG_ABBREV[]_v"; then
	   AC_MSG_WARN([cannot determine how to obtain linking information from C])
	fi],
	                  [AC_MSG_WARN([compilation failed])])
	])])


	_AX_PROG_FC_V
	AC_CACHE_CHECK([for C libraries], ac_cv_[]_AC_LANG_ABBREV[]_libs,
	[if test "x$[]_AC_LANG_PREFIX[]LIBS" != "x"; then
	  ac_cv_[]_AC_LANG_ABBREV[]_libs="$[]_AC_LANG_PREFIX[]LIBS" # Let the user override the test.
	else

	_AX_PROG_FC_V_OUTPUT

	ac_cv_[]_AC_LANG_ABBREV[]_libs=

	# Save positional arguments (if any)
	ac_save_positional="$[@]"

	set X $ac_[]_AC_LANG_ABBREV[]_v_output
	while test $[@%:@] != 1; do
	  shift
	  ac_arg=$[1]
	  case $ac_arg in
	        [[\\/]]*.a | ?:[[\\/]]*.a)
	          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	              ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
	          ;;
	        -bI:*)
	          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	             [_AC_LINKER_OPTION([$ac_arg], ac_cv_[]_AC_LANG_ABBREV[]_libs)])
	          ;;
	          # Ignore these flags.
	        -lang* | -lcrt[[01]].o | -lcrtbegin.o | -lc | -lgcc | -libmil | -LANG:=*)
	          ;;
	        -lkernel32)
	          test x"$CYGWIN" != xyes && ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg"
	          ;;
	        -[[LRuY]])
	          # These flags, when seen by themselves, take an argument.
	          # We remove the space between option and argument and re-iterate
	          # unless we find an empty arg or a new option (starting with -)
	          case $[2] in
	             "" | -*);;
	             *)
	                ac_arg="$ac_arg$[2]"
	                shift; shift
	                set X $ac_arg "$[@]"
	                ;;
	          esac
	          ;;
	        -YP,*)
	          for ac_j in `echo $ac_arg | sed -e 's/-YP,/-L/;s/:/ -L/g'`; do
	            _AC_LIST_MEMBER_IF($ac_j, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	                               [ac_arg="$ac_arg $ac_j"
	                               ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_j"])
	          done
	          ;;
	        -[[lLR]]*)
	          _AC_LIST_MEMBER_IF($ac_arg, $ac_cv_[]_AC_LANG_ABBREV[]_libs, ,
	                             ac_cv_[]_AC_LANG_ABBREV[]_libs="$ac_cv_[]_AC_LANG_ABBREV[]_libs $ac_arg")
	          ;;
	          # Ignore everything else.
	  esac
	done
	# restore positional arguments
	set X $ac_save_positional; shift

	# We only consider "LD_RUN_PATH" on Solaris systems.  If this is seen,
	# then we insist that the "run path" must be an absolute path (i.e. it
	# must begin with a "/").
	case `(uname -sr) 2>/dev/null` in
	   "SunOS 5"*)
	      ac_ld_run_path=`echo $ac_[]_AC_LANG_ABBREV[]_v_output |
	                        sed -n 's,^.*LD_RUN_PATH *= *\(/[[^ ]]*\).*$,-R\1,p'`
	      test "x$ac_ld_run_path" != x &&
	        _AC_LINKER_OPTION([$ac_ld_run_path], ac_cv_[]_AC_LANG_ABBREV[]_libs)
	      ;;
	esac
	fi # test "x$[]_AC_LANG_PREFIX[]LIBS" = "x"
	])
	[]_AC_LANG_PREFIX[]LIBS="$ac_cv_[]_AC_LANG_ABBREV[]_libs"
	AC_SUBST([]_AC_LANG_PREFIX[]LIBS)
fi

########## Check Fortran Compilers ##########
if test "$enable_fortran" = "yes" || test "$enable_f90" = "yes" || test "$enable_saamg" = "yes"; then
  if test -z "$F77"; then
    F77=$FC
  fi
  AC_LANG_PUSH(Fortran 77)
  fflags_save="$FFLAGS"
  AC_PROG_F77([ifort ifx ifc efc xlf95_r xlf90_r xlf_r xlf95 xlf90 xlf ftn frt pgf95 pgf90 pathf95 pathf90 epcf90 f95 fort lf95 gfortran f90 g77])
  FFLAGS="$fflags_save"
  if test "$TARGET" = ""; then
    if test "$enable_mpi" = "yes"; then
      AC_CHECK_PROGS(MPIF77, mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpifrt cmpifc cmpif90c mpf90 mpif90 hf77 mpxlf_r mpxlf mpf77 mpif77, $F77)
      F77=$MPIF77
    fi
  fi
  AC_LANG_POP(Fortran 77)
  if test -n "$F77"; then
    FC=$F77
    MPIFC=$MPIF77
  fi
fi

if test "$enable_f90" = "yes" || test "$enable_saamg" = "yes"; then
  AC_LANG_PUSH(Fortran)
  fcflags_save="$FCFLAGS"
  AC_PROG_FC([ifort ifx ifc efc xlf95_r xlf90_r xlf_r xlf95 xlf90 xlf ftn frt pgf95 pgf90 pathf95 pathf90 epcf90 f95 fort lf95 gfortran f90])
  FCFLAGS="$fcflags_save"
  if test "$TARGET" = ""; then
    if test "$enable_mpi" = "yes"; then
      AC_CHECK_PROGS(MPIFC, mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpifrt cmpifc cmpif90c mpf90 mpif90 hf77, $F77)
      FC=$MPIFC
    fi
  fi
# AC_FC_SRCEXT(F90)
  AC_LANG_POP(Fortran)
fi

########## Check Fortran Flags ########## 
if test "$TARGET" = ""; then
  AC_LANG_PUSH(Fortran 77)
  ac_ext=F
  for ventest in intel:__INTEL_COMPILER ibm:__IBMC__ gnu:__GNUC__ sun:__SUNPRO_F77,__SUNPRO_F90,__SUNPRO_F95 portland:__PGI pathscale:__PATHCC__,__PATHSCALE__; do
    vencpp="defined("`echo $ventest | cut -d: -f2 | sed 's/,/) || defined(/g'`")"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
#if !($vencpp)
      thisisanerror;
#endif
])], [ax_cv_]_AC_LANG_ABBREV[_compiler_vendor=`echo $ventest | cut -d: -f1`; break])
  done
  AC_MSG_CHECKING([for Fortran compiler vendor])
  AC_MSG_RESULT($[ax_cv_]_AC_LANG_ABBREV[_compiler_vendor])
  AC_LANG_POP(Fortran 77)
  case $ax_cv_f77_compiler_vendor in
    gnu) FFLAGS="-O3 -fomit-frame-pointer $FFLAGS" 
         FCFLAGS="-O3 -fomit-frame-pointer $FCFLAGS" ;;
    intel) FFLAGS="-fpp $FFLAGS"
         FCFLAGS="-fpp $FCFLAGS"
#         LDFLAGS="-Vaxlib"
         FCLDFLAGS="-nofor_main" ;;
    portland) FFLAGS="-O3 -fastsse -Mpreprocess $FFLAGS" 
         FCFLAGS="-O3 -fastsse -Mpreprocess $FCFLAGS" 
         FCLDFLAGS="-Mnomain" ;;
    pathscale) FFLAGS="-O3 $FFLAGS" ;;
    ibm) FFLAGS="-O3 -qsuffix=cpp=F -qfixed=72 -w -qarch=auto -qtune=auto $FFLAGS"
         FCFLAGS="-O3 -qsuffix=cpp=F90 -w -qarch=auto -qtune=auto $FCFLAGS" ;;
    sun) case `(uname -m) 2>/dev/null` in
	   "i86pc") FFLAGS="-fpp -xtarget=native -xO5 -dalign $FFLAGS" 
		    FCFLAGS="-fpp -xtarget=native -xO5 -dalign $FCFLAGS"
                    FCLDFLAGS="-xO5" ;;  
	         *) FFLAGS="-fpp -xtarget=native64 -xO5 -dalign $FFLAGS" 
                    FCFLAGS="-fpp -xtarget=native64 -xO5 -dalign $FCFLAGS"
		    FCLDFLAGS="-xO5" ;;  
   	 esac
	 ;;
  esac

  if test x$ax_cv_c_compiler_vendor = x$ax_cv_f77_compiler_vendor; then
    CLIBS=""
  fi
fi

if test "$enable_fortran" = "yes" || test "$enable_f90" = "yes" || test "$enable_saamg" = "yes"; then
  AC_LANG_PUSH(Fortran 77)
  if test "$ax_cv_f77_compiler_vendor" = "sun"; then
    ax_f77_mangling="lower case, underscore, no extra underscore"
  fi
  if test -z "$ax_f77_mangling"; then
    AC_F77_WRAPPERS
    AC_F77_DUMMY_MAIN
  else
    case $ax_f77_mangling in
      "lower case, no underscore, no extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [name])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name])
      ;;
      "lower case, no underscore, extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [name])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## _])
      ;;
      "lower case, underscore, no extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [name ## _])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## _])
      ;;
      "lower case, underscore, extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [name ## _])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [name ## __])
      ;;
      "upper case, no underscore, no extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [NAME])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME])
      ;;
      "upper case, no underscore, extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [NAME])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## _])
      ;;
      "upper case, underscore, no extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [NAME ## _])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## _])
      ;;
      "upper case, underscore, extra underscore")
        AC_DEFINE(_AC_FC[_FUNC(name,NAME)], [NAME ## _])
        AC_DEFINE(_AC_FC[_FUNC_(name,NAME)], [NAME ## __])
      ;;
    esac
  fi
  AC_LANG_POP(Fortran 77)
fi


# Configurations for Build Options

######### Enable Fortran Interface ##########
if test "$enable_fortran" = "yes"; then
  AC_DEFINE([USE_FORTRAN], [1], [Define to 1 to enable Fortran subroutines.])
else
  enable_fortran="no"
  if test "$enable_f90" = "yes"; then
    AC_DEFINE([USE_FORTRAN], [1], [Define to 1 to enable Fortran subroutines.])
    enable_fortran="yes"
  else
    enable_f90="no"
  fi
fi
if test "$enable_fortran" = "yes"; then
  if test "$CC" = "fcc"; then
    AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
  else
    if test "$CC" = "mpifcc"; then
      AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
    fi
  fi
fi
 
######### Enable SA-AMG Preconditioner ##########
if test "$enable_saamg" = "yes"; then
  AC_DEFINE([USE_FORTRAN], [1], [Define to 1 to enable Fortran subroutines.])
  AC_DEFINE([USE_SAAMG], [1], [Define to 1 to enable SA-AMG preconditioner.])
  case $ax_cv_f77_compiler_vendor in
    ibm) FCFLAGS="$FCFLAGS -WF,-DZERO_ORIGIN=1" ;;
    gnu) FCFLAGS="$FCFLAGS -Wp,-DZERO_ORIGIN=1" ;;
    hitachi) FCFLAGS="$FCFLAGS" ;;
    *)   FCFLAGS="$FCFLAGS -DZERO_ORIGIN=1" ;;
  esac
else
  enable_saamg="no"
fi
if test "$enable_fortran" = "yes"; then
  if test "$CC" = "fcc"; then
    AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
  else
    if test "$CC" = "mpifcc"; then
      AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
    fi
  fi
else
  if test "$enable_saamg" = "yes"; then
    if test "$CC" = "fcc"; then
      AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
    else 
      if test "$CC" = "mpifcc"; then
        AC_DEFINE([USE_MAIN__], [1], [Use MAIN__ macro with Fortran.])
      fi
    fi
  fi
fi

########## Enable MPI and OpenMP Directives ##########
if test "$TARGET" = ""; then
  if test -z "$enable_mpi"; then
    enable_mpi="no"
    MPIRUN=""
  else
    if test -z "$MPIRUN"; then	
      AC_CHECK_PROGS(MPIRUN, mpirun mpiexec, mpirun)
    fi
    if test -z "$MPINP"; then
      if test "$MPIRUN" = "mpiexec"; then 
        MPINP="-n"
      else		
        MPINP="-np"
      fi	
    fi
  fi
  if test "$enable_omp" = "yes"; then
    save_cflags="$CFLAGS"
    openmp_flags="-openmp -qopenmp -qsmp=omp -xopenmp -fiopenmp -fopenmp -KOMP -mp -omp -Popenmp none"
    for ompf in $openmp_flags; do
      AC_MSG_CHECKING([for OpenMP flag $ompf])
      case $ompf in
        none) CFLAGS=$save_cflags ;;
           *) CFLAGS="$save_cflags $ompf" ;;
      esac
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <omp.h>]], [[
        #pragma omp parallel
        omp_set_num_threads(1)
      ]])],[OMPFLAG=$ompf;    AC_MSG_RESULT(yes); break],[])
      AC_MSG_RESULT(no)
    done
    CFLAGS=$save_cflags
    if test "$OMPFLAG" = ""; then
      enable_omp="not found"
    else
      if test "$OMPFLAG" = "none"; then
        OMPFLAG=""
      fi
    fi
  else
    case $ax_cv_c_compiler_vendor in
        hitachi) OMPFLAG="-noparallel" ;;
    esac
    case $ax_cv_f77_compiler_vendor in
        hitachi) OMPFLAG="-noparallel" ;;
    esac
    enable_omp="no"

  fi
fi

if test "$enable_mpi" = "yes"; then
  CC=$MPICC
  F77=$MPIF77
  FC=$MPIFC
  CFLAGS="$CFLAGS -DUSE_MPI"
  case $ax_cv_f77_compiler_vendor in
    ibm) FCFLAGS="$FCFLAGS -WF,-DUSE_MPI" 
         FFLAGS="$FFLAGS -WF,-DUSE_MPI"
         ;;
    gnu) FCFLAGS="$FCFLAGS -Wp,-DUSE_MPI"
         FFLAGS="$FFLAGS -Wp,-DUSE_MPI"
         ;;
    *)   FCFLAGS="$FCFLAGS -DUSE_MPI"
         FFLAGS="$FFLAGS -DUSE_MPI"
         ;;
  esac
else
  enable_mpi="no"
fi

if test "$enable_omp" = "yes"; then
  if test "$OMPFLAG" != ""; then
    CFLAGS="$CFLAGS $OMPFLAG"
    FFLAGS="$FFLAGS $OMPFLAG"
    FCFLAGS="$FCFLAGS $OMPFLAG"
    LDFLAGS="$LDFLAGS $OMPFLAG"
  fi
fi

########## Enable x87 FPUs ##########
  AC_MSG_CHECKING([for x87 FPU])
  case $target in
    i*86-* | x86*) AC_MSG_RESULT(yes)
         AC_DEFINE([HAS_X87_FPU], [1], [Define to 1 to enable x87 FPU.])
         ;;
    *)   AC_MSG_RESULT(no) ;;
  esac
#  enable_sse2="no"

if test "$enable_quad" = "yes"; then

########## Enable Intel Streaming SIMD Extensions ##########
#  enable_sse2="yes"
if test "$enable_sse2" = "yes"; then
  AC_MSG_CHECKING([for SSE2])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <emmintrin.h>]], [[
    __m128d xs;
    double hi;
    xs  = _mm_set_pd(1.0,1.0);
    _mm_store_sd(&hi,xs);
  ]])],[AC_MSG_RESULT(yes)
  AC_DEFINE([USE_SSE2], [1], [Define to 1 to enable quad precision operations using SSE2.])
  AC_DEFINE([USE_FMA2_SSE2], [1], [Define to 1 to enable quad precision operations using FMA and SSE2.])
  ],[enable_sse2="no"; AC_MSG_RESULT(no)])
fi
  if test "$TARGET" = ""; then
    if test "$enable_sse2" != "yes"; then
      save_cflags="$CFLAGS"
      case $ax_cv_c_compiler_vendor in
        intel) quad_flags="-mp" ;;
        ibm)   quad_flags="-qstrict none" ;;
        sun)   quad_flags="-fsimple=0 none" ;;
        *)   quad_flags="none" ;;
      esac
      if test "$MPFLAG" != ""; then
        quad_flags="$MPFLAG $quad_flags"
      fi
      for mpf in $quad_flags; do
        case $ax_cv_c_compiler_vendor in
          sun) mpf="$mpf -fns=no" ;;
        esac
        case $mpf in
          none) CFLAGS=$save_cflags ;;
             *) CFLAGS="$save_cflags $mpf" ;;
        esac
        AC_MSG_CHECKING([for mp flag $mpf])
        AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>],[[double r,e;
  		double a=1.0;
		double b=1.1e-16;
		r=a+b;
		e=b-(r-a);
                printf("r=%e e=%e ",r,e);
		if(e==0.0) return 1;
        ]])],
        [MPFLAG=$mpf;    AC_MSG_RESULT(yes); break],
        AC_MSG_RESULT(no))
      done
      CFLAGS=$save_cflags
      if test "$MPFLAG" = ""; then
        enable_quad="not found"
      else
        if test "$MPFLAG" = "none"; then
          MPFLAG=""
        fi
      fi
    fi
  fi
  if test "$enable_quad" = "yes"; then
    AC_DEFINE([USE_QUAD_PRECISION], [1], [Define to 1 to enable quad precision operations.])
    AC_DEFINE([USE_FAST_QUAD_ADD], [1], [Define to 1 to enable quad precision operations using fast quad add.])
    if test "$enable_sse2" = "yes"; then
	  case $ax_cv_c_compiler_vendor in
	    gnu) MPFLAG="-msse2" ;;
	  esac
    fi
  fi
else
  enable_quad="no"
fi

########## Enable SIMD Optimization for Intel Fused Multiply and Add ##########
if test "$enable_fma" = "yes"; then
  AC_DEFINE([HAS_FMA], [1], [Define to 1 to enable quad precision operations using FMA.])
  case $ax_cv_c_compiler_vendor in
    intel) MPFLAG="$MPFLAG -fma" ;;
  esac
else
  enable_fma="no"
fi

########## Enable Vectorization for NEC SX Series ##########
if test "$ax_vector_machine" = "yes"; then
  AC_DEFINE([USE_VEC_COMP], [1], [Define to 1 to use vectorization.])
else
  ax_vector_machine="no"
fi

########## Enable Debugging ###########
if test "$enable_debug" = "yes"; then
  CFLAGS="$CFLAGS -g -D_DEBUG"
  if test "$enable_fortran" = "yes"; then
    FFLAGS="$FFLAGS -g -DDEBUG"
    case $ax_cv_f77_compiler_vendor in
      ibm) FCFLAGS="$FCFLAGS -g -WF,-DDEBUG" ;;
      gnu) FCFLAGS="$FCFLAGS -g -Wp,-DDEBUG" ;;
      *)   FCFLAGS="$FCFLAGS -g -DDEBUG" ;;
    esac
  fi
else
  enable_debug="no"
fi

########## Enable Profiling ##########
if test "$enable_gprof" = "yes"; then
  case $ax_cv_c_compiler_vendor in
    sun)     CFLAGS="$CFLAGS -xpg" 
             FFLAGS="$FFLAGS -xpg"
             FCFLAGS="$FCFLAGS -xpg"
             LDFLAGS="$LDFLAGS -xpg"
             ;;
    fujitsu) CFLAGS="$CFLAGS -p"
             FFLAGS="$FFLAGS -p"
             FCFLAGS="$FCFLAGS -p"
             LDFLAGS="$LDFLAGS -p"
             ;;
    hp)      CFLAGS="$CFLAGS -p"
             FFLAGS="$FFLAGS -p"
             FCFLAGS="$FCFLAGS -p"
             LDFLAGS="$LDFLAGS -p"
             ;;
    gnu)     CFLAGS="$CFLAGS -fno-omit-frame-pointer -pg"
             FFLAGS="$FFLAGS -fno-omit-frame-pointer -pg"
             FCFLAGS="$FCFLAGS -fno-omit-frame-pointer -pg"
             LDFLAGS="$LDFLAGS -fno-omit-frame-pointer -pg"
             ;;
    *)       CFLAGS="$CFLAGS -pg"
             FFLAGS="$FFLAGS -pg"
             FCFLAGS="$FCFLAGS -pg"
             LDFLAGS="$LDFLAGS -pg"
             ;;
  esac
else
  enable_gprof="no"
fi

########## Enable 64bit Binary Executables ##########
if test "$enable_64bit" = "yes"; then
  case $ax_cv_c_compiler_vendor in
    sun)     CFLAGS="$CFLAGS -m64" 
             FFLAGS="$FFLAGS -m64"
             FCFLAGS="$FCFLAGS -m64"
             LDFLAGS="$LDFLAGS -m64"
             ;;
    ibm)     CFLAGS="$CFLAGS -q64" 
             FFLAGS="$FFLAGS -q64"
             FCFLAGS="$FCFLAGS -q64"
             LDFLAGS="$LDFLAGS -q64"
             ;;
  esac
else
  enable_64bit="no"
fi

########## Enable 64bit Integer ##########
if test "$enable_longlong" = "yes"; then

  case $ax_cv_c_compiler_vendor in
    sun)     CFLAGS="$CFLAGS -m64" 
             FFLAGS="$FFLAGS -m64"
             FCFLAGS="$FCFLAGS -m64"
             LDFLAGS="$LDFLAGS -m64"
             ;;
    ibm)     CFLAGS="$CFLAGS -q64" 
             FFLAGS="$FFLAGS -q64"
             FCFLAGS="$FCFLAGS -q64"
             LDFLAGS="$LDFLAGS -q64"
             ;;
  esac
  CFLAGS="$CFLAGS -D_LONG__LONG"

  case $ax_cv_f77_compiler_vendor in
    intel) FFLAGS="$FFLAGS -DLONG__LONG -integer-size 64"
         FCFLAGS="$FCFLAGS -DLONG__LONG -integer-size 64"
	 ;;
    ibm) FFLAGS="$FFLAGS -WF,-DLONG__LONG -qintsize=8"
         FCFLAGS="$FCFLAGS -WF,-DLONG__LONG -qintsize=8"
	 ;;
    fujitsu) FFLAGS="$FFLAGS -DLONG__LONG -CcdII8"
         FCFLAGS="$FCFLAGS -DLONG__LONG -CcdII8"
         ;;
    gnu) FFLAGS="$FFLAGS -Wp,-DLONG__LONG -fdefault-integer-8"
         FCFLAGS="$FCFLAGS -Wp,-DLONG__LONG -fdefault-integer-8"
         ;;
    *)   FFLAGS="$FFLAGS -DLONG__LONG" 
	 FCFLAGS="$FCFLAGS -DLONG__LONG" 
	 ;;
  esac
else
  enable_longlong="no"
fi

########## Enable Long Double (Quadruple) Precision Operations ##########
if test "$enable_longdouble" = "yes"; then

  case $ax_cv_c_compiler_vendor in
    ibm)     CFLAGS="$CFLAGS -qldbl128"
             ;;
  esac
  CFLAGS="$CFLAGS -D_LONG__DOUBLE"

  case $ax_cv_f77_compiler_vendor in
    intel) FFLAGS="$FFLAGS -DLONG__DOUBLE -real-size 128"
         FCFLAGS="$FCFLAGS -DLONG__DOUBLE -real-size 128"
	 ;;
    *)   FFLAGS="$FFLAGS -DLONG__DOUBLE" 
	 FCFLAGS="$FCFLAGS -DLONG__DOUBLE" 
	 ;;
  esac
else
  enable_longdouble="no"
fi


########## Enable Complex Scalar Operations ##########
AC_CHECK_HEADERS([complex.h])
if test "$enable_complex" = "yes"; then
  AC_DEFINE([USE_COMPLEX], [1], [Define to 1 to enable complex scalar operations.])
  if test -z "$HAVE_COMPLEX_H"; then
    CFLAGS="$CFLAGS -D_COMPLEX"
    FFLAGS="$FFLAGS -DCOMPLEX" 
    FCFLAGS="$FCFLAGS -DCOMPLEX" 
  else
    AC_MSG_WARN([no complex.h])
    enable_complex="no"
  fi
else
  enable_complex="no"
fi


########## Enable Quadmath Operations ##########
AC_CHECK_HEADERS([quadmath.h])
if test "$enable_longdouble" = "yes"; then
  if test -z "$HAVE_QUADMATH_H"; then
    CFLAGS="$CFLAGS -lquadmath"
    FFLAGS="$FFLAGS -lquadmath" 
    FCFLAGS="$FCFLAGS -lquadmath"
    AC_MSG_WARN([GCC libquadmath is used if long double is enabled])
  fi
fi


# Configurations for Environment-Dependent Data Sizes

########## Check Size of int ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(int)
else
  AC_DEFINE([SIZEOF_INT], [$ac_cv_sizeof_int], [Test the size of a `int', as computed by sizeof.])
fi

########## Check Size of long int ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(long)
else
  AC_DEFINE([SIZEOF_LONG], [$ac_cv_sizeof_long], [Test the size of a `long', as computed by sizeof.])
fi

########## Check Size of long long int ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(long long)
else
  AC_DEFINE([SIZEOF_LONG_LONG], [$ac_cv_sizeof_long_long], [Test the size of a `long long', as computed by sizeof.])
fi

########## Check Size of float ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(float)
else
  AC_DEFINE([SIZEOF_FLOAT], [$ac_cv_sizeof_float], [Test the size of a `float', as computed by sizeof.])
fi

########## Check Size of double ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(double)
else
  AC_DEFINE([SIZEOF_DOUBLE], [$ac_cv_sizeof_double], [Test the size of a `double', as computed by sizeof.])
fi

########## Check Size of long double ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(long double)
else
  AC_DEFINE([SIZEOF_LONG_DOUBLE], [$ac_cv_sizeof_long_double], [Test the size of a `long double', as computed by sizeof.])
fi

########## Check Size of __float128 ##########
if test -z "$HAVE_QUADMATH_H"; then
  if test "$TARGET" = ""; then
    AC_CHECK_SIZEOF(__float128)
  else
    AC_DEFINE([SIZEOF_LONG_DOUBLE], [$ac_cv_sizeof___float128], [Test the size of a `__float128', as computed by sizeof.])
  fi
fi

########## Check Size of size_t ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(size_t)
else
  AC_DEFINE([SIZEOF_SIZE_T], [$ac_cv_sizeof_size_t], [Test the size of a `size_t', as computed by sizeof.])
fi

########## Check Size of void_p ##########
if test "$TARGET" = ""; then
  AC_CHECK_SIZEOF(void *)
else
  AC_DEFINE([SIZEOF_VOID_P], [$ac_cv_sizeof_void_p], [Test the size of a `void *', as computed by sizeof.])
fi


# Configurations for Autotools

AC_DISABLE_SHARED
LT_INIT
AC_PROG_INSTALL

AM_CONDITIONAL(ENABLE_FORTRAN, test x$enable_fortran = xyes)
AM_CONDITIONAL(ENABLE_F90, test x$enable_f90 = xyes)
if test "$enable_f90" = "yes" ; then
  AM_CONDITIONAL(ENABLE_FORTRAN, test x$enable_f90 = xyes)
fi
AM_CONDITIONAL(ENABLE_SAAMG, test x$enable_saamg = xyes)
AM_CONDITIONAL(ENABLE_MPI, test x"$enable_mpi" = xyes)
AM_CONDITIONAL(ENABLE_OMP, test x"$enable_omp" = xyes)
AM_CONDITIONAL(ENABLE_QUAD, test x$enable_quad = xyes)
AM_CONDITIONAL(ENABLE_SSE2, test x$enable_sse2 = xyes)
AM_CONDITIONAL(DISABLE_TEST, test x$enable_test = xno)
AC_SUBST(enable_mpi)
AC_SUBST(enable_omp)
AC_SUBST(enable_saamg)
AC_SUBST(enable_quad)
AC_SUBST(enable_longdouble)
AC_SUBST(enable_complex)
AC_SUBST(enable_fortran)
AC_SUBST(enable_f90)
AC_SUBST(MPIRUN)
AC_SUBST(MPINP)
AC_SUBST(AMDEFS)
AC_SUBST(SIZEOF_VOID_P, $ac_cv_sizeof_void_p)
AC_SUBST(FCLDFLAGS)


# Configurations for Some Header Files

m4_warn([obsolete],
[The preprocessor macro 'STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  C89 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS([malloc.h stdlib.h string.h sys/time.h])


# Configurations for Compiler Characteristics and Type Definitions

AC_C_CONST
AC_TYPE_SIZE_T


# Configurations for Library Functions

LIBS="-lm $LIBS"


# Output

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/vector/Makefile
                 src/matrix/Makefile
                 src/matvec/Makefile
                 src/array/Makefile
                 src/solver/Makefile
		 src/esolver/Makefile
                 src/precon/Makefile
                 src/system/Makefile
                 src/precision/Makefile
                 src/fortran/Makefile
                 src/fortran/amg/Makefile
                 include/Makefile
                 test/Makefile
                 include/lisf.h
                 ])

DEFS="$DEFSADD $DEFS"
AC_OUTPUT

echo
echo "  Build with OpenMP library                  = $enable_omp"
echo "  Build with MPI library                     = $enable_mpi"
echo "  Enable FORTRAN 77 compatible interface     = $enable_fortran"
echo "  Enable Fortran 90 compatible interface     = $enable_f90"
echo "  Enable SA-AMG preconditioner               = $enable_saamg"
echo "  Enable double-double precision support     = $enable_quad"
if test "$enable_quad" = "yes"; then
echo "  Enable SSE2                                = $enable_sse2"
echo "  Enable FMA                                 = $enable_fma"
fi
echo "  Enable long double precision support       = $enable_longdouble"
echo "  Enable 64bit integer support               = $enable_longlong"
echo "  Enable complex scalar support              = $enable_complex"
if test "$enable_debug" = "yes"; then
echo "  Enable debugging                           = $enable_debug"
fi
echo "  Enable dynamic linking                     = $enable_shared"
echo "  Enable profiling                           = $enable_gprof"
echo
echo "  C compiler        = $CC"
echo "  C flags           = $CFLAGS $DEFS $MPFLAG $AMDEFS"
echo "  C libraries       = $LIBS $CLIBS"

if test "$enable_fortran" = "yes"; then
    echo "  F77 compiler      = $F77"
    echo "  F77 flags         = $FFLAGS"
    echo "  F77 libraries     = $FLIBS"
fi
if test "$enable_f90" = "yes" || test "$enable_saamg" = "yes"; then
    echo "  F90 compiler      = $FC"
    echo "  F90 flags         = $FCFLAGS"
    echo "  F90 libraries     = $FCLIBS"
fi
if test "$enable_mpi" = "yes"; then
    echo "  MPI libraries     = $MPILIBS"
    echo "  MPIRUN script     = $MPIRUN"
fi
echo
