@echo off
rem   Configuration Script for Microsoft Build Engine

:start
rem   Default Configurations

@echo # Configurations for Microsoft Build Engine > Makefile
@echo. >> Makefile
@echo. 
@set prefix=
@set msvcrt=0
@set notest=0
@set omp=0
@set intelc=0
@set icllvm=0
@set fortran=0
@set msmpi32=0
@set msmpi64=0
@set impi=0
@set saamg=0
@set longlong=0
@set longdouble=0
@set complex=0
@set debug=0
@set noifpu=0
@set usercflags=
@set userfflags=
@set userldflags=

:again
rem   Build Options

@if "%1" == "-h" goto usage
@if "%1" == "--help" goto usage
@if "%1" == "--prefix" goto setprefix
@if "%1" == "--msvcrt" goto setmsvcrt
@if "%1" == "--enable-omp" goto setomp
@if "%1" == "--enable-msmpi32" goto setmsmpi32
@if "%1" == "--enable-msmpi64" goto setmsmpi64
@if "%1" == "--enable-impi" goto setimpi
@if "%1" == "--enable-intelc" goto setintelc
@if "%1" == "--enable-icllvm" goto seticllvm
@if "%1" == "--enable-fortran" goto setfortran
@if "%1" == "--enable-saamg" goto setsaamg
@if "%1" == "--enable-longlong" goto setlonglong
@if "%1" == "--enable-longdouble" goto setlongdouble
@if "%1" == "--enable-complex" goto setcomplex
@if "%1" == "--enable-debug" goto setdebug
@if "%1" == "--disable-ifpu" goto setnoifpu
@if "%1" == "--disable-test" goto setnotest
@if "%1" == "--cflags" goto setusercflags
@if "%1" == "--fflags" goto setuserfflags
@if "%1" == "--ldflags" goto setuserldflags
@if "%1" == "" goto genmakefile

:usage

@echo Usage: configure.bat [options]
@echo.	
@echo Options:
@echo.	
@echo.	--prefix PREFIX		Install Lis in directory PREFIX
@echo.	--msvcrt		Use msvcrt*.dll instead of static crt
@echo.	--enable-omp		Build with OpenMP library
@echo.	--enable-msmpi32	Build with 32bit Microsoft MPI library
@echo.	--enable-msmpi64	Build with 64bit Microsoft MPI library
@echo.	--enable-impi		Build with Intel MPI library
@echo.	--enable-intelc		Use Intel C Compiler
@echo.	--enable-icllvm		Use Intel LLVM C Compiler
@echo.	--enable-fortran	Enable Fortran interface
@echo.	--enable-saamg		Enable SA-AMG preconditioner
@echo.	--enable-longlong	Enable 64bit integer support
@echo.	--enable-longdouble	Enable long double support
@echo.	--enable-complex	Complex scalar support (for Intel Compilers)
@echo.	--enable-debug		Enable debugging
@echo.	--disable-ifpu		Disable Intel FPU support
@echo.	--disable-test		Disable building test programs
@echo.	--cflags FLAG		Pass FLAG to C compiler
@echo.	--fflags FLAG		Pass FLAG to Fortran compiler
@echo.	--ldflags FLAG		Pass FLAG to linker
@echo.	
@goto end

:setprefix

@shift
@set prefix=%1
@shift
@goto again

:setmsvcrt

@set msvcrt=1
@shift
@goto again

:setnotest

@set notest=1
@shift
@goto again

:setomp

@set omp=1
@shift
@goto again

:setmsmpi32

@set msmpi32=1
@shift
@goto again

:setmsmpi64

@set msmpi64=1
@shift
@goto again

:setimpi

@set impi=1
@set intelc=1
@shift
@goto again

:setintelc

@set intelc=1
@shift
@goto again

:seticllvm

@set icllvm=1
@set intelc=1
@shift
@goto again

:setfortran

@set fortran=1
@shift
@goto again

:setsaamg

@set fortran=1
@set saamg=1
@shift
@goto again

:setlonglong

@set longlong=1
@shift
@goto again

:setlongdouble

@set longdouble=1
@shift
@goto again

:setcomplex

@set complex=1
@shift
@goto again

:setdebug

@set debug=1
@shift
@goto again

:setnoifpu

@set noifpu=1
@shift
@goto again

:setusercflags

@shift
@set usercflags=%~1
@shift
@goto again

:setuserfflags

@shift
@set userfflags=%~1
@shift
@goto again

:setuserldflags

@shift
@set userldflags=%~1
@shift
@goto again

:genmakefile

@echo # Installation directory >> Makefile
@if not "(%prefix%)" == "()" (
@echo PREFIX=%prefix% >> Makefile
) else (
@echo PREFIX=.. >> Makefile
)
@echo. >> Makefile

@if (%msvcrt%) == (1) (
@echo # Use msvcrt*.dll instead of static runtime library >> Makefile
@echo msvcrt=1 >> Makefile
@echo. >> Makefile
@echo.	Use msvcrt*.dll instead of static crt	= yes
) else (
@echo.	Use msvcrt*.dll instead of static crt	= no
)

@if (%notest%) == (1) (
@echo # Disable building test programs >> Makefile
@echo notest=1 >> Makefile
@echo. >> Makefile
@echo.	Build test programs			= no
) else (
@echo.	Build test programs			= yes
)

@if (%omp%) == (1) (
@echo # Build with OpenMP library >> Makefile
@echo omp=1 >> Makefile
@echo. >> Makefile
@echo.	Build with OpenMP library		= yes
) else (
@echo.	Build with OpenMP library		= no
)

@if (%msmpi32%) == (1) (
@echo # Build with 32bit Microsoft MPI library >> Makefile
@echo msmpi32=1 >> Makefile
@echo. >> Makefile
@echo.	Build with 32bit Microsoft MPI library	= yes
) else (
@echo.	Build with 32bit Microsoft MPI library	= no
)

@if (%msmpi64%) == (1) (
@echo # Build with 64bit Microsoft MPI library >> Makefile
@echo msmpi64=1 >> Makefile
@echo. >> Makefile
@echo.	Build with 64bit Microsoft MPI library	= yes
) else (
@echo.	Build with 64bit Microsoft MPI library	= no
)

@if (%impi%) == (1) (
@echo # Build with Intel MPI library >> Makefile
@echo impi=1 >> Makefile
@echo. >> Makefile
@echo.	Build with Intel MPI library		= yes
) else (
@echo.	Build with Intel MPI library		= no
)

@if (%intelc%) == (1) (
@echo # Use Intel C Compiler >> Makefile
@echo intelc=1 >> Makefile
@echo. >> Makefile
@echo.	Use Intel C Compiler			= yes
) else (
@echo.	Use Intel C Compiler			= no
)

@if (%icllvm%) == (1) (
@echo # Use Intel LLVM C Compiler >> Makefile
@echo icllvm=1 >> Makefile
@echo. >> Makefile
@echo.	Use Intel LLVM C Compiler		= yes
) else (
@echo.	Use Intel LLVM C Compiler		= no
)

@if (%fortran%) == (1) (
@echo # Enable Fortran interface >> Makefile
@echo fortran=1 >> Makefile
@echo. >> Makefile
@echo.	Enable Fortran interface		= yes
) else (
@echo.	Enable Fortran interface		= no
)

@if (%saamg%) == (1) (
@echo # Enable SA-AMG preconditioner >> Makefile
@echo saamg=1 >> Makefile
@echo. >> Makefile
@echo.	Enable SA-AMG preconditioner		= yes
) else (
@echo.	Enable SA-AMG preconditioner		= no
)

@if (%longlong%) == (1) (
@echo # Enable 64bit integer support >> Makefile
@echo longlong=1 >> Makefile
@echo. >> Makefile
@echo.	Enable 64bit integer support		= yes
) else (
@echo.	Enable 64bit integer support		= no
)

@if (%longdouble%) == (1) (
@echo # Enable long double support >> Makefile
@echo longdouble=1 >> Makefile
@echo. >> Makefile
@echo.	Enable long double support		= yes
) else (
@echo.	Enable long double support		= no
)

@if (%complex%) == (1) if (%intelc%) == (1) (
@echo # Enable complex scalar support >> Makefile
@echo complex=1 >> Makefile
@echo. >> Makefile
@echo.	Enable complex scalar support		= yes
) else (
@echo.	Enable complex scalar support		= no
)

@if (%debug%) == (1) (
@echo. >> Makefile
@echo.	Enable debug mode			= yes
) else (
@echo # Disable Debugging >> Makefile
@echo NODEBUG=1 >> Makefile
@echo. >> Makefile
@echo.	Enable debug mode			= no
)

@if (%noifpu%) == (1) (
@echo # Disable Intel FPU support >> Makefile
@echo noifpu=1 >> Makefile
@echo. >> Makefile
@echo.	Disable Intel FPU support		= yes
) else (
@echo.	Disable Intel FPU support		= no
)

@echo # User-defined C flags >> Makefile
@if not "(%usercflags%)" == "()" (
@echo USER_CFLAGS = %usercflags% >> Makefile
) else (
@echo USER_CFLAGS = >> Makefile
)
@echo. >> Makefile

@echo # User-defined Fortran flags >> Makefile
@if not "(%userfflags%)" == "()" (
@echo USER_FFLAGS = %userfflags% >> Makefile
) else (
@echo USER_FFLAGS = >> Makefile
)
@echo. >> Makefile

@echo # User-defined linker flags >> Makefile
@if not "(%userldflags%)" == "()" (
@echo USER_LDFLAGS = %userldflags% >> Makefile
) else (
@echo USER_LDFLAGS = >> Makefile
)
@echo. >> Makefile
) 

@echo. >> Makefile
@type Makefile.in >> Makefile
@echo.

:end
