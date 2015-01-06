@echo off
rem   Configuration script for Microsoft Build Engine

:start
rem   Default settings

@echo # Configuration for Microsoft Build Engine > Makefile
@echo. >> Makefile
@echo. 
@set prefix=
@set omp=0
@set intel=0
@set fortran=0
@set mpi=0
@set mpich=0
@set mpi64=0
@set saamg=0
@set longlong=0
@set longdouble=0
@set debug=0
@set noifpu=0
@set usercflags=
@set userfflags=
@set userldflags=

:again
rem   Handle arguments

@if "%1" == "-h" goto usage
@if "%1" == "--help" goto usage
@if "%1" == "--prefix" goto setprefix
@if "%1" == "--enable-omp" goto setomp
@if "%1" == "--enable-mpi" goto setmpi
@if "%1" == "--enable-mpich" goto setmpich
@if "%1" == "--enable-mpi64" goto setmpi64
@if "%1" == "--enable-intel" goto setintel
@if "%1" == "--enable-fortran" goto setfortran
@if "%1" == "--enable-saamg" goto setsaamg
@if "%1" == "--enable-longlong" goto setlonglong
@if "%1" == "--enable-longdouble" goto setlongdouble
@if "%1" == "--enable-debug" goto setdebug
@if "%1" == "--disable-ifpu" goto setnoifpu
@if "%1" == "--cflags" goto setusercflags
@if "%1" == "--fflags" goto setuserfflags
@if "%1" == "--ldflags" goto setuserldflags
@if "%1" == "" goto genmakefile

:usage

@echo Usage: configure [options]
@echo.	
@echo Options:
@echo.	
@echo.	--prefix PREFIX		Install Lis in directory PREFIX
@echo.	--enable-omp		Build with OpenMP library
@echo.	--enable-mpi		Build with Microsoft MPI library
@echo.	--enable-mpich		Build with MPICH library
@echo.	--enable-mpi64		Build with 64bit MPI library
@echo.	--enable-intel		Use Intel Compiler
@echo.	--enable-fortran	Enable Fortran interface
@echo.	--enable-saamg		Enable SA-AMG preconditioner
@echo.	--enable-longlong	Enable 64bit integer support
@echo.	--enable-longdouble	Enable long double support
@echo.	--enable-debug		Enable debugging
@echo.	--disable-ifpu		Disable Intel FPU support
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

:setomp

@set omp=1
@shift
@goto again

:setmpi

@set mpi=1
@shift
@goto again

:setmpich

@set mpich=1
@shift
@goto again

:setmpi64

@set mpi64=1
@shift
@goto again

:setintel

@set intel=1
@shift
@goto again

:setfortran

@set intel=1
@set fortran=1
@shift
@goto again

:setsaamg

@set intel=1
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
@set usercflags=%1
@shift
@goto again

:setuserfflags

@shift
@set userfflags=%1
@shift
@goto again

:setuserldflags

@shift
@set userldflags=%1
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

@if (%omp%) == (1) (
@echo # Build with OpenMP library >> Makefile
@echo omp=1 >> Makefile
@echo. >> Makefile
@echo.	Build with OpenMP library		= yes
) else (
@echo.	Build with OpenMP library		= no
)

@if (%mpi%) == (1) (
@echo # Build with Microsoft MPI library >> Makefile
@echo mpi=1 >> Makefile
@echo. >> Makefile
@echo.	Build with Microsoft MPI library	= yes
) else (
@echo.	Build with Microsoft MPI library	= no
)

@if (%mpich%) == (1) (
@echo # Build with MPICH library >> Makefile
@echo mpich=1 >> Makefile
@echo. >> Makefile
@echo.	Build with MPICH library		= yes
) else (
@echo.	Build with MPICH library		= no
)

@if (%mpi64%) == (1) (
@echo # Build with 64bit MPI library >> Makefile
@echo mpi64=1 >> Makefile
@echo. >> Makefile
@echo.	Build with 64bit MPI library		= yes
) else (
@echo.	Build with 64bit MPI library		= no
)

@if (%intel%) == (1) (
@echo # Use Intel Compiler >> Makefile
@echo intel=1 >> Makefile
@echo. >> Makefile
@echo.	Use Intel Compiler			= yes
) else (
@echo.	Use Intel Compiler			= no
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

@if (%debug%) == (1) (
@echo # Enable Debugging >> Makefile
@echo debug=1 >> Makefile
@echo. >> Makefile
@echo.	Enable debug mode			= yes
) else (
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
