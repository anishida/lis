:start

@rem Sample batch script to draw residual histories using gnuplot
@rem Make executable files and run the following command:
@rem
@rem > gnuplot1b.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set bindir=..\win
@set srcdir=..\test
@%bindir%\test1.exe %srcdir%\testmat.mtx 1 nul ILU0-BiCG -i bicg -p ilu
@%bindir%\test1.exe %srcdir%\testmat.mtx 1 nul ILU0-CGS -i cgs -p ilu
@%bindir%\test1.exe %srcdir%\testmat.mtx 1 nul ILU0-BiCGSTAB -i bicgstab -p ilu


@rem Draw residual histories using gnuplot.

@gnuplot.exe -e "filenames='ILU0-BiCG ILU0-CGS ILU0-BiCGSTAB'" gnuplot1b.plt

@del ILU0-* 

