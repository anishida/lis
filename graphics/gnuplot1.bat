:start

@rem Sample batch script to draw residual histories using gnuplot
@rem Make executable files and run the following command:
@rem
@rem > gnuplot1.bat
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

@gnuplot.exe gnuplot1.plt

@del ILU0-* 

