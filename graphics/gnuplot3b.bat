:start

@rem Sample batch script to draw Ritz value distribution on complex plane
@rem using gnuplot
@rem Make executable files enabling complex arithmetic and
@rem run the following command:
@rem
@rem > gnuplot3b.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set bindir=..\win
@set srcdir=..\test
@%bindir%\etest5b.exe %srcdir%\testmat2.mtx rvalues.mtx -e ai -ss 20


@rem Draw Ritz value distribution.

@set filename=rvalues.mtx
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%" gnuplot3b.plt

@del rvalues.mtx

