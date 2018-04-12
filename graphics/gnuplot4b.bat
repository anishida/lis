:start

@rem Sample batch script to draw distribution of Ritz values on complex plane
@rem using gnuplot
@rem Make executable files enabling complex arithmetic and
@rem run the following command:
@rem
@rem > gnuplot4b.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set bindir=..\win
@set srcdir=..\test
@%bindir%\etest5b.exe %srcdir%\testmat3.mtx rvalues.mtx -e ai -ss 20


@rem Draw distribution of Ritz values.

@set filename=rvalues.mtx
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%" gnuplot4b.plt

@del rvalues.mtx

