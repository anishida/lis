:start

@rem Sample batch script to draw eigenvalue distribution on complex plane
@rem using gnuplot
@rem Make executable files enabling complex arithmetic and
@rem run the following command:
@rem
@rem > gnuplot4.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set bindir=..\win
@set srcdir=..\test
@%bindir%\etest5.exe %srcdir%\testmat3.mtx evalues.mtx nul nul nul -e ai -ss 20


@rem Draw eigenvalue distribution.

@set filename=evalues.mtx
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%" gnuplot4.plt

@del evalues.mtx

