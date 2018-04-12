:start

@rem Sample batch script to draw eigenvectors using gnuplot
@rem Make executable files and run the following command:
@rem
@rem > gnuplot5.bat
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Run test programs.

@set subspace=3
@set bindir=..\win
@set srcdir=..\test
@%bindir%\etest5.exe %srcdir%\testmat.mtx nul evectors.mtx nul nul -e si -ie ii -ss %subspace%


@rem Draw eigenvectors.

@set filename=evectors.mtx
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%; subspace=%subspace%" gnuplot5.plt

@del evectors.mtx

