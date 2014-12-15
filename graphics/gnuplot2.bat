:start

@rem Sample batch script to draw matrix pattern using gnuplot
@rem Run the following command:
@rem
@rem > gnuplot2.bat matrix_filename
@rem


@rem Set path to the utility. 

@set path="c:\Program Files (x86)\gnuplot\bin";%PATH%


@rem Draw matrix pattern.

@set filename=%1
@for /f "tokens=1-5" %%a in ('findstr /v %% "%filename%"') do @(
    set size=%%a
    set nnz=%%c
    goto :break
    )
:break

@gnuplot.exe -e "filename='%filename%'; size=%size%; nnz=%nnz%" gnuplot2.plt









