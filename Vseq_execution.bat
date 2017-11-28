@echo off
title Vseq .exe
color 30
echo ============================================================================
echo =                                                                          =
echo =                        Vseq Program                                      =
echo =                                                                          =
echo ============================================================================
  echo.                     
echo.
echo Cargando... 
set PROGRAM=%cd%\Programa_R\Vseq_init.R
echo %PROGRAM%
::"C:\Program Files\R\R-3.4.2\bin\i386\Rscript.exe" C:\Users\Proteomica VI\Desktop\Vseq Enrique\Programa_R\Vseq_init.R

"C:\Program Files\R\R-3.4.2\bin\x64\Rscript.exe" "C:\Users\Proteomica VI\Desktop\Vseq_Enrique\Programa_R\Vseq_init.R"
echo Finalizado. Presiona cualquier tecla para salir
::echo %cd%\Programa_R
pause
exit
