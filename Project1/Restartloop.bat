@echo off
:while
 timeout /t 2
 tasklist /FI "IMAGENAME eq Project1.exe" 2>NUL | find /I /N "Project1.exe">NUL
 if NOT "%ERRORLEVEL%" == "0" start "" "C:\Users\hanju\Desktop\Research\NSphdiff_Mixture_Fraction_V\Project1\Project1\Project1.exe"
 goto :while