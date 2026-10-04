@echo off
setlocal

echo ============================================================
echo   Building TMG USER1 DLL - NX 11 / TMG
echo   Intel Fortran 2016 x64
echo ============================================================
echo.

rem ============================================================
rem Configuration
rem ============================================================

set "TMGOPENDLL=tmgopen_user.dll"

set "SDRC_TMG=C:\Program Files\Siemens\NX 11.0\NXCAE_EXTRAS\tmg"

set "IFORTVARS=C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2016.1.146\windows\bin\ifortvars.bat"

rem ============================================================
rem Check TMG installation
rem ============================================================

if not exist "%SDRC_TMG%\if\libtmgopen_ivf.lib" goto TMG_ERROR

if not exist "%SDRC_TMG%\exe" goto TMG_ERROR

echo TMG:
echo   %SDRC_TMG%
echo.

rem ============================================================
rem Initialize Intel Fortran x64 environment
rem ============================================================

if not exist "%IFORTVARS%" goto IFORTVARS_ERROR

echo Initializing Intel Fortran x64 environment...
call "%IFORTVARS%" intel64

if errorlevel 1 goto IFORT_ENV_ERROR

echo.

rem ============================================================
rem Verify compiler and linker
rem ============================================================

where ifort
if errorlevel 1 goto IFORT_ERROR

where lib
if errorlevel 1 goto LIB_ERROR

where link
if errorlevel 1 goto LINK_ERROR

echo.
echo Intel Fortran environment OK.
echo.

rem ============================================================
rem Compiler / linker
rem ============================================================

set "FCOMPILER=ifort"
set "LDCMD=link"
set "LIBCMD=lib"

set "FFLAGS=/nologo /DMAYA_WINNT /f77rtl /Qsave /Qzero /Z7 /check:nobounds /traceback"
set "LDFLAGS=/nologo /dll"

rem ============================================================
rem Optional USERLIB
rem ============================================================

set "USERLIB="

if exist "%SDRC_TMG%\exe\USERLIB.obj" (
    set "USERLIB=%SDRC_TMG%\exe\USERLIB.obj"
)

if exist "%SDRC_TMG%\exe\USERLIB.o" (
    set "USERLIB=%SDRC_TMG%\exe\USERLIB.o"
)

rem ============================================================
rem Clean old temporary files
rem ============================================================

del /q TMGDataInterface.obj 2>nul
del /q User1Function.obj 2>nul
del /q TMGFuncInterface.obj 2>nul
del /q TMGSubrInterface.obj 2>nul
del /q USER1.obj 2>nul
del /q libtmgopen.lib 2>nul
del /q tmgopen.exp 2>nul
del /q tmgopen.lib 2>nul
del /q err 2>nul

rem ============================================================
rem Copy TMG interface library
rem ============================================================

echo Copying TMG interface library...

copy /Y "%SDRC_TMG%\if\libtmgopen_ivf.lib" "libtmgopen.lib"

if errorlevel 1 goto COPY_ERROR

echo.

rem ============================================================
rem Extract required interface objects
rem ============================================================

echo Extracting TMG interface objects...

%LIBCMD% /NOLOGO /EXTRACT:TMGDataInterface.obj libtmgopen.lib

if errorlevel 1 goto EXTRACT_ERROR

%LIBCMD% /NOLOGO /EXTRACT:User1Function.obj libtmgopen.lib

if errorlevel 1 goto EXTRACT_ERROR

%LIBCMD% /NOLOGO /EXTRACT:TMGFuncInterface.obj libtmgopen.lib

if errorlevel 1 goto EXTRACT_ERROR

%LIBCMD% /NOLOGO /EXTRACT:TMGSubrInterface.obj libtmgopen.lib

if errorlevel 1 goto EXTRACT_ERROR

echo.
echo TMG interface objects extracted successfully.
echo.

rem ============================================================
rem Check USER1 source
rem ============================================================

if not exist USER1.f goto USER1_SOURCE_ERROR

rem ============================================================
rem Compile USER1.f
rem ============================================================

echo ============================================================
echo Compiling USER1.f
echo ============================================================
echo.

echo %FCOMPILER% %FFLAGS% /c /W0 USER1.f /FoUSER1.obj

%FCOMPILER% %FFLAGS% /c /W0 USER1.f /FoUSER1.obj >compile_output.txt 2>&1

if errorlevel 1 goto COMPILE_ERROR

if not exist USER1.obj goto COMPILE_ERROR

echo.
echo USER1.f compiled successfully.
echo.

rem ============================================================
rem Optional JAXA library
rem ============================================================

set "TMGJAXA="

if exist "%SDRC_TMG%\exe\crunchTmg.lib" (
    set "TMGJAXA=%SDRC_TMG%\exe\crunchTmg.lib"
)

rem ============================================================
rem Link USER1 DLL
rem ============================================================

echo ============================================================
echo Linking %TMGOPENDLL%
echo ============================================================
echo.

%LDCMD% %LDFLAGS% ^
 /OUT:%TMGOPENDLL% ^
 TMGDataInterface.obj ^
 TMGFuncInterface.obj ^
 USER1.obj ^
 User1Function.obj ^
 TMGSubrInterface.obj ^
 %TMGJAXA% ^
 %USERLIB% ^
 Kernel32.lib >link_output.txt 2>&1

if errorlevel 1 goto LINK_ERROR2

if not exist "%TMGOPENDLL%" goto LINK_ERROR2

echo.
echo ============================================================
echo SUCCESS
echo ============================================================
echo.
echo Created:
echo   %CD%\%TMGOPENDLL%
echo.

dir "%TMGOPENDLL%"

goto CLEANUP_SUCCESS


rem ============================================================
rem ERROR HANDLERS
rem ============================================================

:TMG_ERROR
echo.
echo ERROR: TMG installation was not found.
echo Expected:
echo %SDRC_TMG%\if\libtmgopen_ivf.lib
goto FAILED

:IFORTVARS_ERROR
echo.
echo ERROR: Intel Fortran environment script was not found:
echo %IFORTVARS%
goto FAILED

:IFORT_ENV_ERROR
echo.
echo ERROR: Intel Fortran environment initialization failed.
goto FAILED

:IFORT_ERROR
echo.
echo ERROR: ifort is not available after initializing Intel environment.
goto FAILED

:LIB_ERROR
echo.
echo ERROR: Intel lib.exe is not available.
goto FAILED

:LINK_ERROR
echo.
echo ERROR: linker is not available.
goto FAILED

:COPY_ERROR
echo.
echo ERROR: Could not copy libtmgopen_ivf.lib.
goto FAILED

:EXTRACT_ERROR
echo.
echo ERROR: Could not extract required TMG interface object.
goto FAILED

:USER1_SOURCE_ERROR
echo.
echo ERROR: USER1.f was not found in:
echo %CD%
goto FAILED

:COMPILE_ERROR
echo.
echo ============================================================
echo USER1 COMPILATION FAILED
echo ============================================================
echo.
if exist compile_output.txt type compile_output.txt
goto FAILED

:LINK_ERROR2
echo.
echo ============================================================
echo DLL LINKING FAILED
echo ============================================================
echo.
if exist link_output.txt type link_output.txt
goto FAILED

:FAILED
echo.
echo ============================================================
echo BUILD FAILED
echo ============================================================
echo.
pause
exit /b 1


:CLEANUP_SUCCESS
echo.
echo Cleaning temporary build files...

del /q TMGDataInterface.obj 2>nul
del /q User1Function.obj 2>nul
del /q TMGFuncInterface.obj 2>nul
del /q TMGSubrInterface.obj 2>nul
del /q USER1.obj 2>nul
del /q libtmgopen.lib 2>nul
del /q tmgopen.exp 2>nul
del /q tmgopen.lib 2>nul

echo.
echo Build completed successfully.
echo.
pause
exit /b 0