@echo off
setlocal
title Leonard Analysis

set "BIN=%~dp0leonard.exe"

if not exist "%BIN%" (
    echo Error: leonard.exe not found in %~dp0
    pause
    exit /b 1
)

"%BIN%" %*
set "ERR=%ERRORLEVEL%"

if not "%ERR%"=="0" (
    echo(
    echo -------------------------------------------------------------
    echo Leonard encountered an error (exit code: %ERR%).
    pause
)
exit /b %ERR%
