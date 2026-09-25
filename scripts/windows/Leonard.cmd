@echo off
setlocal enabledelayedexpansion
title Leonard Analysis

set "BIN_DIR=%~dp0"
set "BIN=%BIN_DIR%leonard.exe"
set "FILE=%~1"

if not exist "%BIN%" (
    echo Error: leonard.exe not found in %BIN_DIR%
    pause
    exit /b 1
)

:: If double-clicked without dropping a file, prompt with native Windows File Dialog
if "%FILE%"=="" (
    echo No workspace dropped. Opening file picker...
    for /f "usebackq delims=" %%I in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $f = New-Object System.Windows.Forms.OpenFileDialog; $f.Filter = 'FlowJo Workspace (*.wsp)|*.wsp|All Files (*.*)|*.*'; $f.Title = 'Select FlowJo Workspace'; if ($f.ShowDialog() -eq 'OK') { $f.FileName }"`) do (
        set "FILE=%%I"
    )
)

if "%FILE%"=="" (
    echo No workspace selected. Exiting.
    exit /b 0
)

:: Switch working directory to the folder containing the workspace
for %%F in ("%FILE%") do set "WORKDIR=%%~dpF"
cd /d "%WORKDIR%"

echo Running Leonard on: "%FILE%"
echo -------------------------------------------------------------
"%BIN%" "%FILE%"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo -------------------------------------------------------------
    echo Leonard encountered an error (exit code: %ERRORLEVEL%).
    pause
)
