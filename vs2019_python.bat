@echo off
REM Step 1: Set up cl
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"

REM Step 2: Re-add Python (adjust to your install location)
set PATH=C:\Users\YourName\AppData\Local\Programs\Python\Python311\;%PATH%

REM Optional: If you use Conda instead
REM call C:\Users\YourName\miniconda3\condabin\conda.bat activate base

REM Launch VS Code
code .
