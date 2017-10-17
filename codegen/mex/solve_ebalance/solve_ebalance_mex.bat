@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2016b
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2016b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=solve_ebalance_mex
set MEX_NAME=solve_ebalance_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for solve_ebalance > solve_ebalance_mex.mki
echo COMPILER=%COMPILER%>> solve_ebalance_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> solve_ebalance_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> solve_ebalance_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> solve_ebalance_mex.mki
echo LINKER=%LINKER%>> solve_ebalance_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> solve_ebalance_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> solve_ebalance_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> solve_ebalance_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> solve_ebalance_mex.mki
echo BORLAND=%BORLAND%>> solve_ebalance_mex.mki
echo OMPFLAGS= >> solve_ebalance_mex.mki
echo OMPLINKFLAGS= >> solve_ebalance_mex.mki
echo EMC_COMPILER=mingw64>> solve_ebalance_mex.mki
echo EMC_CONFIG=optim>> solve_ebalance_mex.mki
"C:\Program Files\MATLAB\R2016b\bin\win64\gmake" -B -f solve_ebalance_mex.mk
