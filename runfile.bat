@echo off
setlocal enabledelayedexpansion

:: Ensure repo is up to date
:: git stash
:: git pull

:: Move to scripts directory
cd %~dp0scripts

:: Load input variables
call input_batch.bat

:: Backup files
md %backup%
robocopy /s /e %in% %backup%

:: Process ppmooring with python script
%pythonenv% %script% %in%

:: Move input files
robocopy /move /s /e %in% %Level0%

:: Push changes to remote repository
:: git add --all
:: git commit -m "Auto Upload"
:: git push