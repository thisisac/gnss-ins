@echo off
cd /d %~dp0\..
echo Run GNSSINS pseudo-online-fast demo with first 50 GNSS epochs.
"src\GNSSINS\x64\Debug\GNSSINS.exe" "output\demo_output" 50 pseudo-online-fast
pause
