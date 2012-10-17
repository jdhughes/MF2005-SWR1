rem ***Sample 01
cd SWRSample01\
mf2005-swr_x64.exe SWRSample01.01min.nam
mf2005-swr_x64.exe SWRSample01.nam
cd ..\Python\
python SWRSample01.py

rem ***Sample 02
cd ..\SWRSample02\
mf2005-swr_x64.exe SWRSample02.nam
cd ..\Python\
python SWRSample02.py

rem ***Sample 03
cd ..\SWRSample03\
mf2005-swr_x64.exe SWRSample03.nam
cd ..\Python\
python SWRSample03.py

rem ***Sample 04
cd ..\SWRSample04\
mf2005-swr_x64.exe SWRSample04.nam
cd ..\Python\
python SWRSample04.py

rem ***Sample 05
cd ..\SWRSample05\
mf2005-swr_x64.exe SWRSample05.nam
cd ..\Python\
python SWRSample05.py

rem END OF BATCH FILE
cd ..\

pause
