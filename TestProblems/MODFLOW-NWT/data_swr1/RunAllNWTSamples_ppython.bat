rem ***Sample 01
cd SWRSample01\
mfnwt-SWR.exe SWRSample01.01min.nam
mfnwt-SWR.exe SWRSample01.nam
cd ..\Python\
python SWRSample01.py

rem ***Sample 02
cd ..\SWRSample02\
mfnwt-SWR.exe SWRSample02.nam
cd ..\Python\
python SWRSample02.py

rem ***Sample 03
cd ..\SWRSample03\
mfnwt-SWR.exe SWRSample03.nam
cd ..\Python\
python SWRSample03.py

rem ***Sample 04
cd ..\SWRSample04\
mfnwt-SWR.exe SWRSample04.nam
cd ..\Python\
python SWRSample04.py

rem ***Sample 05
cd ..\SWRSample05\
mfnwt-SWR.exe SWRSample05-nwt.nam
cd ..\Python\
python SWRSample05.py

rem END OF BATCH FILE
cd ..\

pause
