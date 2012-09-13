import sys
import string
import numpy as np
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

#import arrayUtil as au
import MFBinaryClass as mfb 
reload(mfb)


def set_sizexaxis(a,fmt,sz):
    success = 0
    x = a.get_xticks()
#    print x
    xc = np.chararray(len(x), itemsize=16)
    for i in range(0,len(x)):
        text = fmt % ( x[i] )
        xc[i] = string.strip(string.ljust(text,16))
#    print xc
    a.set_xticklabels(xc, size=sz)
    success = 1
    return success

def set_sizeyaxis(a,fmt,sz):
    success = 0
    y = a.get_yticks()
#    print y
    yc = np.chararray(len(y), itemsize=16)
    for i in range(0,len(y)):
        text = fmt % ( y[i] )
        yc[i] = string.strip(string.ljust(text,16))
#    print yc
    a.set_yticklabels(yc, size=sz)
    success = 1
    return success


#--main script
OutputToScreen = False
try:
    if sys.argv[1] == '-s':
        OutputToScreen = True
except:
    OutputToScreen = False


#--SWR results
SWR_file = '..\\SWRSample02\\Results\\SWRSample02_Stage.stg'
SWRObj = mfb.SWR_Record(0,SWR_file)
itime  = SWRObj.get_item_number('totim')
istage = SWRObj.get_item_number('stage')
ce1 = SWRObj.get_gage(57)
nt =  np.shape(ce1)[0]
swr_t = np.zeros(nt) 
swr_s57 = np.zeros(nt) 
for n  in range(0,nt):
    swr_t[n] = ce1[n,itime] / ( 60. * 60. )
    swr_s57[n] = ce1[n,istage]

SWRObj = mfb.SWR_Record(0,SWR_file)
ce1 = SWRObj.get_gage(61)
nt =  np.shape(ce1)[0]
swr_s61 = np.zeros(nt) 
for n  in range(0,nt):
    swr_s61[n] = ce1[n,istage]

SWRObj = mfb.SWR_Record(0,SWR_file)
ce1 = SWRObj.get_gage(64)
nt =  np.shape(ce1)[0]
swr_s64 = np.zeros(nt) 
for n  in range(0,nt):
    swr_s64[n] = ce1[n,istage]

#--SWIFT2D Results
infile = '..\\SWRSample02\\SWIFT2D\\Processed\\c3r6_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_time = np.zeros(rows, float)
swift2d_s57 = np.zeros(rows, float)
#--parse s57 data
for i in range(0,rows):
    ipos = i + 1
    swift2d_time[i] = float( swift2d[ipos,2] ) * 2.0 / 60.0
    swift2d_s57[i] = float( swift2d[ipos,3] )

infile = '..\\SWRSample02\\SWIFT2D\\Processed\\c7r6_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_s61 = np.zeros(rows, float)
#--parse s61 data
for i in range(0,rows):
    ipos = i + 1
    swift2d_s61[i] = float( swift2d[ipos,3] )

infile = '..\\SWRSample02\\SWIFT2D\\Processed\\c10r6_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_s64 = np.zeros(rows, float)
#--parse s64 data
for i in range(0,rows):
    ipos = i + 1
    swift2d_s64[i] = float( swift2d[ipos,3] )

mpl.rcParams['legend.fontsize'] = 6
fig = mpl.pyplot.figure(figsize=(3, 3), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.125,right=0.975,bottom=0.10,top=0.975)
plta = fig.add_subplot(111)
plta.plot(swift2d_time, swift2d_s57, 'o', linewidth=0, markeredgecolor='b', markerfacecolor='w', markersize=4, markevery=2)
plta.plot(swift2d_time, swift2d_s61, 'o', linewidth=0, markeredgecolor='r', markerfacecolor='w', markersize=4, markevery=2)
plta.plot(swift2d_time, swift2d_s64, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=4, markevery=2)
plta.plot(swr_t, swr_s57, 'b', linewidth=1)
plta.plot(swr_t, swr_s61, 'r', linewidth=1)
plta.plot(swr_t, swr_s64, 'k', linewidth=1)
plta.set_xlabel('Time,  hours', size='8')    
plta.set_ylabel(r'Stage,  $m$', size='8')
plta.set_xlim(0,4)
plta.set_ylim(1.2,2.40001)
i = set_sizexaxis(plta,'%16.6g','6')
i = set_sizeyaxis(plta,'%16.6g','6')
#-legend for BH and 10 min
leg = plta.legend(('SWIFT2D Obs. 1', 'SWIFT2D Obs. 2', 'SWIFT2D Obs. 3', 'SWR1 Obs. 1', 'SWR1 Obs. 2', 'SWR1 Obs. 3'), loc='upper center',ncol=2,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 
# matplotlib.text.Text instances
for l in leg.get_texts():
    l.set_fontsize(6)    # the legend text fontsize


outfigpng = '..\Figures\SWR01Sample02.png'
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
    show()
else:
    print 'no screen output requested...'

   