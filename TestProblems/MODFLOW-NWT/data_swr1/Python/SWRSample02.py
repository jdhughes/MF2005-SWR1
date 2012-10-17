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
#	print x
	xc = np.chararray(len(x), itemsize=16)
	for i in range(0,len(x)):
		text = fmt % ( x[i] )
		xc[i] = string.strip(string.ljust(text,16))
#	print xc
	a.set_xticklabels(xc, size=sz)
	success = 1
	return success

def set_sizeyaxis(a,fmt,sz):
	success = 0
	y = a.get_yticks()
#	print y
	yc = np.chararray(len(y), itemsize=16)
	for i in range(0,len(y)):
		text = fmt % ( y[i] )
		yc[i] = string.strip(string.ljust(text,16))
#	print yc
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

mpl.rcParams['font.sans-serif']          = 'Arial'
mpl.rcParams['font.serif']               = 'Times'
mpl.rcParams['font.cursive']             = 'Zapf Chancery'
mpl.rcParams['font.fantasy']             = 'Comic Sans MS'
mpl.rcParams['font.monospace']           = 'Courier New'
mpl.rcParams['pdf.compression']          = 0
mpl.rcParams['pdf.fonttype']             = 42

ticksize = 10
mpl.rcParams['legend.fontsize']  = 7
mpl.rcParams['axes.labelsize']   = 12
mpl.rcParams['xtick.labelsize']  = ticksize
mpl.rcParams['ytick.labelsize']  = ticksize

fig = mpl.pyplot.figure(figsize=(4, 4), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.125,right=0.975,bottom=0.10,top=0.975)
plta = fig.add_subplot(111)
plta.plot(swift2d_time, swift2d_s57, 'o', linewidth=0, markeredgecolor='b', markerfacecolor='w', markersize=4, markevery=24, label='SWIFT2D observation location 1')
plta.plot(swift2d_time, swift2d_s61, 'o', linewidth=0, markeredgecolor='r', markerfacecolor='w', markersize=4, markevery=24, label='SWIFT2D observation location 2')
plta.plot(swift2d_time, swift2d_s64, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=4, markevery=24, label='SWIFT2D observation location 3')
plta.plot(swr_t, swr_s57, 'b', linewidth=1, label='SWR1 observation location 1')
plta.plot(swr_t, swr_s61, 'r', linewidth=1, label='SWR1 observation location 2')
plta.plot(swr_t, swr_s64, 'k', linewidth=1, label='SWR1 observation location 3')
plta.set_xlabel('Time, hours')	
plta.set_ylabel('Stage, m')
plta.set_xlim(0,30)
plta.set_ylim(1.2,2.40001)
i = set_sizexaxis(plta,'%16.6g',ticksize)
i = set_sizeyaxis(plta,'%16.6g',ticksize)
#-legend for BH and 10 min
#leg = plta.legend(('SWIFT2D Observation Location 1', 'SWIFT2D Observation Location 2', 'SWIFT2D Observation Location 3', 'SWR1 Observation Location 1', 'SWR1 Observation Location 2', 'SWR1 Observation Location 3'), loc='upper center',ncol=2,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
handles, labels = plta.get_legend_handles_labels()
leg = plta.legend((handles[0],handles[3],handles[1],handles[4],handles[2],handles[5]),
                  (labels[0],labels[3],labels[1],labels[4],labels[2],labels[5]),
                  loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

#--save figures
outfigpng = '..\Figures\SWR01Sample02.png'
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
	show()
else:
	print 'no screen output requested...'

   