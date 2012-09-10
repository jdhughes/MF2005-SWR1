import sys
import string
import pylab
from pylab import *
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
#	a.set_xticklabels(xc, size=sz)
	a.set_xticklabels(xc)
	for tl in a.get_xticklabels():
	    tl.set_size(sz)
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
#	a.set_yticklabels(yc, size=sz)
	a.set_yticklabels(yc)
	for tl in a.get_yticklabels():
	    tl.set_size(sz)
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
SWR_file = '..\\SWRSample03\\Results\\SWRSample03_Stage.stg'
SWRObj = mfb.SWR_Record(0,SWR_file)
itime  = SWRObj.get_item_number('totim')
istage = SWRObj.get_item_number('stage')
ce1 = SWRObj.get_gage(39)
nt =  np.shape(ce1)[0]
swr_t = np.zeros(nt) 
swr_s01 = np.zeros(nt) 
swr_r01 = np.zeros(nt) 
for n  in range(0,nt):
	swr_t[n] = ce1[n,itime] / ( 60. * 60. )
	swr_s01[n] = ce1[n,istage]

SWRObj = mfb.SWR_Record(0,SWR_file)
ce1 = SWRObj.get_gage(72)
nt =  np.shape(ce1)[0]
swr_s02 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s02[n] = ce1[n,istage]

SWRObj = mfb.SWR_Record(0,SWR_file)
ce1 = SWRObj.get_gage(105)
nt =  np.shape(ce1)[0]
swr_s03 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s03[n] = ce1[n,istage]

#--SWIFT2D Results
infile = '..\\SWRSample03\\SWIFT2D\\Processed\\c5r7_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_time = zeros(rows, float)
swift2d_s01 = zeros(rows, float)
#--parse s01 data
for i in range(0,rows):
	ipos = i + 1
	swift2d_time[i] = float( swift2d[ipos,2] ) * 2.0 / 60.0
	swift2d_s01[i] = float( swift2d[ipos,3] )

infile = '..\\SWRSample03\\SWIFT2D\\Processed\\c8r7_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_s02 = zeros(rows, float)
#--parse s02 data
for i in range(0,rows):
	ipos = i + 1
	swift2d_s02[i] = float( swift2d[ipos,3] )

infile = '..\\SWRSample03\\SWIFT2D\\Processed\\c11r7_stage.csv'
swift2d = np.genfromtxt(infile,delimiter=',')
datasize = swift2d.shape
rows = datasize[0] - 1
cols = datasize[1]
swift2d_s03 = zeros(rows, float)
#--parse s03 data
for i in range(0,rows):
	ipos = i + 1
	swift2d_s03[i] = float( swift2d[ipos,3] )

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


fig = pylab.figure(figsize=(4, 4), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.125,right=0.975,bottom=0.10,top=0.975)
plta = fig.add_subplot(111)
plta.plot(swift2d_time, swift2d_s01, 'o', linewidth=0, markeredgecolor='b', markerfacecolor='w', markersize=4, markevery=12, label='SWIFT2D observation location 1')
plta.plot(swift2d_time, swift2d_s02, 'o', linewidth=0, markeredgecolor='r', markerfacecolor='w', markersize=4, markevery=12, label='SWIFT2D observation location 2')
plta.plot(swift2d_time, swift2d_s03, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=4, markevery=12, label='SWIFT2D observation location 3')
plta.plot(swr_t, swr_s01, 'b', linewidth=1, label='SWR1 observation location 1')
plta.plot(swr_t, swr_s02, 'r', linewidth=1, label='SWR1 observation location 2')
plta.plot(swr_t, swr_s03, 'k', linewidth=1, label='SWR1 observation location 3')
plta.set_xlabel('Time, hours')	
plta.set_ylabel('Stage, m')
plta.set_xlim(0,160)
plta.set_ylim(0,1.25)
i = set_sizexaxis(plta,'%16.6g',ticksize)
i = set_sizeyaxis(plta,'%16.6g',ticksize)
#-legend
#leg = plta.legend(('SWIFT2D Obs. 1', 'SWIFT2D Obs. 2', 'SWIFT2D Obs. 3', 'SWR1 Obs. 1', 'SWR1 Obs. 2', 'SWR1 Obs. 3'), loc='upper center',ncol=2,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
handles, labels = plta.get_legend_handles_labels()
leg = plta.legend((handles[0],handles[3],handles[1],handles[4],handles[2],handles[5]),
                  (labels[0],labels[3],labels[1],labels[4],labels[2],labels[5]),
                  loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 


outfigpng = '..\\Figures\\SWR01Sample03.png'
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
	show()
else:
	print 'no screen output requested...'

   