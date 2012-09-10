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


#--HEC-RAS Results
infile = '..\\SWRSample04\\HEC-RAS\\SWRSample04_HECRAS.csv'
hec = np.genfromtxt(infile,delimiter=',')
datasize = hec.shape
rows = datasize[0] - 2
cols = datasize[1]
hec_time = zeros(rows, float)
hec_q01 = zeros(rows, float)
hec_q02 = zeros(rows, float)
hec_q03 = zeros(rows, float)
hec_q04 = zeros(rows, float)
hec_q05 = zeros(rows, float)
hec_s06 = zeros(rows, float)
hec_s07 = zeros(rows, float)
#--parse HEC-RAS data
for i in range(0,rows):
	ipos = i + 2
	hec_time[i] = float( hec[ipos,0] )
	hec_q01[i] = float( hec[ipos,1] )
	hec_q02[i] = float( hec[ipos,2] )
	hec_q03[i] = float( hec[ipos,3] )
	hec_q04[i] = float( hec[ipos,4] )
	hec_q05[i] = float( hec[ipos,5] )
	hec_s06[i] = float( hec[ipos,6] )
	hec_s07[i] = float( hec[ipos,7] )

#-SWR Results
SWR_file = '..\\SWRSample04\\Results\\SWRSample04_GroupFlow.flow'
SWRObj = mfb.SWR_Record(-1,SWR_file)
itime  = SWRObj.get_item_number('totim')
istage = SWRObj.get_item_number('stage')
ilat   = SWRObj.get_item_number('qlatflow')
irai   = SWRObj.get_item_number('rain')
iqe    = SWRObj.get_item_number('qeflow')
ibc    = SWRObj.get_item_number('qbcflow')
ce = SWRObj.get_gage(1)
nt =  np.shape(ce)[0]
swr_t = np.zeros(nt) 
swr_q01 = np.zeros(nt) 
swr_l01 = np.zeros(nt) 
for n  in range(0,nt):
	swr_t[n] = ce[n,itime] / ( 60. * 60. )
	swr_q01[n] = ce[n,ilat]
	swr_l01[n] = ce[n,ilat]

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce = SWRObj.get_gage(5)
nt =  np.shape(ce)[0]
swr_l05 = np.zeros(nt) 
for n  in range(0,nt):
	swr_l05[n] = ce[n,ilat]

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce = SWRObj.get_gage(14)
nt =  np.shape(ce)[0]
swr_q02 = np.zeros(nt) 
swr_s06 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s06[n] = ce[n,istage]
	swr_q02[n] = -ce[n,iqe]

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce = SWRObj.get_gage(18)
nt =  np.shape(ce)[0]
swr_q03 = np.zeros(nt) 
swr_s07 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s07[n] = ce[n,istage]
	swr_q03[n] = -ce[n,ibc]

SWR_file = '..\\SWRSample04\\Results\\SWRSample04_Velocity.vel'
SWRObj = mfb.SWR_Record(-2,SWR_file)
iflow = SWRObj.get_item_number('flow')
ce = SWRObj.get_gage(rec_num=4,iconn=5)
#print ce
nt =  np.shape(ce)[0]
swr_q04 = np.zeros(nt) 
for n  in range(0,nt):
	swr_q04[n] = -ce[n,iflow]

SWRObj = mfb.SWR_Record(-2,SWR_file)
ce = SWRObj.get_gage(rec_num=5,iconn=6)
#print ce
nt =  np.shape(ce)[0]
swr_q05 = np.zeros(nt) 
for n  in range(0,nt):
	swr_q05[n] = -ce[n,iflow]


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

fig = pylab.figure(figsize=(8, 8), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.075,right=0.975,bottom=0.075,top=0.975)

plta = fig.add_subplot(221)
#plta.semilogx(t, inf, 'r', linewidth=4)
plta.semilogx(hec_time, hec_q01, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=6, label='HEC-RAS observation location 1')
plta.semilogx(hec_time, hec_q02, 'o', linewidth=0, markeredgecolor='b', markerfacecolor='w', markersize=6, label='HEC-RAS observation location 2')
plta.semilogx(hec_time, hec_q03, 'o', linewidth=0, markeredgecolor='r', markerfacecolor='w', markersize=6, label='HEC-RAS observation location 3')
plta.semilogx(swr_t, swr_q01, 'k', linewidth=2, label='SWR1 observation location 1')
plta.semilogx(swr_t, swr_q02, 'b', linewidth=2, label='SWR1 observation location 2')
plta.semilogx(swr_t, swr_q03, 'r', linewidth=2, label='SWR1 observation location 3')
text(0.05,0.925,'A',transform = plta.transAxes,size=12)
plta.set_xlabel('Time, hours')	
plta.set_ylabel('Flow, m$^3$/s')
plta.set_xlim(1,100)
i = set_sizexaxis(plta,'%16.6g',ticksize)
i = set_sizeyaxis(plta,'%16.6g',ticksize)
#-legend
#leg = plta.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
handles, labels = plta.get_legend_handles_labels()
leg = plta.legend((handles[0],handles[3],handles[1],handles[4],handles[2],handles[5]),
                  (labels[0],labels[3],labels[1],labels[4],labels[2],labels[5]),
                  loc = 'upper left', bbox_to_anchor = (0.01, 0.925), 
                  ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltb = fig.add_subplot(223)
pltb.semilogx(hec_time, hec_q05, 'o', linewidth=0, markeredgecolor='b', markerfacecolor='w', markersize=6, label='HEC-RAS observation location 5')
pltb.semilogx(hec_time, hec_q04, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=6, label='HEC-RAS observation location 4')
pltb.semilogx(swr_t, swr_q05, 'b', linewidth=2, label='SWR1 observation location 5')
pltb.semilogx(swr_t, swr_q04, 'k', linewidth=2, label='SWR1 observation location 4')
text(0.05,0.925,'B',transform = pltb.transAxes,size=12)
pltb.set_xlabel('Time, hours')	
pltb.set_ylabel(r'Flow, m$^3$/s')
pltb.set_xlim(1,100)
i = set_sizexaxis(pltb,'%16.6g',ticksize)
i = set_sizeyaxis(pltb,'%16.6g',ticksize)
#-legend
handles, labels = pltb.get_legend_handles_labels()
#print handles
#print labels
leg = pltb.legend((handles[1],handles[3],handles[0],handles[2]), 
                  (labels[1],labels[3],labels[0],labels[2]), 
                  loc = 'upper left', bbox_to_anchor = (0.01, 0.925), 
                  ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltc = fig.add_subplot(222)
pltc.semilogx(hec_time, hec_s06, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=6)
pltc.semilogx(swr_t, swr_s06, 'k', linewidth=2)
text(0.05,0.925,'C',transform = pltc.transAxes,size=12)
pltc.set_xlabel('Time, hours')	
pltc.set_ylabel(r'Stage, m')
pltc.set_xlim(1,100)
i = set_sizexaxis(pltc,'%16.6g',ticksize)
i = set_sizeyaxis(pltc,'%16.6g',ticksize)
#-legend
leg = pltc.legend(('HEC-RAS observation location 6','SWR1 observation location 6'), 
                   loc = 'upper left', bbox_to_anchor = (0.01, 0.925), 
                   ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltd = fig.add_subplot(224)
pltd.semilogx(hec_time, hec_s07, 'o', linewidth=0, markeredgecolor='k', markerfacecolor='w', markersize=6)
pltd.semilogx(swr_t, swr_s07, 'k', linewidth=2)
#pltd.semilogx(swr_t, swr_s, 'ko', linewidth=0)
text(0.05,0.925,'D',transform = pltd.transAxes,size=12)
pltd.set_xlabel('Time, hours')	
pltd.set_ylabel(r'Stage, m')
pltd.set_xlim(1,100)
i = set_sizexaxis(pltd,'%16.6g',ticksize)
i = set_sizeyaxis(pltd,'%16.6g',ticksize)
#-legend
leg = pltd.legend(('HEC-RAS observation location 7','SWR1 observation location 7'), 
                  loc = 'upper left', bbox_to_anchor = (0.01, 0.925), 
                  ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

fig2 = pylab.figure(figsize=(4, 4), facecolor='w')
#fig2.subplots_adjust(wspace=0.1,hspace=0.1,left=0.075,right=0.975,bottom=0.075,top=0.975)

plt2 = fig2.add_subplot(111)
plt2.plot(swr_t, swr_l01, 'b', linewidth=2)
plt2.plot(swr_t, swr_l05, 'r', linewidth=2)
plt2.set_xlabel('Time, hours')	
plt2.set_ylabel('Inflow rate, m$^3$/s')
plt2.set_xlim(0,25)
i = set_sizexaxis(plt2,'%16.6g',ticksize)
i = set_sizeyaxis(plt2,'%16.6g',ticksize)
#-legend
leg = plt2.legend(('Specified inflow rate reach 1','Specified inflow rate reach 5'), loc='upper left',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 


outfigpng = '..\Figures\SWR01Sample04.png'
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

outfigpng = '..\Figures\SWR01Sample04_Inflow.png'
fig2.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
	show()
else:
	print 'no screen output requested...'

   