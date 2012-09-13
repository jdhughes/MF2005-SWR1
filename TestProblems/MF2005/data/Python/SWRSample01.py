import sys
import string
import numpy as np
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

#import arrayUtil as au
import MFBinaryClass as mfb 
reload(mfb)

def BedientHuber():
	d = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], np.float)
	s = np.array([0, 1, 2, 3, 4, 5, 6, 7.5, 10.5, 12, 13.5, 20, 22], np.float)
	q = np.array([0, 15, 32, 55, 90, 125, 158, 185, 210, 230, 250, 270, 290], np.float)
	si = np.zeros(len(q))
	inf = np.array([0, 0, 60, 120, 180, 240, 300, 360, 320, 280, 240, 200, 160, 120, 80, 40, 0, 0, 0, 0, 0, 0, 0], np.float)
	tmin = np.linspace(0, 220, 23)
	t = np.zeros(len(tmin))
	q1 = np.zeros(len(tmin))
	s1 = np.zeros(len(tmin))
	d1 = np.zeros(len(tmin))
	
	#-convert units
	for n in range(0,len(d)):
		#-ft to m
		d[n] = d[n] / 3.28081
		#-ac-ft to m3
		s[n] = s[n] * 43560. / pow(3.28081,3)
		#-ft3/s to m3/s
		q[n] = q[n] / pow(3.28081,3)

	for n in range(0,len(t)):
		#-min to s
		t[n] = tmin[n] * 60.
		#-ft3/s to m3/s
		inf[n] = inf[n] / pow(3.28081,3)
	
	s0 = 0.0
	q0 = 0.0
	for n in range(1,len(t)):
		dt = t[n] - t[n-1]
#-calculate si for current dt
		for k in range(0,len(d)):
			si[k] = ( s[k] / dt ) + q[k]
#-calculate backward in time outflow, storage, and depth
# at the end of the timestep
		rhs = inf[n] + s0 / dt
		q1[n] = xinterp(q,si,rhs)
		s1[n] = xinterp(s,q,q1[n])
		d1[n] = xinterp(d,s,s1[n])
		s0 = s1[n]

	return tmin,inf,q1,s1,d1
	
def xinterp(a,b,x):
	nlen = len(b)
	if x <= b[0]:
		t = a[0]
		return t
	if x > b[nlen-1]:
		t = a[nlen-1]
		return t
	for n in range(1,nlen):
		if b[n] == x:
			t = a[n]
			return t
		if b[n] > x:
			f = ( x - b[n-1] ) / ( b[n] - b[n-1] )
			t = a[n-1] + f * ( a[n] - a[n-1] )
			return t

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

fig = mpl.pyplot.figure(figsize=(8, 8), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.075,right=0.975,bottom=0.075,top=0.975)

SWR_file = '..\\SWRSample01\\Results\\SWRSample01_PoolFlow.flow'
SWRObj = mfb.SWR_Record(-1,SWR_file)
itime  = SWRObj.get_item_number('totim')
idt    = SWRObj.get_item_number('dt')
istage = SWRObj.get_item_number('stage')
ibc    = SWRObj.get_item_number('qbcflow')
ce1 = SWRObj.get_gage(1)
nt =  np.shape(ce1)[0]
swr_t = np.zeros(nt) 
swr_s = np.zeros(nt) 
swr_q = np.zeros(nt) 
swr_cq = np.zeros(nt) 
for n  in range(0,nt):
	swr_t[n] = ce1[n,itime] / 60.
	swr_s[n] = ce1[n,istage]
	swr_q[n] = -ce1[n,ibc]

for n  in range(1,nt):
	swr_cq[n] = swr_cq[n-1] - ce1[n,ibc] * ce1[n,idt]

SWR_file = '..\\SWRSample01\\Results\\SWRSample01_PoolFlow.01min.flow'
SWRObj = mfb.SWR_Record(-1,SWR_file)
ce101 = SWRObj.get_gage(1)
nt01 =  np.shape(ce101)[0]
swr_t01 = np.zeros(nt01) 
swr_s01 = np.zeros(nt01) 
swr_q01 = np.zeros(nt01) 
swr_cq01 = np.zeros(nt01) 
for n  in range(0,nt01):
	swr_t01[n] = ce101[n,itime] / 60.
	swr_s01[n] = ce101[n,istage]
	swr_q01[n] = -ce101[n,ibc]

for n  in range(1,nt01):
	swr_cq01[n] = swr_cq01[n-1] - ce101[n,ibc] * ce101[n,idt]

t,inf,q,s,d = BedientHuber()

plta = fig.add_subplot(221)
plta.plot(t, inf, 'r', linewidth=4)
plta.plot(swr_t, swr_q, 'k', linewidth=2)
plta.plot(t, q, 'bo', linewidth=0)
plta.text(0.05,0.925,'A',transform = plta.transAxes,size=12)
plta.set_xlabel('Time, minutes')	
plta.set_ylabel(r'Flow, m$^3$/s')
i = set_sizexaxis(plta,'%16.6g',ticksize)
i = set_sizeyaxis(plta,'%16.6g',ticksize)
#-legend for BH and 10 min
leg = plta.legend(('Specified inflow rate', 'SWR1 outflow - 10 minute SWR1 time step', 'Bedient and Huber (1988) outflow'), loc='upper right',labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltb = fig.add_subplot(223)
pltb.plot(swr_t, swr_s, 'k', linewidth=2)
pltb.plot(t, d+1.52, 'bo', linewidth=0)
pltb.text(0.05,0.925,'B',transform = pltb.transAxes,size=12)
pltb.set_xlabel('Time, minutes')	
pltb.set_ylabel('Stage, m')
i = set_sizexaxis(pltb,'%16.6g',ticksize)
i = set_sizeyaxis(pltb,'%16.6g',ticksize)

pltc = fig.add_subplot(222)
pltc.plot(swr_t01, swr_cq01 / 1000., 'g', linewidth=2)
pltc.plot(swr_t, swr_cq / 1000., 'ko', linewidth=0)
pltc.text(0.05,0.925,'C',transform = pltc.transAxes,size=12)
pltc.set_xlabel('Time, minutes')	
pltc.set_ylabel(r'Cumulative outflow x 1,000, m$^3$')
i = set_sizexaxis(pltc,'%16.6g',ticksize)
i = set_sizeyaxis(pltc,'%16.6g',ticksize)
#-legend for BH and 10 min
leg = pltc.legend(('1 minute SWR1 time step', '10 minute SWR1 time step'), loc='lower right',labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltd = fig.add_subplot(224)
pltd.plot(swr_t01, swr_s01, 'g', linewidth=2)
pltd.plot(swr_t, swr_s, 'ko', linewidth=0)
pltd.text(0.05,0.925,'D',transform = pltd.transAxes,size=12)
pltd.set_xlabel('Time, minutes')	
pltd.set_ylabel(r'Stage, m')
i = set_sizexaxis(pltd,'%16.6g',ticksize)
i = set_sizeyaxis(pltd,'%16.6g',ticksize)

#

plta.set_xlim(0,200)
pltb.set_xlim(0,200)
pltc.set_xlim(0,200)
pltd.set_xlim(0,200)


outfigpng = '..\\Figures\\SWR01Sample01.png'
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
	show()
else:
	print 'no screen output requested...'

   