import sys
import string
import numpy as np
import matplotlib as mpl
from matplotlib.font_manager import FontProperties

import MFBinaryClass as mfb 
reload(mfb)

def averageq(a,ni):
	nt = len(a)
	b = np.zeros(nt)
	for n  in range(ni,nt-ni):
		t = 0.0
		for nn in range(-ni,ni):
			ipos = n - nn
			t = t + a[ipos]
		b[n] = t / ( float(ni+1) * 2.0 )
	return b

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
FigDir = '..\\Figures\\'
BaseFigName = 'SWR01Sample05'
FigExt = '.png'
narg = len(sys.argv)
iarg = 0
if narg > 1:
    while iarg < narg-1:
        iarg += 1
        basearg = sys.argv[iarg].lower()
        if basearg == '-s':
		    OutputToScreen == True
            #print 'command line arg: -s OutputToScreen = ', OutputToScreen
        elif basearg == '-o':
            try:
                iarg += 1
                BaseFigName = sys.argv[iarg]
                print 'command line arg: -o BaseFigName = ', BaseFigName
            except:
                print 'cannot parse command line arg: -o BaseFigName'

#--ASCII Rainfall Data
inperday2mmperhr = 25.4 / 24.0
cfil_rain = '..\\SWRSample05\\Rainfall\\DailyRainfall.csv'
rain0 = np.loadtxt(cfil_rain,usecols=[2,4],delimiter=',',skiprows=1)
datasize = rain0.shape
ndays = datasize[0]

rain_t   = np.zeros( (ndays), np.float )
rain_r01 = np.zeros( (ndays), np.float )
rain_r02 = np.zeros( (ndays), np.float )
for it in range(0,ndays):
    rain_t[it]   = float(it) + 1.0
    rain_r01[it] = rain0[it,0] * inperday2mmperhr
    rain_r02[it] = rain0[it,1] * inperday2mmperhr

cfil_bc = '..\\SWRSample05\\Rainfall\\WetlandStage.csv'
bc0 = np.loadtxt(cfil_bc,usecols=[1],delimiter=',',skiprows=1)
bc_s01 = np.zeros( (ndays), np.float )
for it in range(0,ndays):
    bc_s01[it] =  bc0[it]

#--ASCII cross-section data
infile = '..\\SWRSample05\\IrregularCrossSection_Reach17.dat'
cs = np.genfromtxt(infile,comments='#')
datasize = cs.shape
rows = datasize[0]
cols = datasize[1]
csr17x = np.zeros(rows, float)
csr17z = np.zeros(rows, float)
#--parse cross-section data
for i in range(0,rows):
	csr17x[i] = float( cs[i,0] )
	csr17z[i] = float( cs[i,1] )

infile = '..\\SWRSample05\\IrregularCrossSection_Reach18.dat'
cs = np.genfromtxt(infile,comments='#')
datasize = cs.shape
rows = datasize[0]
cols = datasize[1]
csr18x = np.zeros(rows, float)
csr18z = np.zeros(rows, float)
#--parse cross-section data
for i in range(0,rows):
	csr18x[i] = float( cs[i,0] )
	csr18z[i] = float( cs[i,1] )


#--SWR results
#--stages and flows
SWR_file = '..\\SWRSample05\\Results\\SWRSample05_GroupFlow.flow'
SWRObj = mfb.SWR_Record(-1,SWR_file)
itime  = SWRObj.get_item_number('totim')
istage = SWRObj.get_item_number('stage')
iuz    = SWRObj.get_item_number('quzflow')
ibf    = SWRObj.get_item_number('qbflow')
iqe    = SWRObj.get_item_number('qeflow')
ibc    = SWRObj.get_item_number('qbcflow')
ce1 = SWRObj.get_gage(1)
nt =  np.shape(ce1)[0]
swr_t = np.zeros(nt) 
swr_pct = np.zeros(nt) 
nt =  np.shape(ce1)[0]
swr_q05 = np.zeros(nt) 
for n  in range(0,nt):
	swr_t[n] = ce1[n,itime] #* 24.
	swr_pct[n] = 100.0 * float( nt - n ) / float( nt )
	swr_q05[n] = -ce1[n,iqe] / 86400.

#--add other inflow to reach group 3
SWRObj = mfb.SWR_Record(-1,SWR_file)
ce1 = SWRObj.get_gage(2)
for n  in range(0,nt):
	swr_q05[n] -= ce1[n,iqe] / 86400.


SWRObj = mfb.SWR_Record(-1,SWR_file)
ce1 = SWRObj.get_gage(3)
nt =  np.shape(ce1)[0]
swr_s01 = np.zeros(nt) 
swr_q06 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s01[n] = ce1[n,istage]
	swr_q06[n] = -ce1[n,iqe] / 86400.

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce1 = SWRObj.get_gage(4)
nt =  np.shape(ce1)[0]
swr_s02 = np.zeros(nt) 
swr_q07 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s02[n] = ce1[n,istage]
	swr_q07[n] = -ce1[n,iqe] / 86400.

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce1 = SWRObj.get_gage(6)
nt =  np.shape(ce1)[0]
swr_s03 = np.zeros(nt) 
swr_q08 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s03[n] = ce1[n,istage]
	swr_q08[n] = -ce1[n,ibc] / 86400.

SWRObj = mfb.SWR_Record(-1,SWR_file)
ce1 = SWRObj.get_gage(7)
nt =  np.shape(ce1)[0]
swr_s04 = np.zeros(nt) 
for n  in range(0,nt):
	swr_s04[n] = ce1[n,istage]

#--total runoff
swr_ro = np.zeros(nt) 
swr_bf = np.zeros(nt) 
for irg in range(3,8):
    SWRObj = mfb.SWR_Record(-1,SWR_file)
    rg1 = SWRObj.get_gage(irg)
    for n  in range(0,nt):
	    swr_ro[n] += rg1[n,iuz] / 86400.
	    swr_bf[n] += rg1[n,ibf] / 86400.
    
for n  in range(0,nt):
    swr_q05[n] = max( 1e-5, swr_q05[n] )
    swr_q06[n] = max( 1e-5, swr_q06[n] )
    swr_q07[n] = max( 1e-5, swr_q07[n] )
    swr_q08[n] = max( 1e-5, swr_q08[n] )
    swr_ro[n]  = max( 1e-5, swr_ro[n]  )

#--aquifer-reach exchanges
iobs01 = 14
iobs02 = 16
iobs03 = 18
iobs04 = 19
#iextbf = 15
SWR_file = '..\\SWRSample05\\Results\\SWRSample05_QAQ.bin'
SWRObj = mfb.SWR_Record(1,SWR_file)
itime  = SWRObj.get_item_number('totim')
iextbf = SWRObj.get_item_number('aq-rchflow')
ce1 = SWRObj.get_gage(iobs01)
#print 'qaq01\n', ce1
nt =  np.shape(ce1)[0]
swr_qaqt = np.zeros(nt) 
swr_qaq01 = np.zeros(nt) 
for n  in range(0,nt):
	swr_qaqt[n] = ce1[n,itime] #* 24.
	swr_qaq01[n] = ce1[n,iextbf] #/ 86400.


SWRObj = mfb.SWR_Record(1,SWR_file)
ce1 = SWRObj.get_gage(iobs02)
#print 'qaq02\n', ce1
nt =  np.shape(ce1)[0]
swr_qaq02 = np.zeros(nt) 
for n  in range(0,nt):
	swr_qaq02[n] = ce1[n,iextbf] #/ 86400.

SWRObj = mfb.SWR_Record(1,SWR_file)
ce1 = SWRObj.get_gage(iobs03)
#print 'qaq03 layer 1\n', ce1
nt =  np.shape(ce1)[0]
swr_qaq03 = np.zeros(nt) 
for n  in range(0,nt):
	swr_qaq03[n] = ce1[n,iextbf] #/ 86400.

SWRObj = mfb.SWR_Record(2,SWR_file)
ce1 = SWRObj.get_gage(iobs03)
#print 'qaq03 layer 2\n', ce1
nt =  np.shape(ce1)[0]
for n  in range(0,nt):
	swr_qaq03[n] = swr_qaq03[n] + ce1[n,iextbf] #/ 86400.

SWRObj = mfb.SWR_Record(1,SWR_file)
ce1 = SWRObj.get_gage(iobs04)
#print 'qaq04 layer 1\n', ce1
nt =  np.shape(ce1)[0]
swr_qaq04 = np.zeros(nt) 
for n  in range(0,nt):
	swr_qaq04[n] = ce1[n,iextbf] #/ 86400.

#SWRObj = mfb.SWR_Record(2,SWR_file)
#ce1 = SWRObj.get_gage(iobs04)
##print 'qaq04 layer 2\n', ce1
#nt =  np.shape(ce1)[0]
#for n  in range(0,nt):
#	swr_qaq04[n] = swr_qaq04[n] + ce1[n,iextbf] #/ 86400.

#--MODFLOW results
#--heads
ncol = 6
nrow = 6
nlay = 2

iobs1node = mfb.icrl_from_kij(1,3,4,nlay,nrow,ncol)
iobs2node = mfb.icrl_from_kij(1,4,5,nlay,nrow,ncol)
iobs3node = mfb.icrl_from_kij(1,5,6,nlay,nrow,ncol)
iobs4node = mfb.icrl_from_kij(1,6,5,nlay,nrow,ncol)
#print 'MODFLOW observation location', iobs1node, iobs2node, iobs3node, iobs4node

HEAD_file = '..\\SWRSample05\\Results\\SWRSample05.hds'

HEADObj = mfb.MODFLOW_Head(nlay,ncol,nrow,HEAD_file)
mf1 = HEADObj.get_gage(iobs1node)
nmft =  np.shape(mf1)[0]
mf_t = np.zeros(nmft) 
mf_h01 = np.zeros(nmft) 
t = np.zeros(nmft) 
for n  in range(0,nmft):
	mf_t[n] = mf1[n,0]
	t[n] = mf1[n,1]

dryv = t.min()  #-1.00000002e+30
mf_h01 = np.ma.masked_equal(t,dryv)

HEADObj = mfb.MODFLOW_Head(nlay,ncol,nrow,HEAD_file)
mf1 = HEADObj.get_gage(iobs2node)
mf_h02 = np.zeros(nmft) 
for n  in range(0,nmft):
	t[n] = mf1[n,1]
mf_h02 = np.ma.masked_equal(t,dryv)

HEADObj = mfb.MODFLOW_Head(nlay,ncol,nrow,HEAD_file)
mf1 = HEADObj.get_gage(iobs3node)
mf_h03 = np.zeros(nmft) 
for n  in range(0,nmft):
	t[n] = mf1[n,1]
mf_h03 = np.ma.masked_equal(t,dryv)

HEADObj = mfb.MODFLOW_Head(nlay,ncol,nrow,HEAD_file)
mf1 = HEADObj.get_gage(iobs4node)
mf_h04 = np.zeros(nmft) 
for n  in range(0,nmft):
	t[n] = mf1[n,1]
mf_h04 = np.ma.masked_equal(t,dryv)

#-Make figures
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

xr = [0, 365]

fig = mpl.pyplot.figure(figsize=(8, 8), facecolor='w')
fig.subplots_adjust(wspace=0.25,hspace=0.25,left=0.075,right=0.975,bottom=0.075,top=0.975)

plta = fig.add_subplot(221)
plta.plot(swr_t, swr_s01, 'b', linewidth=1)
plta.plot(swr_t, swr_s02, 'g', linewidth=1)
plta.plot(swr_t, swr_s03, 'm', linewidth=1)
plta.plot(swr_t, swr_s04, 'k', linewidth=1)
plta.text(0.05,0.925,'A',transform = plta.transAxes,size=12)
plta.set_xlabel('Time, days')	
plta.set_ylabel('Stage, m')
plta.set_xlim(0,350)
plta.set_ylim(0.4,1.6)
i = set_sizexaxis(plta,'%16.6g',ticksize)
i = set_sizeyaxis(plta,'%16.6g',ticksize)
#-legend
leg = plta.legend(('Observation location 1', 'Observation location 2', 'Observation location 3', 'Observation location 4'), loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltb = fig.add_subplot(223)
width = 0.5
pltb.plot(swr_t, swr_q08, 'r', linewidth=1, label='Observation location 8')
pltb.text(0.05,0.925,'B',transform = pltb.transAxes,size=12)
pltb.set_xlabel('Time, days')	
pltb.set_ylabel('Discharge, m$^3$/s')
pltb.set_xlim(0,350)
pltb.set_ylim(0,8.0)
i = set_sizexaxis(pltb,'%16.6g',ticksize)
i = set_sizeyaxis(pltb,'%16.6g',ticksize)
#-legend
leg = pltb.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

pltc = fig.add_subplot(222)
pltc.plot(mf_t, mf_h01, 'b', linewidth=1, label='Observation location 1')
pltc.plot(mf_t, mf_h02, 'g', linewidth=1, label='Observation location 2')
pltc.plot(mf_t, mf_h03, 'm', linewidth=1, label='Observation location 3')
pltc.plot(mf_t, mf_h04, 'k', linewidth=1, label='Observation location 4')
pltc.text(0.05,0.925,'C',transform = pltc.transAxes,size=12)
pltc.set_xlabel('Time, days')	
pltc.set_ylabel('Groundwater head, m')
pltc.set_xlim(0,350)
pltc.set_ylim(0.4,1.6)
i = set_sizexaxis(pltc,'%16.6g',ticksize)
i = set_sizeyaxis(pltc,'%16.6g',ticksize)
#-legend
leg = pltc.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 


pltd = fig.add_subplot(224)
pltd.plot(swr_qaqt, swr_qaq01, 'b', linewidth=1, label='Observation location 1')
pltd.plot(swr_qaqt, swr_qaq02, 'g', linewidth=1, label='Observation location 2')
pltd.plot(swr_qaqt, swr_qaq03, 'm', linewidth=1, label='Observation location 3')
pltd.plot(swr_qaqt, swr_qaq04, 'k', linewidth=1, label='Observation location 4')
pltd.plot([0,368], [0,0], '0.5', linewidth=0.5, label='_Zero')
pltd.text(0.05,0.925,'D',transform = pltd.transAxes,size=12)
pltd.set_xlabel('Time, days')	
pltd.set_ylabel('Aquifer-reach exchange, m$^3$/d')
pltd.set_xlim(0,350)
pltd.set_ylim(-750,1000)
i = set_sizexaxis(pltd,'%16.6g',ticksize)
i = set_sizeyaxis(pltd,'%16.6g',ticksize)
#-legend
leg = pltd.legend(loc='lower left',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1)
leg._drawFrame=False 

#--runoff and flow duration
fig2 = mpl.pyplot.figure(figsize=(8, 4), facecolor='w')
fig2.subplots_adjust(wspace=0.25,hspace=0.30,left=0.100,right=0.975,bottom=0.100,top=0.900)
ax = fig2.add_subplot(211)
ax.plot(swr_t, swr_q08, 'r', linewidth=1, label='Observation location 8')
ax.plot(swr_t, swr_q05, 'm', linewidth=1, label='Observation location 5')
ax.plot(swr_t, swr_ro, 'c', linewidth=1, label='Total UZF1 runoff')
ax.set_xlim(0,350)
ax.text(0.025,0.85,'A',transform = ax.transAxes,size=12)
ax.set_xlabel('Time, days')	
ax.set_ylabel('Discharge, m$^3$/s')
i = set_sizexaxis(ax,'%16.7g',ticksize)
i = set_sizeyaxis(ax,'%16.7g',ticksize)
#-legend
leg = ax.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=2.0,numpoints=1)
leg._drawFrame=False 

ax = fig2.add_subplot(212)
ax.semilogy(swr_pct, np.sort(swr_q05), 'm', linewidth=1, label='Observation location 5')
ax.semilogy(swr_pct, np.sort(swr_q06), 'k', linewidth=1, label='Observation location 6')
ax.semilogy(swr_pct, np.sort(swr_q07), 'g', linewidth=1, label='Observation location 7')
ax.semilogy(swr_pct, np.sort(swr_q08), 'r', linewidth=1, label='Observation location 8')
ax.semilogy(swr_pct, np.sort(swr_ro),  'c', linewidth=1, label='Total UZF1 runoff')
ax.set_ylim(1e-5,100.0)
ax.text(0.025,0.85,'B',transform = ax.transAxes,size=12)
ax.set_xlabel('Percentage of time value was equalled or exceeded')	
ax.set_ylabel('Discharge, m$^3$/s')
i = set_sizexaxis(ax,'%16.8g',ticksize)
i = set_sizeyaxis(ax,'%16.8g',ticksize)
#-legend
leg = ax.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=2.0,numpoints=1)
leg._drawFrame=False 

#--rainfall and boundary condition
fig3 = mpl.pyplot.figure(figsize=(8, 6), facecolor='w')
fig3.subplots_adjust(wspace=0.25,hspace=0.25,left=0.075,right=0.975,bottom=0.075,top=0.975)
ax = fig3.add_subplot(311)
width = 1.0
ax.bar(rain_t, rain_r01, width=width, align='center', color='b', linewidth=0, label='Rainfall zone 1')
ax.set_xlim(0,350)
ax.set_ylim(0.0,3.5)
ax.text(0.025,0.85,'A',transform = ax.transAxes,size=12)
ax.set_ylabel('Rainfall rate, mm/h')
i = set_sizexaxis(ax,'%16.6g',ticksize)
i = set_sizeyaxis(ax,'%16.6g',ticksize)
#-legend
leg = ax.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=2.0,numpoints=1)
leg._drawFrame=False 

ax = fig3.add_subplot(312)
ax.bar(rain_t, rain_r02, width=width, align='center', color='r', linewidth=0, label='Rainfall zone 2')
ax.set_xlim(0,350)
ax.set_ylim(0.0,3.5)
ax.text(0.025,0.85,'B',transform = ax.transAxes,size=12)
ax.set_ylabel('Rainfall rate, mm/h')
i = set_sizexaxis(ax,'%16.6g',ticksize)
i = set_sizeyaxis(ax,'%16.6g',ticksize)
#-legend
leg = ax.legend(loc='upper right',ncol=1,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=2.0,numpoints=1)
leg._drawFrame=False 

ax = fig3.add_subplot(313)
ax.plot(rain_t, bc_s01, 'k', linewidth=1, label='Boundary stage')
ax.set_xlim(0,350)
ax.set_ylim(1.5,1.7)
ax.text(0.025,0.85,'C',transform = ax.transAxes,size=12)
ax.set_xlabel('Time, days')	
ax.set_ylabel('Boundary stage, m')
i = set_sizexaxis(ax,'%16.6g',ticksize)
i = set_sizeyaxis(ax,'%16.6g',ticksize)

fig4 = mpl.pyplot.figure(figsize=(8, 4.00), facecolor='w')
fig4.subplots_adjust(wspace=0.25,hspace=0.25,left=0.075,right=0.975,bottom=0.2,top=0.8)
ax = fig4.add_subplot(1,2,1)
ax.plot(csr17x, csr17z, 'ko', linewidth=2, label='Station-elevation data')
ax.plot(csr17x, csr17z, 'k-', linewidth=2, label='Inferred cross section profile')
ax.set_xlabel('Horizontal distance, m')	
ax.set_ylabel('Elevation, m')
ct = 'Cross section 6\nReach 17'
mpl.pyplot.title(ct)
ax.set_xlim(-1,51)
ax.set_ylim(0,3.5)
i = set_sizexaxis(ax,'%16.6g',ticksize)
i = set_sizeyaxis(ax,'%16.6g',ticksize)

ax = fig4.add_subplot(1,2,2)
ax.plot(csr18x, csr18z, 'ko', linewidth=2, label='Station-elevation data')
ax.plot(csr18x, csr18z, 'k-', linewidth=2, label='Inferred cross section profile')
ax.set_xlabel('Horizontal distance, m')	
ax.set_ylabel('Elevation, m')
ct = 'Cross section 7\nReach 18'
mpl.pyplot.title(ct)
ax.set_xlim(-1,51)
ax.set_ylim(0,3.5)
i = set_sizexaxis(ax,'%16.6g',ticksize)
i = set_sizeyaxis(ax,'%16.6g',ticksize)

ax = mpl.pyplot.subplot(1,2,1)
##-legend
leg = ax.legend(loc='center',ncol=2,labelspacing=0.25,columnspacing=1,handletextpad=0.5,handlelength=0.5,numpoints=1, bbox_to_anchor=(1.0, -0.25)  )
leg._drawFrame=False 

#-write figure
outfigpng = FigDir + BaseFigName + FigExt
fig.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

#-write figure
outfigpng = FigDir + BaseFigName + '_RO_FlowDuration' + FigExt
fig2.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

#-write figure
outfigpng = FigDir + BaseFigName + '_Rainfall' + FigExt
fig3.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

#-write figure
outfigpng = '..\\Figures\\SWR01Sample05_CrossSections.png'
fig4.savefig(outfigpng,dpi=300)
print 'created...', outfigpng

if OutputToScreen:
	show()
else:
	print 'no screen output requested...'

   