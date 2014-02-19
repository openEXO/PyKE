#!/usr/bin/env python

import kepmsg, kepstat
import numpy, scipy, pylab
import scipy.stats
from numpy import *
from pylab import *
from matplotlib import *

# -----------------------------------------------------------
# clean up x-axis of plot

def cleanx(time,logfile,verbose):

    status = 0
    try:
        time0 = float(int(time[0] / 100) * 100.0)
        if time0 < 2.4e6: time0 += 2.4e6
        timeout = time - time0
        label = 'BJD $-$ %d' % time0
    except:
        txt = 'ERROR -- KEPPLOT.CLEANX: cannot calculate plot scaling in x dimension'
        status = kepmsg.err(logfile,txt,verbose)
        label = ''

    return timeout, label, status

# -----------------------------------------------------------
# clean up y-axis of plot

def cleany(signal,cadence,logfile,verbose):

    status = 0
    try:
        signal /= cadence
	nrm = len(str(int(signal.max())))-1
	signal = signal / 10**nrm
	label = '10$^%d$ e$^-$ s$^{-1}$' % nrm
    except:
        txt = 'ERROR -- KEPPLOT.CLEANY: cannot calculate plot scaling in y dimension'
        status = kepmsg.err(logfile,txt,verbose)
        label = ''

    return signal, label, status

# -----------------------------------------------------------
# plot limits

def limits(x,y,logfile,verbose):

    status = 0
    try:
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xr = xmax - xmin
        yr = ymax - ymin
        x = insert(x,[0],[x[0]]) 
        x = append(x,[x[-1]])
        y = insert(y,[0],[0.0]) 
        y = append(y,0.0)
    except:
        txt = 'ERROR -- KEPPLOT.LIMITS: cannot calculate plot limits'
        status = kepmsg.err(logfile,txt,verbose)

    return x, y, xmin,  xmax, xr, ymin, ymax, yr, status

# -----------------------------------------------------------
# plot look and feel

def define(labelsize,ticksize,logfile,verbose):

    status = 0
    try:
        rc('text', usetex=True)
#        rc('font',**{'family':'sans-serif','sans-serif':['sans-serif']})
        params = {'backend': 'png',
                  'axes.linewidth': 2.5,
                  'axes.labelsize': labelsize,
                  'axes.font': 'sans-serif',
                  'axes.fontweight' : 'bold',
                  'text.fontsize': 12,
                  'legend.fontsize': 12,
                  'xtick.labelsize': ticksize,
                  'ytick.labelsize': ticksize}
        rcParams.update(params)
    except:
        pass

    return status

# -----------------------------------------------------------
# intensity scale limits of 1d array

def intScale1D(image,imscale):

    seterr(all="ignore") 
    nstat = 2; work2 = []
    work1 = array(sort(image),dtype=float32)
    for i in range(len(work1)):
        if 'nan' not in str(work1[i]).lower():
            work2.append(work1[i])
    work2 = array(work2,dtype=float32)
    if int(float(len(work2)) / 10 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 10 + 0.5)
    zmin = median(work2[:nstat])
    zmax = median(work2[-nstat:])
    if imscale == 'logarithmic':
        if zmin < 0.0: zmin = 100.0
        image = log10(image)
        zmin = log10(zmin)
        zmax = log10(zmax)
    if (imscale == 'squareroot'):
        if zmin < 0.0: zmin = 100.0
        image = sqrt(image)
        zmin = sqrt(zmin)
        zmax = sqrt(zmax)

    return image, zmin, zmax

# -----------------------------------------------------------
# intensity scale limits of 2d array

def intScale2D(image,imscale):

    seterr(all="ignore") 
    nstat = 2
    work1 = numpy.array([],dtype='float32')
    (ysiz,xsiz) = numpy.shape(image)
    for i in range(ysiz):
        for j in range(xsiz):
            if numpy.isfinite(image[i,j]) and image[i,j] > 0.0:
                work1 = numpy.append(work1,image[i,j])
    work2 = array(sort(work1))
    if int(float(len(work2)) / 1000 + 0.5) > nstat:
        nstat = int(float(len(work2)) / 1000 + 0.5)
    zmin = median(work2[:nstat])
    zmax = median(work2[-nstat:])
    if imscale == 'logarithmic':
        image = log10(image)
        zmin = log10(zmin)
        zmax = log10(zmax)
    if (imscale == 'squareroot'):
        image = sqrt(image)
        zmin = sqrt(zmin)
        zmax = sqrt(zmax)

    return image, zmin, zmax


# ------------------------------------------
# plot mask borders in CCD coordinates

def borders(maskimg,xdim,ydim,pixcoord1,pixcoord2,bit,lcolor,lstyle,lwidth):

    for i in range(1,ydim):
        for j in range(1,xdim):
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

# corner cases 

    for j in range(ydim):
        try:
            if kepstat.bitInBitmap(maskimg[j,0],bit) and not kepstat.bitInBitmap(maskimg[j-1,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[j+1,0],bit) and kepstat.bitInBitmap(maskimg[j,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) + 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j,0],bit):
            x = array([pixcoord1[0,j],pixcoord1[0,j]]) - 0.5
            try:
                y = array([pixcoord2[0,j],pixcoord2[0,j+1]]) - 0.5
            except:
                y = array([pixcoord2[0,j-1],pixcoord2[0,j]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j,xdim-1],bit):
            x = array([pixcoord1[xdim-1,j],pixcoord1[xdim-1,j]]) + 0.5
            try:
                y = array([pixcoord2[xdim-1,j],pixcoord2[xdim-1,j+1]]) - 0.5
            except:
                y = array([pixcoord2[xdim-1,j-1],pixcoord2[xdim-1,j]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
    for i in range(xdim):
        try:
            if kepstat.bitInBitmap(maskimg[0,i],bit) and not kepstat.bitInBitmap(maskimg[0,i-1],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) - 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[0,i+1],bit) and kepstat.bitInBitmap(maskimg[0,i],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) + 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0,i],bit):
            try:
                x = array([pixcoord1[i,0],pixcoord1[i+1,0]]) - 0.5
            except:
                x = array([pixcoord1[i-1,0],pixcoord1[i,0]]) + 0.5
            y = array([pixcoord2[i,0],pixcoord2[i,0]]) - 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim-1,i],bit):
            try:
                x = array([pixcoord1[i,ydim-1],pixcoord1[i+1,ydim-1]]) - 0.5
            except:
                x = array([pixcoord1[i-1,ydim-1],pixcoord1[i,ydim-1]]) - 0.5
            y = array([pixcoord2[i,ydim-1],pixcoord2[i,ydim-1]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            
    if kepstat.bitInBitmap(maskimg[ydim-1,xdim-1],bit):
        x = array([pixcoord1[xdim-2,ydim-1],pixcoord1[xdim-1,ydim-1]]) + 0.5
        y = array([pixcoord2[xdim-1,ydim-1],pixcoord2[xdim-1,ydim-1]]) + 0.5
        pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0,xdim-1],bit):
        x = array([pixcoord1[xdim-1,0],pixcoord1[xdim-1,0]]) + 0.5
        y = array([pixcoord2[xdim-1,0],pixcoord2[xdim-1,1]]) - 0.5
        pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    return


# ------------------------------------------
# plot mask borders in CCD coordinates

def PrfBorders(maskimg,xdim,ydim,pixcoord1,pixcoord2,bit,lcolor,lstyle,lwidth):

    for i in range(1,ydim):
        for j in range(1,xdim):
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                pylab.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i-1,j],bit):
                x = array([pixcoord1[j-1,i],pixcoord1[j,i]]) + 0.5
                y = array([pixcoord2[j,i],pixcoord2[j,i]]) - 0.5
                pylab.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if kepstat.bitInBitmap(maskimg[i,j],bit) and not kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                pylab.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            if not kepstat.bitInBitmap(maskimg[i,j],bit) and kepstat.bitInBitmap(maskimg[i,j-1],bit):
                x = array([pixcoord1[j,i],pixcoord1[j,i]]) - 0.5
                y = array([pixcoord2[j,i-1],pixcoord2[j,i]]) + 0.5
                pylab.plot(x*50,y*50,color=lcolor,linestyle=lstyle,linewidth=lwidth)

# corner cases 

    for j in range(ydim):
        try:
            if kepstat.bitInBitmap(maskimg[j,0],bit) and not kepstat.bitInBitmap(maskimg[j-1,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[j+1,0],bit) and kepstat.bitInBitmap(maskimg[j,0],bit):
                x = array([pixcoord1[0,j],pixcoord1[1,j]]) - 0.5
                y = array([pixcoord2[0,j],pixcoord2[0,j]]) + 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        if kepstat.bitInBitmap(maskimg[j,0],bit):
            x = array([pixcoord1[0,j],pixcoord1[0,j]]) - 0.5
            try:
                y = array([pixcoord2[0,j],pixcoord2[0,j+1]]) - 0.5
            except:
                y = array([pixcoord2[0,j-1],pixcoord2[0,j]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[j,xdim-1],bit):
            x = array([pixcoord1[xdim-1,j],pixcoord1[xdim-1,j]]) + 0.5
            try:
                y = array([pixcoord2[xdim-1,j],pixcoord2[xdim-1,j+1]]) - 0.5
            except:
                y = array([pixcoord2[xdim-1,j-1],pixcoord2[xdim-1,j]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
    for i in range(xdim):
        try:
            if kepstat.bitInBitmap(maskimg[0,i],bit) and not kepstat.bitInBitmap(maskimg[0,i-1],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) - 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        try:
            if not kepstat.bitInBitmap(maskimg[0,i+1],bit) and kepstat.bitInBitmap(maskimg[0,i],bit):
                x = array([pixcoord1[i,0],pixcoord1[i,0]]) + 0.5
                y = array([pixcoord2[i,0],pixcoord2[i,1]]) - 0.5
                pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth) 
        except:
            pass
        if kepstat.bitInBitmap(maskimg[0,i],bit):
            try:
                x = array([pixcoord1[i,0],pixcoord1[i+1,0]]) - 0.5
            except:
                x = array([pixcoord1[i-1,0],pixcoord1[i,0]]) + 0.5
            y = array([pixcoord2[i,0],pixcoord2[i,0]]) - 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
        if kepstat.bitInBitmap(maskimg[ydim-1,i],bit):
            try:
                x = array([pixcoord1[i,ydim-1],pixcoord1[i+1,ydim-1]]) - 0.5
            except:
                x = array([pixcoord1[i-1,ydim-1],pixcoord1[i,ydim-1]]) - 0.5
            y = array([pixcoord2[i,ydim-1],pixcoord2[i,ydim-1]]) + 0.5
            pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)
            
    if kepstat.bitInBitmap(maskimg[ydim-1,xdim-1],bit):
        x = array([pixcoord1[xdim-2,ydim-1],pixcoord1[xdim-1,ydim-1]]) + 0.5
        y = array([pixcoord2[xdim-1,ydim-1],pixcoord2[xdim-1,ydim-1]]) + 0.5
        pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    if kepstat.bitInBitmap(maskimg[0,xdim-1],bit):
        x = array([pixcoord1[xdim-1,0],pixcoord1[xdim-1,0]]) + 0.5
        y = array([pixcoord2[xdim-1,0],pixcoord2[xdim-1,1]]) - 0.5
        pylab.plot(x,y,color=lcolor,linestyle=lstyle,linewidth=lwidth)

    return
