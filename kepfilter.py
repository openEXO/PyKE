
import numpy, sys, time, pyfits, pylab
from pyfits import *
from pylab import *
from matplotlib import *
from math import *
import kepio, kepmsg, kepkey, kepfunc, kepstat

def kepfilter(infile,outfile,datacol,function,cutoff,passband,plot,plotlab,
              clobber,verbose,logfile,status,cmdLine=False): 

## startup parameters

    status = 0
    numpy.seterr(all="ignore") 
    labelsize = 24
    ticksize = 16
    xsize = 16
    ysize = 6
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2

## log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPFILTER -- '
    call += 'infile='+infile+' '
    call += 'outfile='+outfile+' '
    call += 'datacol='+str(datacol)+' '
    call += 'function='+str(function)+' '
    call += 'cutoff='+str(cutoff)+' '
    call += 'passband='+str(passband)+' '
    plotit = 'n'
    if (plot): plotit = 'y'
    call += 'plot='+plotit+ ' '
    call += 'plotlab='+str(plotlab)+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

## start time

    kepmsg.clock('KEPFILTER started at',logfile,verbose)

## test log file

    logfile = kepmsg.test(logfile)

## clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
	    message = 'ERROR -- KEPFILTER: ' + outfile + ' exists. Use clobber=yes'
	    status = kepmsg.err(logfile,message,verbose)

## open input file

    if status == 0:
        instr, status = kepio.openfits(infile,'readonly',logfile,verbose)
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)
    if status == 0:
        try:
            work = instr[0].header['FILEVER']
            cadenom = 1.0
        except:
            cadenom = cadence

## fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

## read table structure

    if status == 0:
	table, status = kepio.readfitstab(infile,instr[1],logfile,verbose)

# read time and flux columns

    if status == 0:
        barytime, status = kepio.readtimecol(infile,table,logfile,verbose)
        flux, status = kepio.readsapcol(infile,table,logfile,verbose)

# filter input data table

    if status == 0:
        try:
            nanclean = instr[1].header['NANCLEAN']
        except:
            naxis2 = 0
            for i in range(len(table.field(0))):
                if (numpy.isfinite(barytime[i]) and numpy.isfinite(flux[i]) and flux[i] != 0.0):
                    table[naxis2] = table[i]
                    naxis2 += 1
            instr[1].data = table[:naxis2]
            comment = 'NaN cadences removed from data'
            status = kepkey.new('NANCLEAN',True,comment,instr[1],outfile,logfile,verbose)

## read table columns

    if status == 0:
        intime, status = kepio.readtimecol(infile,instr[1].data,logfile,verbose)
    if status == 0:
	indata, status = kepio.readfitscol(infile,instr[1].data,datacol,logfile,verbose)
    if status == 0:
        intime = intime + bjdref
        indata = indata / cadenom

## define data sampling

    if status == 0:
        tr = 1.0 / (cadence / 86400)
        timescale = 1.0 / (cutoff / tr)

## define convolution function

    if status == 0:
        if function == 'boxcar':
            filtfunc = numpy.ones(numpy.ceil(timescale))
        elif function == 'gauss':
            timescale /= 2
            dx = numpy.ceil(timescale * 10 + 1)
            filtfunc = kepfunc.gauss()
            filtfunc = filtfunc([1.0,dx/2-1.0,timescale],linspace(0,dx-1,dx))
        elif function == 'sinc':
            dx = numpy.ceil(timescale * 12 + 1)
            fx = linspace(0,dx-1,dx)
            fx = fx - dx / 2 + 0.5
            fx /= timescale
            filtfunc = numpy.sinc(fx)
        filtfunc /= numpy.sum(filtfunc)

## pad time series at both ends with noise model

    if status == 0:
        ave, sigma  = kepstat.stdev(indata[:len(filtfunc)])
        padded = append(kepstat.randarray(np.ones(len(filtfunc)) * ave,
                                          np.ones(len(filtfunc)) * sigma), indata)
        ave, sigma  = kepstat.stdev(indata[-len(filtfunc):])
        padded = append(padded, kepstat.randarray(np.ones(len(filtfunc)) * ave,
                                                  np.ones(len(filtfunc)) * sigma))

## convolve data

    if status == 0:
        convolved = convolve(padded,filtfunc,'same')

## remove padding from the output array

    if status == 0:
        if function == 'boxcar':
            outdata = convolved[len(filtfunc):-len(filtfunc)]
        else:
            outdata = convolved[len(filtfunc):-len(filtfunc)]
            

## subtract low frequencies

    if status == 0 and passband == 'high':
        outmedian = median(outdata)
        outdata = indata - outdata + outmedian

## comment keyword in output file

    if status == 0:
        status = kepkey.history(call,instr[0],outfile,logfile,verbose)

## clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
        if intime0 < 2.4e6: intime0 += 2.4e6
	ptime = intime - intime0
	xlab = 'BJD $-$ %d' % intime0

## clean up y-axis units

    if status == 0:
        pout = indata * 1.0
        pout2 = outdata * 1.0
	nrm = len(str(int(numpy.nanmax(pout))))-1
	pout = pout / 10**nrm
	pout2 = pout2 / 10**nrm
	ylab = '10$^%d$ %s' % (nrm, plotlab)

## data limits

	xmin = ptime.min()
	xmax = ptime.max()
	ymin = numpy.nanmin(pout)
	ymax = numpy.nanmax(pout)
	xr = xmax - xmin
	yr = ymax - ymin
        ptime = insert(ptime,[0],[ptime[0]]) 
        ptime = append(ptime,[ptime[-1]])
        pout = insert(pout,[0],[0.0]) 
        pout = append(pout,0.0)
        pout2 = insert(pout2,[0],[0.0]) 
        pout2 = append(pout2,0.0)

## plot light curve

    if status == 0 and plot:
        try:
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
            print 'ERROR -- KEPFILTER: install latex for scientific plotting'
            status = 1
    if status == 0 and plot:
        pylab.figure(figsize=[xsize,ysize])
        pylab.clf()

## plot filtered data

        ax = pylab.axes([0.06,0.1,0.93,0.87])
        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        labels = ax.get_yticklabels()
        setp(labels, 'rotation', 90, fontsize=12)
        pylab.plot(ptime,pout,color='#ff9900',linestyle='-',linewidth=lwidth)
        fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
        if passband == 'low':
            pylab.plot(ptime[1:-1],pout2[1:-1],color=lcolor,linestyle='-',linewidth=lwidth)
        else:
            pylab.plot(ptime,pout2,color=lcolor,linestyle='-',linewidth=lwidth)
            fill(ptime,pout2,color=lcolor,linewidth=0.0,alpha=falpha)
	xlabel(xlab, {'color' : 'k'})
	ylabel(ylab, {'color' : 'k'})
	xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0: 
            ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            ylim(1.0e-10,ymax+yr*0.01)
        pylab.grid()
        
# render plot

        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
## write output file

    if status == 0:
        for i in range(len(outdata)):
            instr[1].data.field(datacol)[i] = outdata[i]
        instr.writeto(outfile)
    
## close input file

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

## end time

    if (status == 0):
	    message = 'KEPFILTER completed at'
    else:
	    message = '\nKEPFILTER aborted at'
    kepmsg.clock(message,logfile,verbose)

## main
if '--shell' in sys.argv:
    import argparse
    
    parser = argparse.ArgumentParser(description='Low bandpass or high bandpass signal filtering')
    parser.add_argument('--shell', action='store_true', help='Are we running from the shell?')

    parser.add_argument('infile', help='Name of input file', type=str)
    parser.add_argument('outfile', help='Name of FITS file to output', type=str)

    parser.add_argument('--datacol', default='SAP_FLUX', help='Name of data column', type=str)
    parser.add_argument('--function', default='boxcar', help='The bandpass convolution function', type=str, 
        choices=['boxcar','gauss','sinc'])

    parser.add_argument('--cutoff', default=1.0, help='Characteristic frequency cutoff of filter [1/days]', type=float)
    parser.add_argument('--passband', help='low- or high-bandpass filter', type=str, choices=['low','high'])
    
    parser.add_argument('--plot', action='store_true', help='Plot result?')



    parser.add_argument('--clobber', action='store_true', help='Overwrite output file?')
    parser.add_argument('--verbose', action='store_true', help='Write to a log file?')
    parser.add_argument('--logfile', '-l', help='Name of ascii log file', default='kepcotrend.log', dest='logfile', type=str)
    parser.add_argument('--status', '-e', help='Exit status (0=good)', default=0, dest='status', type=int)


    args = parser.parse_args()
    
    cmdLine=True

    kepfilter(args.infile,args.outfile,args.datacol,args.function,args.cutoff,args.passband,args.plot,'e$^-$ s$^{-1}$',
              args.clobber,args.verbose,args.logfile,args.status,cmdLine)
    

else:
    from pyraf import iraf
    parfile = iraf.osfn("kepler$kepfilter.par")
    t = iraf.IrafTaskFactory(taskname="kepfilter", value=parfile, function=kepfilter)
