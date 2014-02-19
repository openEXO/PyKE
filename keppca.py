from pyraf import iraf
import pylab, numpy, pyfits, scipy
from scipy import stats
from scipy.stats import stats
from pylab import *
from matplotlib import *
from numpy import *
from pyfits import *
import kepio, kepmsg, kepstat, kepkey
import sys, time

# -----------------------------------------------------------
# core code

def keppca(infile,maskfile,outfile,components,clobber,verbose,logfile,status): 

# startup parameters

    cmdLine=False
    status = 0
    labelsize = 32
    ticksize = 18
    xsize = 16
    ysize = 10
    lcolor = '#0000ff'
    lwidth = 1.0
    fcolor = '#ffff00'
    falpha = 0.2
    seterr(all="ignore") 

# log the call 

    hashline = '----------------------------------------------------------------------------'
    kepmsg.log(logfile,hashline,verbose)
    call = 'KEPPCA -- '
    call += 'infile='+infile+' '
    call += 'maskfile='+maskfile+' '
    call += 'outfile='+outfile+' '
    call += 'components='+components+' '
    overwrite = 'n'
    if (clobber): overwrite = 'y'
    call += 'clobber='+overwrite+ ' '
    chatter = 'n'
    if (verbose): chatter = 'y'
    call += 'verbose='+chatter+' '
    call += 'logfile='+logfile
    kepmsg.log(logfile,call+'\n',verbose)

# start time

    kepmsg.clock('KEPPCA started at',logfile,verbose)

# test log file

    logfile = kepmsg.test(logfile)

# clobber output file

    if clobber: status = kepio.clobber(outfile,logfile,verbose)
    if kepio.fileexists(outfile): 
        message = 'ERROR -- KEPPCA: ' + outfile + ' exists. Use --clobber'
        status = kepmsg.err(logfile,message,verbose)

# open input file

    status = 0
    instr = pyfits.open(infile,mode='readonly',memmap=True)
    if status == 0:
        tstart, tstop, bjdref, cadence, status = kepio.timekeys(instr,infile,logfile,verbose,status)

# fudge non-compliant FITS keywords with no values

    if status == 0:
        instr = kepkey.emptykeys(instr,file,logfile,verbose)

# input file data

    if status == 0:
        cards0 = instr[0].header.ascardlist()
        cards1 = instr[1].header.ascardlist()
        cards2 = instr[2].header.ascardlist()
        table = instr[1].data[:]
        maskmap = copy(instr[2].data)

# open TPF FITS file

    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, barytime, status = \
            kepio.readTPF(infile,'TIME',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, tcorr, status = \
            kepio.readTPF(infile,'TIMECORR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, cadno, status = \
            kepio.readTPF(infile,'CADENCENO',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, fluxpixels, status = \
            kepio.readTPF(infile,'FLUX',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, errpixels, status = \
            kepio.readTPF(infile,'FLUX_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, flux_bkg, status = \
            kepio.readTPF(infile,'FLUX_BKG',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, flux_bkg_err, status = \
            kepio.readTPF(infile,'FLUX_BKG_ERR',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, qual, status = \
            kepio.readTPF(infile,'QUALITY',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, pcorr1, status = \
            kepio.readTPF(infile,'POS_CORR1',logfile,verbose)
    if status == 0:
        kepid, channel, skygroup, module, output, quarter, season, \
            ra, dec, column, row, kepmag, xdim, ydim, pcorr2, status = \
            kepio.readTPF(infile,'POS_CORR2',logfile,verbose)

# read mask definition file

    if status == 0 and 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        maskx = array([],'int')
        masky = array([],'int')
        lines, status = kepio.openascii(maskfile,'r',logfile,verbose)
        for line in lines:
            line = line.strip().split('|')
            if len(line) == 6:
                y0 = int(line[3])
                x0 = int(line[4])
                line = line[5].split(';')
                for items in line:
                    try:
                        masky = numpy.append(masky,y0 + int(items.split(',')[0]))
                        maskx = numpy.append(maskx,x0 + int(items.split(',')[1]))
                    except:
                        continue
        status = kepio.closeascii(lines,logfile,verbose)
        if len(maskx) == 0 or len(masky) == 0:
            message = 'ERROR -- KEPPCA: ' + maskfile + ' contains no pixels.'
            status = kepmsg.err(logfile,message,verbose)

# subimage physical WCS data

    if status == 0:
        crpix1p = cards2['CRPIX1P'].value
        crpix2p = cards2['CRPIX2P'].value
        crval1p = cards2['CRVAL1P'].value
        crval2p = cards2['CRVAL2P'].value
        cdelt1p = cards2['CDELT1P'].value
        cdelt2p = cards2['CDELT2P'].value

# define new subimage bitmap...

    if status == 0 and 'aper' not in maskfile.lower() and maskfile.lower() != 'all':
        aperx = numpy.array([],'int')
        apery = numpy.array([],'int')
        aperb = numpy.array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperx = numpy.append(aperx,crval1p + (j + 1 - crpix1p) * cdelt1p)
                apery = numpy.append(apery,crval2p + (i + 1 - crpix2p) * cdelt2p)
                if maskmap[i,j] == 0:
                    aperb = numpy.append(aperb,0)
                else:
                    aperb = numpy.append(aperb,1)
                    maskmap[i,j] = 1
                    for k in range(len(maskx)):
                        if aperx[-1] == maskx[k] and apery[-1] == masky[k]:
                            aperb[-1] = 3
                            maskmap[i,j] = 3

# ...or use old subimage bitmap

    if status == 0 and 'aper' in maskfile.lower():
        aperb = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                aperb = numpy.append(aperb,maskmap[i,j])

# ...or use all pixels

    if status == 0 and maskfile.lower() == 'all':
        aperb = array([],'int')
        for i in range(maskmap.shape[0]):
            for j in range(maskmap.shape[1]):
                if maskmap[i,j] == 0:
                    aperb = numpy.append(aperb,0)
                else:
                    aperb = numpy.append(aperb,3)
                    maskmap[i,j] = 3

# legal mask defined?

    if status == 0:
        if len(aperb) == 0:
            message = 'ERROR -- KEPPCA: no legal pixels within the subimage are defined.'
            status = kepmsg.err(logfile,message,verbose)
        
# identify principal components to be combined

    if status == 0:
        pcaout = []
        txt = components.strip().split(',')
        for work1 in txt:
            try:
                pcaout.append(int(work1.strip()))
            except:
                work2 = work1.strip().split('-')
                try:
                    for work3 in range(int(work2[0]),int(work2[1]) + 1):
                        pcaout.append(work3)
                except:
                    message = 'ERROR -- KEPPCA: cannot understand principal component list requested'
                    status = kepmsg.err(logfile,message,verbose)
    if status == 0:
        pcaout = set(sort(pcaout))

# flux pixel array size

    if status == 0:
        ntim = 0
        time = numpy.array([],dtype='float64')
        timecorr = numpy.array([],dtype='float32')
        cadenceno = numpy.array([],dtype='int')
        pixseries = numpy.array([],dtype='float32')
        errseries = numpy.array([],dtype='float32')
        bkgseries = numpy.array([],dtype='float32')
        berseries = numpy.array([],dtype='float32')
        quality = numpy.array([],dtype='float32')
        pos_corr1 = numpy.array([],dtype='float32')
        pos_corr2 = numpy.array([],dtype='float32')
        nrows = numpy.size(fluxpixels,0)
        npix = numpy.size(fluxpixels,1)

# remove NaN timestamps

        for i in range(nrows):
            if qual[i] == 0 and \
                    numpy.isfinite(barytime[i]) and \
                    numpy.isfinite(fluxpixels[i,ydim*xdim/2]) and \
                    numpy.isfinite(fluxpixels[i,1+ydim*xdim/2]):
                ntim += 1
                time = numpy.append(time,barytime[i])
                timecorr = numpy.append(timecorr,tcorr[i])
                cadenceno = numpy.append(cadenceno,cadno[i])
                pixseries = numpy.append(pixseries,fluxpixels[i])
                errseries = numpy.append(errseries,errpixels[i])
                bkgseries = numpy.append(bkgseries,flux_bkg[i])
                berseries = numpy.append(berseries,flux_bkg_err[i])
                quality = numpy.append(quality,qual[i])
                pos_corr1 = numpy.append(pos_corr1,pcorr1[i])
                pos_corr2 = numpy.append(pos_corr2,pcorr2[i])
        pixseries = numpy.reshape(pixseries,(-1,npix))
        errseries = numpy.reshape(errseries,(-1,npix))
        bkgseries = numpy.reshape(bkgseries,(-1,npix))
        berseries = numpy.reshape(berseries,(-1,npix))

# dummy columns for output file

    if status == 0:
        pdc_flux = numpy.empty(len(time)); pdc_flux[:] = numpy.nan
        pdc_flux_err = numpy.empty(len(time)); pdc_flux_err[:] = numpy.nan
        psf_centr1 = numpy.empty(len(time)); psf_centr1[:] = numpy.nan
        psf_centr1_err = numpy.empty(len(time)); psf_centr1_err[:] = numpy.nan
        psf_centr2 = numpy.empty(len(time)); psf_centr2[:] = numpy.nan
        psf_centr2_err = numpy.empty(len(time)); psf_centr2_err[:] = numpy.nan
        mom_centr1 = numpy.empty(len(time)); mom_centr1[:] = numpy.nan
        mom_centr1_err = numpy.empty(len(time)); mom_centr1_err[:] = numpy.nan
        mom_centr2 = numpy.empty(len(time)); mom_centr2[:] = numpy.nan
        mom_centr2_err = numpy.empty(len(time)); mom_centr2_err[:] = numpy.nan

# subtract mean over time from each pixel in the mask

    if status == 0:
        nmask = 0
        for i in range(npix):
            if aperb[i] == 3:
                nmask += 1
        work1 = numpy.zeros((len(pixseries),nmask))
        nmask = -1
        for i in range(npix):
            if aperb[i] == 3:
                nmask += 1
                maskedFlux = numpy.ma.masked_invalid(pixseries[:,i])
                pixMean = numpy.mean(maskedFlux)
                if numpy.isfinite(pixMean):
                    work1[:,nmask] = maskedFlux - pixMean
                else:
                    work1[:,nmask] = numpy.zeros((ntim))

# calculate covariance matrix

    if status == 0:
        work2 = work1.T
        covariance = numpy.cov(work2)

# determine eigenfunctions and eigenvectors of the covariance matrix
        
    if status == 0:
        [latent,coeff] = numpy.linalg.eig(covariance)

# projection of the data in the new space

    if status == 0:
        score = numpy.dot(coeff.T,work2).T

# construct new table data

    if status == 0:
        sap_flux = numpy.array([],'float32')
        sap_flux_err = numpy.array([],'float32')
        sap_bkg = numpy.array([],'float32')
        sap_bkg_err = numpy.array([],'float32')
        for i in range(len(time)):
            work1 = numpy.array([],'float64')
            work2 = numpy.array([],'float64')
            work3 = numpy.array([],'float64')
            work4 = numpy.array([],'float64')
            work5 = numpy.array([],'float64')
            for j in range(len(aperb)):
                if (aperb[j] == 3):
                    work1 = numpy.append(work1,pixseries[i,j])
                    work2 = numpy.append(work2,errseries[i,j])
                    work3 = numpy.append(work3,bkgseries[i,j])
                    work4 = numpy.append(work4,berseries[i,j])
            sap_flux = numpy.append(sap_flux,kepstat.sum(work1))
            sap_flux_err = numpy.append(sap_flux_err,kepstat.sumerr(work2))
            sap_bkg = numpy.append(sap_bkg,kepstat.sum(work3))
            sap_bkg_err = numpy.append(sap_bkg_err,kepstat.sumerr(work4))
        sap_mean = scipy.stats.stats.nanmean(sap_flux)

# coadd principal components

    if status == 0:
        pca_flux = numpy.zeros((len(sap_flux)))
        for i in range(nmask):
            if (i + 1) in pcaout:
                pca_flux = pca_flux + score[:,i]
        pca_flux += sap_mean

# construct output primary extension

    if status == 0:
        hdu0 = pyfits.PrimaryHDU()
        for i in range(len(cards0)):
            if cards0[i].key not in hdu0.header.ascardlist().keys():
                hdu0.header.update(cards0[i].key, cards0[i].value, cards0[i].comment)
            else:
                hdu0.header.ascardlist()[cards0[i].key].comment = cards0[i].comment
        status = kepkey.history(call,hdu0,outfile,logfile,verbose)
        outstr = HDUList(hdu0)

# construct output light curve extension

    if status == 0:
        col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
        col2 = Column(name='TIMECORR',format='E',unit='d',array=timecorr)
        col3 = Column(name='CADENCENO',format='J',array=cadenceno)
        col4 = Column(name='SAP_FLUX',format='E',array=sap_flux)
        col5 = Column(name='SAP_FLUX_ERR',format='E',array=sap_flux_err)
        col6 = Column(name='SAP_BKG',format='E',array=sap_bkg)
        col7 = Column(name='SAP_BKG_ERR',format='E',array=sap_bkg_err)
        col8 = Column(name='PDCSAP_FLUX',format='E',array=pdc_flux)
        col9 = Column(name='PDCSAP_FLUX_ERR',format='E',array=pdc_flux_err)
        col10 = Column(name='SAP_QUALITY',format='J',array=quality)
        col11 = Column(name='PSF_CENTR1',format='E',unit='pixel',array=psf_centr1)
        col12 = Column(name='PSF_CENTR1_ERR',format='E',unit='pixel',array=psf_centr1_err)
        col13 = Column(name='PSF_CENTR2',format='E',unit='pixel',array=psf_centr2)
        col14 = Column(name='PSF_CENTR2_ERR',format='E',unit='pixel',array=psf_centr2_err)
        col15 = Column(name='MOM_CENTR1',format='E',unit='pixel',array=mom_centr1)
        col16 = Column(name='MOM_CENTR1_ERR',format='E',unit='pixel',array=mom_centr1_err)
        col17 = Column(name='MOM_CENTR2',format='E',unit='pixel',array=mom_centr2)
        col18 = Column(name='MOM_CENTR2_ERR',format='E',unit='pixel',array=mom_centr2_err)
        col19 = Column(name='POS_CORR1',format='E',unit='pixel',array=pos_corr1)
        col20 = Column(name='POS_CORR2',format='E',unit='pixel',array=pos_corr2)
        cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11, \
                            col12,col13,col14,col15,col16,col17,col18,col19,col20])
        hdu1 = new_table(cols)
        hdu1.header.update('TTYPE1','TIME','column title: data time stamps')
        hdu1.header.update('TFORM1','D','data type: float64')
        hdu1.header.update('TUNIT1','BJD - 2454833','column units: barycenter corrected JD')
        hdu1.header.update('TDISP1','D12.7','column display format')
        hdu1.header.update('TTYPE2','TIMECORR','column title: barycentric-timeslice correction')
        hdu1.header.update('TFORM2','E','data type: float32')
        hdu1.header.update('TUNIT2','d','column units: days')
        hdu1.header.update('TTYPE3','CADENCENO','column title: unique cadence number')
        hdu1.header.update('TFORM3','J','column format: signed integer32')
        hdu1.header.update('TTYPE4','SAP_FLUX','column title: aperture photometry flux')
        hdu1.header.update('TFORM4','E','column format: float32')
        hdu1.header.update('TUNIT4','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE5','SAP_FLUX_ERR','column title: aperture phot. flux error')
        hdu1.header.update('TFORM5','E','column format: float32')
        hdu1.header.update('TUNIT5','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE6','SAP_BKG','column title: aperture phot. background flux')
        hdu1.header.update('TFORM6','E','column format: float32')
        hdu1.header.update('TUNIT6','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE7','SAP_BKG_ERR','column title: ap. phot. background flux error')
        hdu1.header.update('TFORM7','E','column format: float32')
        hdu1.header.update('TUNIT7','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE8','PDCSAP_FLUX','column title: PDC photometry flux')
        hdu1.header.update('TFORM8','E','column format: float32')
        hdu1.header.update('TUNIT8','e-/s','column units: electrons per second')
        hdu1.header.update('TTYPE9','PDCSAP_FLUX_ERR','column title: PDC flux error')
        hdu1.header.update('TFORM9','E','column format: float32')
        hdu1.header.update('TUNIT9','e-/s','column units: electrons per second (1-sigma)')
        hdu1.header.update('TTYPE10','SAP_QUALITY','column title: aperture photometry quality flag')
        hdu1.header.update('TFORM10','J','column format: signed integer32')
        hdu1.header.update('TTYPE11','PSF_CENTR1','column title: PSF fitted column centroid')
        hdu1.header.update('TFORM11','E','column format: float32')
        hdu1.header.update('TUNIT11','pixel','column units: pixel')
        hdu1.header.update('TTYPE12','PSF_CENTR1_ERR','column title: PSF fitted column error')
        hdu1.header.update('TFORM12','E','column format: float32')
        hdu1.header.update('TUNIT12','pixel','column units: pixel')
        hdu1.header.update('TTYPE13','PSF_CENTR2','column title: PSF fitted row centroid')
        hdu1.header.update('TFORM13','E','column format: float32')
        hdu1.header.update('TUNIT13','pixel','column units: pixel')
        hdu1.header.update('TTYPE14','PSF_CENTR2_ERR','column title: PSF fitted row error')
        hdu1.header.update('TFORM14','E','column format: float32')
        hdu1.header.update('TUNIT14','pixel','column units: pixel')
        hdu1.header.update('TTYPE15','MOM_CENTR1','column title: moment-derived column centroid')
        hdu1.header.update('TFORM15','E','column format: float32')
        hdu1.header.update('TUNIT15','pixel','column units: pixel')
        hdu1.header.update('TTYPE16','MOM_CENTR1_ERR','column title: moment-derived column error')
        hdu1.header.update('TFORM16','E','column format: float32')
        hdu1.header.update('TUNIT16','pixel','column units: pixel')
        hdu1.header.update('TTYPE17','MOM_CENTR2','column title: moment-derived row centroid')
        hdu1.header.update('TFORM17','E','column format: float32')
        hdu1.header.update('TUNIT17','pixel','column units: pixel')
        hdu1.header.update('TTYPE18','MOM_CENTR2_ERR','column title: moment-derived row error')
        hdu1.header.update('TFORM18','E','column format: float32')
        hdu1.header.update('TUNIT18','pixel','column units: pixel')
        hdu1.header.update('TTYPE19','POS_CORR1','column title: col correction for vel. abbern')
        hdu1.header.update('TFORM19','E','column format: float32')
        hdu1.header.update('TUNIT19','pixel','column units: pixel')
        hdu1.header.update('TTYPE20','POS_CORR2','column title: row correction for vel. abbern')
        hdu1.header.update('TFORM20','E','column format: float32')
        hdu1.header.update('TUNIT20','pixel','column units: pixel')
        hdu1.header.update('EXTNAME','LIGHTCURVE','name of extension')
        for i in range(len(cards1)):
            if (cards1[i].key not in hdu1.header.ascardlist().keys() and
                cards1[i].key[:4] not in ['TTYP','TFOR','TUNI','TDIS','TDIM','WCAX','1CTY',
                                          '2CTY','1CRP','2CRP','1CRV','2CRV','1CUN','2CUN',
                                          '1CDE','2CDE','1CTY','2CTY','1CDL','2CDL','11PC',
                                          '12PC','21PC','22PC']):
                hdu1.header.update(cards1[i].key, cards1[i].value, cards1[i].comment)
        outstr.append(hdu1)

# construct output mask bitmap extension

    if status == 0:
        hdu2 = ImageHDU(maskmap)
        for i in range(len(cards2)):
            if cards2[i].key not in hdu2.header.ascardlist().keys():
                hdu2.header.update(cards2[i].key, cards2[i].value, cards2[i].comment)
            else:
                hdu2.header.ascardlist()[cards2[i].key].comment = cards2[i].comment
        outstr.append(hdu2)

# construct principal component table

    if status == 0:
        cols = []
        for i in range(nmask):
            colname = 'PC' + str(i + 1)
            col = Column(name=colname,format='E',unit='e-/s',array=score[:,i])
            cols.append(col)
        hdu3 = new_table(ColDefs(cols))
        hdu3.header.update('EXTNAME','PRINCIPAL_COMPONENTS','name of extension')
        for i in range(nmask):
            hdu3.header.update('TTYPE' + str(i + 1),'PC' + str(i + 1),'column title: principal component number' + str(i + 1))
            hdu3.header.update('TFORM' + str(i + 1),'E','column format: float32')
            hdu3.header.update('TUNIT' + str(i + 1),'e-/s','column units: electrons per sec')
        outstr.append(hdu3)

# write output file

    if status == 0:
        outstr.writeto(outfile,checksum=True)

# close input structure

    if status == 0:
        status = kepio.closefits(instr,logfile,verbose)	    

# plotting defaults

    if status == 0:
        plotLatex = True
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
            plotLatex = False
    if status == 0:
        pylab.figure(figsize=[xsize,ysize])
        pylab.clf()

# clean up x-axis unit

    if status == 0:
	intime0 = float(int(tstart / 100) * 100.0)
	ptime = time + bjdref - intime0
	xlab = 'BJD $-$ %d' % intime0

# clean up y-axis units

    if status == 0:
        pout = copy(score)
	nrm = len(str(int(pout.max())))-1
	pout = pout / 10**nrm
	ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

# data limits

	xmin = ptime.min()
	xmax = ptime.max()
	ymin = pout.min()
	ymax = pout.max()
	xr = xmax - xmin
	yr = ymax - ymin

# plot window

        ax = pylab.axes([0.06,0.54,0.93,0.43])

# force tick labels to be absolute rather than relative

        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        pylab.setp(labels, 'rotation', 90)
        pylab.setp(pylab.gca(),xticklabels=[])

# plot principal components

        for i in range(nmask):
            pylab.plot(ptime,pout[:,i],linestyle='-',linewidth=lwidth)
        if not plotLatex:
            ylab = '10**%d electrons/sec' % nrm
        ylabel(ylab, {'color' : 'k'})
        grid()

# plot ranges

        pylab.xlim(xmin-xr*0.01,xmax+xr*0.01)
        pylab.ylim(ymin-yr*0.01,ymax+yr*0.01)

# plot output data

        ax = pylab.axes([0.06,0.09,0.93,0.43])

# force tick labels to be absolute rather than relative

        pylab.gca().xaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
        pylab.gca().yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))

# rotate y labels by 90 deg

        labels = ax.get_yticklabels()
        setp(labels, 'rotation', 90)

# clean up y-axis units

    if status == 0:
        pout = copy(pca_flux)
	nrm = len(str(int(pout.max())))-1
	pout = pout / 10**nrm
	ylab = '10$^%d$ e$^-$ s$^{-1}$' % nrm

# data limits

	ymin = pout.min()
	ymax = pout.max()
	yr = ymax - ymin
        ptime = numpy.insert(ptime,[0],[ptime[0]]) 
        ptime = numpy.append(ptime,[ptime[-1]])
        pout = numpy.insert(pout,[0],[0.0]) 
        pout = numpy.append(pout,0.0)

# plot time coadded principal component series

        pylab.plot(ptime[1:-1],pout[1:-1],color=lcolor,linestyle='-',linewidth=lwidth)
        pylab.fill(ptime,pout,color=fcolor,linewidth=0.0,alpha=falpha)
	pylab.xlabel(xlab, {'color' : 'k'})
        pylab.ylabel(ylab, {'color' : 'k'})
        pylab.grid()

# plot ranges

        pylab.xlim(xmin-xr*0.01,xmax+xr*0.01)
        if ymin >= 0.0: 
            pylab.ylim(ymin-yr*0.01,ymax+yr*0.01)
        else:
            pylab.ylim(1.0e-10,ymax+yr*0.01)

# render plot

        if cmdLine: 
            pylab.show()
        else: 
            pylab.ion()
            pylab.plot([])
            pylab.ioff()
	
# stop time

    if status == 0:
        kepmsg.clock('KEPPCA ended at',logfile,verbose)

    return

# -----------------------------------------------------------
# main

parfile = iraf.osfn("kepler$keppca.par")
t = iraf.IrafTaskFactory(taskname="keppca", value=parfile, function=keppca)
