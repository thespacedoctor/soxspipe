
from astropy.nddata import CCDData
from matplotlib import pyplot as plt
from astropy import units as u
from fundamentals.logs import emptyLogger
from copy import copy
from astropy.stats import sigma_clipped_stats

frame = CCDData.read('/Volumes/My Passport/SOXS/nir/raw_frames/SOXS_SLT_WAVE_NIR_122_0034.fits', hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                            hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

from soxspipe.commonutils.toolkit import quicklook_image

log = emptyLogger()
#quicklook_image(log=log, CCDObject=frame, show=True, ext='data', stdWindow=3, title="Background Light", surfacePlot=False)

mean, median, std = sigma_clipped_stats(frame, sigma=50.0, stdfunc="mad_std", cenfunc="median", maxiters=3)



#Reading the content of the fits table in resources/static_calibrations/soxs/SOXS_MPH_TAB_NIR_Xe.fits
from astropy.table import Table

table = Table.read('./soxspipe/resources/static_calibrations/soxs/SOXS_SPH_TAB_NIR_Xe.fits', hdu=1)

palette = copy(plt.cm.viridis)
palette.set_bad("#dc322f", 1.0)
stdWindow=10
vmax = median + stdWindow * 0.5 * std
vmin = median - stdWindow * 0.5 * std
print(table)
fig = plt.figure(figsize=(12, 5))




plt.imshow(frame.data, origin='lower', vmin=vmin, vmax=vmax, cmap=palette, alpha=1, aspect='auto')
plt.scatter(table['detector_x'],table['detector_y'], s=10, c='red', marker='x')
#plt.scatter(table['detector_x']+0.0,table['detector_y']-0.0, s=10, c='red', marker='x')
plt.show()


#table['detector_x'] = table['detector_x'] - 1.5
#table['detector_y'] = table['detector_y'] - 4.0

#table.write('./soxspipe/resources/static_calibrations/soxs/SOXS_SPH_TAB_NIR_Xe.fits', format='fits', overwrite=True)


