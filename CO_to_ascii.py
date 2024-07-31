# -*- coding: utf-8 -*-
##############################################################################
#                                                                            #
##############################################################################

'''
Max Planck Institute for Extraterrestrial Physics
L-G Group

Nicolas Kurtovic
Contact: kurtovic at mpe.mpg.de

The following functions are used to support the extraction, handling, and
rewriting of visibility tables in the CASA environment. 

Some of these functions are adapted from SIMIO (Kurtovic 2024a).

If you use this file, please cite Kurtovic & Pinilla 2024.
https://ui.adsabs.harvard.edu/abs/2024A%26A...687A.188K/abstract

'''

import sys
import os
import numpy as np


##############################################################################

################################################################################
#                           FUNCTIONS FOR GAS HANDLING                         #
################################################################################


def _shift_uv(u, v, Re, Im, dRa, dDec):
    '''
    Applies spatial shift to the visibilities
    
    Args:
        - u: (numpy array) Array containing the ``u`` coordinates of the
                    visibilities.
        - v: (numpy array) Array containing the ``v`` coordinates of the
                    visibilities.
        - Re: (numpy array) Array containing the ``Real`` part of the
                    visibilities.
        - Im: (numpy array) Array containing the ``Imaginary`` part of the
                    visibilities.
        - dRa: (float) Value to shift the visiblities in RA, in **radians**.
        - dDec: (float) Value to shift the visiblities in Dec, in **radians**.
    '''
    # Calculate the value of the visibilities
    intensity = (Re + Im * 1j) * np.exp(2j * np.pi * (u * -dRa + v * -dDec))
    # Separate in real and imaginary part
    Re_shift = np.real(intensity)
    Im_shift = np.imag(intensity)
    # Return
    return u, v, Re_shift, Im_shift


def _rotate_uv(u, v, pa):
    '''
    Applies a rotation matrix to the visibilities, given a rotation angle
    of 'pa'.

    Args:
        - u: (numpy array) Array containing the ``u`` coordinates of the
                    visibilities.
        - v: (numpy array) Array containing the ``v`` coordinates of the
                    visibilities.
        - pa: (float) Value to rotate the visibilities, in **radians**.
    '''
    # Calculate rotation components
    cos_rot = np.cos(pa)
    sin_rot = np.sin(pa)
    # Rotate coordinates
    u_rot = u * cos_rot - v * sin_rot
    v_rot = u * sin_rot + v * cos_rot
    # Return
    return u_rot, v_rot


def _deproject_uv(u, v, inc):
    '''
    Applies an inverse inclination to the input visibilities.
    
    Args:
        - u: (numpy array) Array containing the ``u`` coordinates of the
                    visibilities.
        - v: (numpy array) Array containing the ``v`` coordinates of the
                    visibilities.
        - inc: (float) Value to incline the visibilities, in **radians**.
    '''
    # Apply inclination
    u_inc = u * np.cos(inc)
    # Return
    return u_inc, v


def _project_uv(u, v, inc):
    '''
    Applies an inclination to the input visibilities.
    
    Args:
        - u: (numpy array) Array containing the ``u`` coordinates of the
                    visibilities.
        - v: (numpy array) Array containing the ``v`` coordinates of the
                    visibilities.
        - inc: (float) Value to incline the visibilities, in **radians**.
    '''
    # Apply inclination
    u_inc = u / np.cos(inc)    
    # Return
    return u_inc, v


def handle_uv(u, v, Re, Im, inc, pa, dRa, dDec, inverse=False):
    '''
    Projects, deprojects, rotates, and shifts visibility points based on the
    input parameters.
    
    Args:
        - u: (numpy array) Array containing the ``u`` coordinates of the
                    visibilities.
        - v: (numpy array) Array containing the ``v`` coordinates of the
                    visibilities.
        - Re: (numpy array) Array containing the ``Real`` part of the
                    visibilities.
        - Im: (numpy array) Array containing the ``Imaginary`` part of the
                    visibilities.
        - inc: (float) Value to incline the visibilities, in **degrees**.
        - pa: (float) Value to rotate the visibilities, in **degrees**.
        - dRa: (float) Value to shift the visiblities in RA, in **arcsec**.
        - dDec: (float) Value to shift the visiblities in Dec, in **arcsec**.
        - inverse: (bool) Set ``False`` to deproject, set ``True`` to project.
                    Default: False.
    '''
    # Convert to radian
    pa   *= np.pi / 180.
    inc  *= np.pi / 180.
    dDec *= np.pi / (180. * 3600)
    dRa  *= np.pi / (180. * 3600)

    # If deprojection
    if not inverse:
        # Calculate shifted visibilities
        u, v, Re_shift, Im_shift = _shift_uv(u, v, Re, Im, dRa, dDec)
        # Correct rotation
        u_rot, v_rot = _rotate_uv(u, v, pa)
        # Correct inclination
        u_inc, v_inc = _deproject_uv(u_rot, v_rot, inc)
        # Return
        return u_inc, v_inc, Re_shift, Im_shift

    # If projection
    if inverse:
        # Apply inclination
        u_inc, v_inc = _project_uv(u, v, inc)
        # Apply rotation
        u_rot, v_rot = _rotate_uv(u_inc, v_inc, -pa)
        # Calculate shifted visibilities
        u_rot, v_rot, Re_shift, Im_shift = _shift_uv(u_rot, v_rot, Re, Im, dRa, dDec)
        # Return
        return u_rot, v_rot, Re_shift, Im_shift


def get_chanvel(ms_file, rest_freq):
    '''
    Returns channel velocities in m/s
    '''
    # Use CASA table tools to get frequencies
    tb.open(ms_file+"/SPECTRAL_WINDOW")
    freqs = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.close()
    
    vel = ((freqs / rest_freq) - 1) * qa.constants('c')['value']
    return -vel


def save_chan_info(name, chan_num, chan_vel):
    '''
    Writes the information of the channels in csv
    '''
    import csv
    aux_file = csv.writer(open('uvtables/info_channels.csv', 'w'))
    for key, val in dict(zip(list(chan_num), list(chan_vel))).items():
        aux_file.writerow([key, val])


def load_chan_info(name):
    '''
    Writes the information of the channels in csv
    '''
    import csv
    # Load velocities and channels
    chan_vel = []
    chan_num = []
    with open('uvtables/info_channels.csv') as aux_file:
        for i in csv.reader(aux_file):
            chan_num.append(i[0])
            chan_vel.append(i[1])
    return chan_num, chan_vel


def get_freqchan(ms_file):
    '''
    Gets the central frequency of each channel of the measurement set.
    The measurement set must have only one spw.
    '''
    # Use CASA table tools to get frequencies
    tb.open(ms_file+'/SPECTRAL_WINDOW')
    freqs = tb.getcol('CHAN_FREQ')
    tb.close()
    # Return
    return freqs


def change_freq(ms_file, over_freq):
    '''
    Gets the central frequency of each channel of the measurement set.
    The measurement set must have only one spw.
    '''
    # Use CASA table tools to get frequencies
    tb.open(ms_file+'/SPECTRAL_WINDOW', nomodify=False)
    tb.putcol('CHAN_FREQ', over_freq)
    tb.flush()
    tb.close()
    # Return
    return True


def ms_to_ascii(ms_cont,ascii_file,with_flags):
    '''
    Modified version of functions from Laura Perez to extract visibilities.
    '''
#    if with_flags: print "WARNING: flagged data may be included, unset with_flags=False."

    # Use CASA table tools to get columns of UVW, DATA, WEIGHT, etc.
    tb.open(ms_cont)
    data      = tb.getcol("DATA")
    uvw       = tb.getcol("UVW")
    weights   = tb.getcol("WEIGHT")
    wspectrum = tb.getcol("WEIGHT_SPECTRUM")
    ant1      = tb.getcol("ANTENNA1")
    ant2      = tb.getcol("ANTENNA2")
    flags     = tb.getcol("FLAG")
    tb.close()
    
    # Use CASA ms tools to get the channel/spw info
    ms.open(ms_cont)
    spw_info = ms.getspectralwindowinfo()
    nchan = spw_info["0"]["NumChan"]
    npol = spw_info["0"]["NumCorr"]
    ms.close()

    # Use CASA table tools to get frequencies
    tb.open(ms_cont+"/SPECTRAL_WINDOW")
    freqs = tb.getcol("CHAN_FREQ")
    tb.close()
    
    # convert spatial frequencies from m to lambda 
    uu = uvw[0,:]*freqs/qa.constants('c')['value']
    vv = uvw[1,:]*freqs/qa.constants('c')['value']
    Ruv = np.sqrt(uu[0]**2 + vv[0]**2)

    # check to see whether the polarizations are already averaged
    data = np.squeeze(data)
    wspectrum = np.squeeze(wspectrum)
    flags = np.squeeze(flags)
    
    if npol==1:
        Re = data.real
        Im = data.imag
        wspec = wspectrum

    elif npol==2:
        # polarization averaging
        Re_xx = data[0,:].real
        Re_yy = data[1,:].real
        Im_xx = data[0,:].imag
        Im_yy = data[1,:].imag
        if nchan==1:
            wspectrum_xx = weights[0,:]
            wspectrum_yy = weights[1,:]
        else:
            wspectrum_xx = wspectrum[0,:,:]
            wspectrum_yy = wspectrum[1,:,:]            
            flags = flags[0,:]*flags[1,:]
            # - weighted averages
        with np.errstate(divide='ignore', invalid='ignore'):
            Re = np.where((wspectrum_xx + wspectrum_yy) != 0, (Re_xx*wspectrum_xx + Re_yy*wspectrum_yy) / (wspectrum_xx + wspectrum_yy), 0.)
            Im = np.where((wspectrum_xx + wspectrum_yy) != 0, (Im_xx*wspectrum_xx + Im_yy*wspectrum_yy) / (wspectrum_xx + wspectrum_yy), 0.)
            wspec = (wspectrum_xx + wspectrum_yy)

    elif npol==4:
        print ('npol is '+str(npol))
        # polarization averaging
        Re_xx = data[0,:].real
        Re_xy = data[1,:].real
        Re_yx = data[2,:].real
        Re_yy = data[3,:].real
        Im_xx = data[0,:].imag
        Im_xy = data[1,:].imag
        Im_yx = data[2,:].imag
        Im_yy = data[3,:].imag
        if nchan==1:
            wspectrum_xx = weights[0,:]
            wspectrum_xy = weights[1,:]
            wspectrum_yx = weights[2,:]
            wspectrum_yy = weights[3,:]
        else:
            wspectrum_xx = weights[0,:,:]
            wspectrum_xy = weights[1,:,:]
            wspectrum_yx = weights[2,:,:]
            wspectrum_yy = weights[3,:,:]
            flags = flags[0,:]*flags[1,:]*flags[2,:]*flags[3,:]
            # - weighted averages
        with np.errstate(divide='ignore', invalid='ignore'):
            #wspec = wspectrum_xx + wspectrum_xy + wspectrum_yx + wspectrum_yy
            wspec = wspectrum_xx + wspectrum_yy
            zero_cond = (wspec != 0)
            #Re = np.where(zero_cond, (Re_xx*wspectrum_xx + Re_xy*wspectrum_xy + \
            #                          Re_yx*wspectrum_yx + Re_yy*wspectrum_yy) / wspec, 0.)
            #Im = np.where(zero_cond, (Im_xx*wspectrum_xx + Im_xy*wspectrum_xy + \
            #                          Im_yx*wspectrum_yx + Im_yy*wspectrum_yy) / wspec, 0.)
            Re = np.where(zero_cond, (Re_xx*wspectrum_xx + Re_yy*wspectrum_yy) / wspec, 0.)
            Im = np.where(zero_cond, (Im_xx*wspectrum_xx + Im_yy*wspectrum_yy) / wspec, 0.)

    # toss out the autocorrelation placeholders
    xc = np.where(ant1 != ant2)[0]

    # check if there's only a single channel
    # else, make 1-d arrays since channel info is not important now:
    if nchan==1:
        data_real = Re[xc]
        data_imag = Im[xc]
        data_flags = flags[:,xc]
        data_wspec = wspec[xc]
        data_uu = uu[:,xc].flatten()
        data_vv = vv[:,xc].flatten()
    else: 
        data_real = Re[:,xc].flatten()
        data_imag = Im[:,xc].flatten()
        data_flags = flags[:,xc].flatten()
        data_wspec = wspec[:,xc].flatten()
        data_uu = uu[:,xc].flatten()
        data_vv = vv[:,xc].flatten()

    # remove flagged data if user wants to:
    if with_flags == False:
        if np.any(data_flags): 
            print ("flagged data not included!")
            data_uu = data_uu[np.logical_not(data_flags)]
            data_vv = data_vv[np.logical_not(data_flags)]
            data_real = data_real[np.logical_not(data_flags)]
            data_imag = data_imag[np.logical_not(data_flags)]
            data_wspec =data_wspec[np.logical_not(data_flags)]
            
    # Write the ascii file (Format will be: u,v,Re,Im,We)
    out_file = ascii_file
    np.savetxt(out_file,np.transpose([data_uu, data_vv, \
                                      data_real, data_imag, data_wspec]), \
                                     fmt="%+.15e")


def ascii_to_ms(ascii_file, ms_file, new_ms_file):
    '''
    Modified version of functions from Laura Perez to insert visibilities.
    '''

    # Read ascii file
    model_uu, model_vv, \
    model_real, model_imag, model_wspec = np.loadtxt(ascii_file,unpack=True)

    # Use CASA table tools to get columns of UVW, DATA, WEIGHT, etc.
    tb.open(ms_file)
    data      = tb.getcol("DATA")
    ant1      = tb.getcol("ANTENNA1")
    ant2      = tb.getcol("ANTENNA2")
    tb.close()

    # Use CASA ms tools to get the channel/spw info
    ms.open(ms_file)
    spw_info = ms.getspectralwindowinfo()
    nchan = spw_info["0"]["NumChan"]
    npol = spw_info["0"]["NumCorr"]
    ms.close()
    
    # Figure out if there is autocorrelation data
    ac = np.where(ant1 == ant2)[0] # autocorrelation
    xc = np.where(ant1 != ant2)[0] # correlation
    
    # Check that ascii_file data has the same shape as the ms_file data
    if len(model_real) != nchan*len(xc): 
        print ("Problem found, array dimensions do not match!")
        print ('Will not create a model ms file, ' + \
               'probably need to include flagged data.')
        return

    # New data has to have the same shape as original data
    data_array = np.zeros((npol, nchan, ant1.shape[0])).astype(complex)
    
    # Put the model data on both polarizations
    data_array[0, :, xc] = model_real.reshape(nchan,len(xc)).transpose() + 1j * model_imag.reshape(nchan,len(xc)).transpose()
    data_array[1, :, xc] = model_real.reshape(nchan,len(xc)).transpose() + 1j * model_imag.reshape(nchan,len(xc)).transpose()
    # There is no model for the autocorrelation data
    data_array[0, :, ac] = 0 + 0j
    data_array[1, :, ac] = 0 + 0j

    # Duplicate ms_file and remove any extra columns
    split(vis=ms_file, outputvis=new_ms_file, datacolumn='data')
    tb.open(new_ms_file, nomodify=False)
    tb.putcol("DATA", data_array)
    tb.flush()
    tb.close()

    print ('New ms file is: ' + new_ms_file)

