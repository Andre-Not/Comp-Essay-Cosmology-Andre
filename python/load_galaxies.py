"""A collection of dictionaries containing data for various galaxies.

This module exists to support other modules by organizing data in a quasi-standardized manner.
"""
#Data to fit to for each galaxy to be used in workshop

###############################
########## Imports ############
###############################

import dataPython         as dp
import numpy              as np
import scipy.interpolate  as inter

###############################
########### NGC5533 ###########
###############################

NGC5533 = {
    # Load data from files for Noordermeer 2008 band and fitted curves
    
    # 'raw' in the sense that I haven't split everything up. Idk there's probably a better way to name these
    'raw_total'        : dp.getXYdata('NGC5533/noord-120kpc-total.txt'     ),
    'raw_blackhole'    : dp.getXYdata('NGC5533/noord-120kpc-blackhole.txt' ),
    'raw_bulge'        : dp.getXYdata('NGC5533/noord-120kpc-bulge.txt'     ),
    'raw_disk'         : dp.getXYdata('NGC5533/noord-120kpc-disk.txt'      ),
    'raw_halo'         : dp.getXYdata('NGC5533/noord-120kpc-halo.txt'      ),
    'raw_gas'          : dp.getXYdata('NGC5533/noord-120kpc-gas.txt'       ),
    'raw_band_btm'     : dp.getXYdata('NGC5533/noord-120kpc-bottomband.txt'),
    'raw_band_top'     : dp.getXYdata('NGC5533/noord-120kpc-topband.txt'   ),
    
    # Get data from 100kpc file
    'measured_data'    : dp.getXYdata_wXYerr('NGC5533/100kpc_data.txt'),
    
    # Some constants
    'i'                : 52,    # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)  
    'D_Mpc'            : 54.3,  # Distance [Mpc]
    'Mabs'             : -21.66,# B-band absolute magnitude of NGC 5533 (Norrdermeer, Van Der Hulst, 2007)
    'check'            : True
}
"""Data for galaxy NGC5533. Parameters and data measurements in this dictionary are from [Noordermeer2007]_ unless noted otherwise.

:keys:
    blackhole: [dict] Further information pertaining to the black hole component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the black hole component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the black hole component.
        
        Mbh: [float] Mass of the black hole.
        
        r: [array] Radial (kpc) values of the black hole component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the black hole component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the black hole component.

        v: [array] Velocity (km/s) values of the black hole component.
    
    bulge: [dict] Further information pertaining to the bulge component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bulge component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the bulge component.

        Lb: [float] Luminosity (solar luminosities) of the bulge.

        ML: [float] Mass-to-light ratio of the bulge.

        n: [float] Concentration parameter.

        q: [float] Intrinsic axis ratio.

        r: [array] Radial (kpc) values of the bulge component.

        re_arcsec: [float] Effective radius (arcsec).

        re_rad: [float] Effective radius (radians).

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the bulge component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the bulge component.

        v: [array] Velocity (km/s) values of the bulge component.

    check: [bool] A flag indicating that the necessary data for the "Literature Search Check" in `07_Bonus_Bulge_Rotation_Curve.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_ is available.

    disk: [dict] Further information pertaining to the disk component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the disk component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the disk component.
        
        r: [array] Radial (kpc) values of the disk component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the disk component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the disk component.

        v: [array] Velocity (km/s) values of the disk component.

    D_Mpc: [float] The distance to the galaxy in Megaparsecs.

    D_kpc: [float] The distance to the galaxy in kiloparsecs.

    galaxyname: [string] Name of the galaxy.

    gas: [dict] Further information pertaining to the gas component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the gas component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the gas component.
        
        r: [array] Radial (kpc) values of the gas component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the gas component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the gas component.

        v: [array] Velocity (km/s) values of the gas component.

    halo: [dict] Further information pertaining to the halo component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the halo component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the halo component.
        
        r: [array] Radial (kpc) values of the halo component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the halo component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the halo component.

        v: [array] Velocity (km/s) values of the halo component.

    i: [int] The inclination angle, in degrees, of the galaxy.

    m_r_errors: [array] Errors of measured radii (kpc).

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).

    Mabs: [float] The absolute magnitude of the galaxy.

    massBH: [float] Mass of central black hole (solar masses).

    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ex: [list] Float values representing the error on the radial (kpc) data.

        ey: [list] Float values representing the error on the velocity (km/s) data.

    n_band_btm: [scipy.interpolate._bsplines.BSpline] Spline of the bottom side of the band.

    n_band_top: [scipy.interpolate._bsplines.BSpline] Spline of the top side of the band.

    n_cb: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_ct: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the top side of the error band.

    n_kb: [int] Degree of spline returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_kt: [int] Degree of spline returned by using scipy.interpolate.splrep on the top side of the error band.

    n_r_btmband: [array] Radial values of the bottom side of the error band.

    n_r_topband: [array] Radial values of the top side of the error band.

    n_tb: [array] Vector of knots returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_tt: [array] Vector of knots returned by using scipy.interpolate.splrep on the top side of the error band.

    n_v_bandwidth: [array] Width of the error band, in km/s.

    n_v_btmband: [array] Velocity values of the bottom side of the error band.

    n_v_btmband: [array] Velocity values of the top side of the error band.

    raw_band_btm: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_band_top: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_blackhole: [dict] Data representing the black hole's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_bulge: [dict] Data representing the bulge's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_disk: [dict] Data representing the disk's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_gas: [dict] Data representing the gas' contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_halo: [dict] Data representing the halo's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.    
    
    raw_total: [dict] Data representing the total theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    rc: [float] Core radius (kpc).

    rho0: [float] Central mass density (solar mass/kpc^3).

    sources: [dict] Printing helper variables indicating units and sources for certain variables.
        D: [string] Units and source for distance to galaxy.
        
        i: [string] Units and source for inclination angle.

        Lb: [string] Units and source for bulge luminosity.
        
        Mabs: [string] Units and source for absolute magnitude.

        ML: [string] Units and source for bulge mass-light ratio.

        n: [string] Units and source for concentration parameter.

        re: [string] Units and source for effective radius.

        q: [string] Units and source for intrinsic axis ratio.

    total: [dict] Further information pertaining to the total theoretical rotation curve.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the total curve.
        
        k: [int] Degree of the spline returned by using scipy.interpolate.splrep on the total curve.

        r: [array] Radial (kpc) values of the total curve.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the total curve.

        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the total curve.

        v: [array] Velocity (km/s) values of the total curve.

"""

# Parameters ########################
NGC5533['galaxyname'] = 'NGC 5533'   # NGC catalog number of the galaxy
NGC5533['rho0'] = 0.31e9       # central mass density (in solar mass/kpc^3), Source: Noordermeer (2007)     
NGC5533['rc'] = 1.4            # core radius (in kpc), Source: Noordermeer (2007)
NGC5533['massbh'] = 2.7e9      # mass of central black hole (in solar masses), Source: Noordermeer (2007)

#Organize 100kpc data
NGC5533['m_radii']      = np.asarray(NGC5533['measured_data']['xx'])
NGC5533['m_velocities'] = np.asarray(NGC5533['measured_data']['yy'])
NGC5533['m_r_errors']   = np.asarray(NGC5533['measured_data']['ex'])
NGC5533['m_v_errors']   = np.asarray(NGC5533['measured_data']['ey'])

#Organize band data
NGC5533['n_r_btmband']   = np.asarray(NGC5533['raw_band_btm']['xx'])
NGC5533['n_v_btmband']   = np.asarray(NGC5533['raw_band_btm']['yy'])
NGC5533['n_r_topband']   = np.asarray(NGC5533['raw_band_top']['xx'])
NGC5533['n_v_topband']   = np.asarray(NGC5533['raw_band_top']['yy'])
NGC5533['n_v_bandwidth'] = (NGC5533['n_v_topband'] - NGC5533['n_v_btmband'])/2
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][0::28] #For weights, v_errors and band must line up.
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][1:]
    
# Smoothing
NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'] = inter.splrep(NGC5533['n_r_btmband'],NGC5533['n_v_btmband'])
NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'] = inter.splrep(NGC5533['n_r_topband'],NGC5533['n_v_topband'])
NGC5533['n_band_btm'] = inter.BSpline(NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'])
NGC5533['n_band_top'] = inter.BSpline(NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'])

# Total Curve #######################
NGC5533['total'] = {
    'r' : np.asarray(NGC5533['raw_total']['xx']),
    'v' : np.asarray(NGC5533['raw_total']['yy'])
}
NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'] = inter.splrep(NGC5533['total']['r'], NGC5533['total']['v'])
NGC5533['total']['spline'] = inter.BSpline(NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'])

# Black Hole ########################
NGC5533['blackhole'] = {
    'r'  : np.asarray(NGC5533['raw_blackhole']['xx']),
    'v'  : np.asarray(NGC5533['raw_blackhole']['yy']),
    'Mbh': 2.7e9
}
NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'] = inter.splrep(NGC5533['blackhole']['r'], NGC5533['blackhole']['v'])
NGC5533['blackhole']['spline'] = inter.BSpline(NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'])

# Bulge #############################
NGC5533['bulge'] = {
    'r' : np.asarray(NGC5533['raw_bulge']['xx']),
    'v' : np.asarray(NGC5533['raw_bulge']['yy']),
    'n'      : 2.7,     # Concentration parameter [unitless] (Noordermeer, Van Der Hulst, 2007)
    'q'      : 0.33,    # Intrinsic axis ratio [unitless] (Noordermeer, 2008)
    're_arcsec' : 9.9,     # Effective radius [arcsec] (Noordermeer, Van Der Hulst, 2007)
    'Lb'     : 3.6e10,  # Bulge luminosity [solar luminosity] (Calculated from absolute magnitude of the bulge)
    'ML'     : 2.8      # Mass-to-light ratio of bulge [unitless] (Noordermeer, 2008)
}
NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'] = inter.splrep(NGC5533['bulge']['r'], NGC5533['bulge']['v'])
NGC5533['bulge']['spline'] = inter.BSpline(NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'])
NGC5533['bulge']['re_rad'] = NGC5533['bulge']['re_arcsec'] * (np.pi / (3600*180)) #arcsec to rad
NGC5533['D_kpc'] = NGC5533['D_Mpc'] * 1000 #Mpc to kpc
NGC5533['bulge']['re_kpc'] = NGC5533['bulge']['re_rad'] * NGC5533['D_kpc'] #Effective radius kpc
NGC5533['sources'] = {
    'Mabs': ' [unitless] (Source #2, Table A4)',
    'n'   : ' [unitless] (Source #2, Table A4)',
    'q'   : ' [unitless] (Source #1, Table 1)',
    'i'   : ' [degrees] (Source #2, Table A2)',
    're'  : ' [kpc] (Source #2, Table A4)',
    'Lb'  : ' [Lsun] (Source #2, calculated from absolute magnitude)',
    'ML'  : ' [unitless] (Source #1, Table 1)',
    'D'   : ' [Mpc] (Source #2, Table 1)'
}

# Disk ##############################
NGC5533['disk'] = {
    'r' : np.asarray(NGC5533['raw_disk']['xx']),
    'v' : np.asarray(NGC5533['raw_disk']['yy'])
}
NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'] = inter.splrep(NGC5533['disk']['r'], NGC5533['disk']['v'])
NGC5533['disk']['spline'] = inter.BSpline(NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'])

# Halo ##############################
NGC5533['halo'] = {
    'r' : np.asarray(NGC5533['raw_halo']['xx']),
    'v' : np.asarray(NGC5533['raw_halo']['yy'])
}
NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'] = inter.splrep(NGC5533['halo']['r'], NGC5533['halo']['v'])
NGC5533['halo']['spline'] = inter.BSpline(NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'])

# Gas ###############################
NGC5533['gas'] = {
    'r' : np.asarray(NGC5533['raw_gas']['xx']),
    'v' : np.asarray(NGC5533['raw_gas']['yy'])
}
NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'] = inter.splrep(NGC5533['gas']['r'], NGC5533['gas']['v'])
NGC5533['gas']['spline'] = inter.BSpline(NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'])

