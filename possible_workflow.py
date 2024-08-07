#!/usr/bin/env python
# coding: utf-8

# Imports


import settings
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from lasairmod import LasairError, lasair_client as lasair
import sys, time
import json
import scipy
from sympy import *
from ztfquery import lightcurve, query, metasearch
import scipy.stats as st
from scipy.stats import ks_2samp
from pandas.core.computation.check import NUMEXPR_INSTALLED
import csv
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import glob



# Defining paths and creating directories

path = os.getcwd()

directories = ['Lightcurves', 'CMDs', 'AlertStats']
for directory in directories:         # create output directories if not present
    if not os.path.exists(directory):
        os.makedirs(directory)



# Defining Functions


def diff_to_app_mag(ref_mag, alert_mag, sign): 
   
    """Function converts from difference magnitudes to apparent magnitudes. 
    
    Parameters
    ----------
    ref_mag : float
        Magnitude of reference source.
    alert_mag : float
        Difference magnitude of alert.
    sign : str
        Sign of the difference magnitude where
        t = positive difference magnitude = brightening
        f = negative difference magnitude = faintening.
        

    Returns
    -------
    app_mag
        Apparent magnitude as type ``float``.
    
    """
    if sign =='t': 
        app_mag = ((-2.5)*log((10**((-.4)*ref_mag)) + (10**((-.4)*alert_mag)),10))
    else:
        app_mag = ((-2.5)*log((10**((-.4)*ref_mag)) - (10**((-.4)*alert_mag)),10)) 
    
    return float(app_mag)




def app_mag_err_prop(ref_mag, alert_mag, sign, ref_err, alert_err):
    
    """Function propagates error when converting from difference magnitudes to apparent magnitudes. 
    
    Parameters
    ----------
    ref_mag : float
        Magnitude of reference source.
    
    alert_mag : float
        Difference magnitude of alert.
    
    sign : str
        Sign of the difference magnitude where
        t = positive difference magnitude = brightening
        f = negative difference magnitude = faintening.
    
    ref_err : float
        1-sigma uncertainty in reference magnitude.
   
    alert_err : float
        1-sigma uncertainty in alert magnitude.
        

    Returns
    -------
    app_err
        Apparent magnitude error as type ``float``.
    
    """
    
    
    r, a= symbols('r a', real=True) # Defines variables r (reference mag) and a (alert mag)
    
    if sign =='t': 
        f = ((-2.5)*log((10**((-.4)*r)) + (10**((-.4)*a)),10))
    else:
        f = ((-2.5)*log((10**((-.4)*r)) - (10**((-.4)*a)),10)) 
    
    
    # f = diff_to_app_mag(ref_mag, alert_mag, sign)
    
    d_ref = diff(f, r).subs([(a, alert_mag), (r, ref_mag)])         # Partial derivative of f with respect to r
    d_alert = diff(f, a).subs([(a, alert_mag), (r, ref_mag)])       # Partial derivative of f with respect to a
    
    app_err = float(sqrt((((d_ref)**2)*((ref_err)**2))+(((d_alert)**2)*((alert_err)**2))))   # Error propagation formula
    
    return app_err


# Query for Disappearing Stars


mjdnow = str(Time.now().mjd)
jdnow = str(Time.now().jd)
days = str(1) # number of days 

selected = '*'


tables = 'objects,sherlock_classifications'


conditions = """
objects.objectId=sherlock_classifications.objectId
AND (objects.sgscore1 > 0.9)
AND (sherlock_classifications.classification != "SN")
AND (sherlock_classifications.classification != "NT")
AND (sherlock_classifications.classification != "AGN")
AND (objects.ncand >= 10)
AND (sherlock_classifications.catalogue_table_name LIKE "%gaia%")
AND (objects.objectId LIKE "ZTF24%")
AND sherlock_classifications.separationArcsec < 0.5
AND ((objects.sgmag1 < 16)
   OR (objects.srmag1 < 16))
AND ISNULL(objects.ncandgp)
AND ("""+jdnow+"""- objects.jdmax) < """+days+"""

"""

L = lasair(settings.API_TOKEN, endpoint = "https://lasair-ztf.lsst.ac.uk/api")

try:
    v4 = L.query(selected, tables, conditions)
except LasairError as e:
    print(e)
    
print('Query returned ' + str(len(v4)) + ' candidates in the past '+ days + ' days')


# Pulling info from Gaia catalog of nearby stars

#reading csv file
df = pd.read_csv(path + '/gaia.tsv',comment="#",delimiter=";")

#dropping empty rows and rows with units rather than numbers
gaia = df.drop([0,1])

#calculating BP-RP color and absolute G magnitude for Gaia catalog of nearby stars

BP = pd.to_numeric(gaia['BPmag'], errors='coerce') #BP magnitude
RP = pd.to_numeric(gaia['RPmag'], errors='coerce') #RP magnitude
color = BP - RP #BP-RP color
gaia_gmag = pd.to_numeric(gaia['Gmag'], errors='coerce') #apparent g magnitudes
plx = pd.to_numeric(gaia['Plx'], errors= 'coerce')/1000 #parallax in arcsec
d = 1/plx #distance in parsecs
gaia_GMAG = gaia_gmag - 5*(np.log10(d/10)) #absolute g magnitudes using distance modulus and parallax


#Pulling latest csv file 


# List all files in the directory
CSVs = glob.glob(path + '/AlertStats/*')

# Sort files lexicographically
CSVs.sort()

# Get the lexicographically greatest file
latest_CSV = CSVs[-1] if len(CSVs) > 0 else None

# read the file if it exists
if latest_CSV:
    latest_stats = pd.read_csv(latest_CSV, delimiter=",")


# Plotting Light Curves from Alert Packet and ZTF archive 


# Pulling alert packets from stars that passed Lasair filter

Dips= L.objects([row['objectId'] for row in v4])
lc = {}
candidates = []
trash = []
cand_color = []
trash_color = []
cand_gmag = []
trash_gmag = []

broker = path+ '/AlertStats/stats_' + mjdnow+ '.csv'

headers = ["Object", "g_X2", "r_X2", "i_X2", "g_5-95_X2", "r_5-95_X2", "i_5-95_X2", "g_KSpvalue", "r_KSpvalue", "i_KSpvalue", 
           "g_depth", "r_depth", "i_depth", "nhist", "g_nhist", "r_nhist", "i_nhist", "nalert", "g_nalert", "r_nalert", "i_nalert", 
           "ramean", "decmean", "gaiara", "gaiadec", "g_med", "r_med", "i_med", "disc_mjd", "latest_mjd", "gaia_sourceid", 
           "gaia_app_gmag", "gaia_plx", "gaia_abs_gmag", "gaia_BP-RP", "num"]

units = ["", "", "", "", "", "", "", "", "", "", 
           "mag", "mag", "mag", "", "", "", "", "", "", "", "", 
           "deg", "deg", "deg", "deg", "mag", "mag", "mag", "", "", "", 
           "mag", "mas", "mag", "mag", ""]   
    
    
with open(broker, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(headers)
    writer.writerow(units)
    

for obj in Dips:
    try:
        lc[obj['objectId']] = {'candidates': obj['candidates']}
    except TypeError:
        print('no')
    
    # Creating dictionary of variables
    variables = {}
    
    # Pulling historic ZTF light curve
    ztf = lightcurve.LCQuery.from_position(obj['objectData']['ramean'], obj['objectData']['decmean'], 1)
        
    # Defining arrays of ZTF historic magnitudes for each filter
    ztf_gsamp = np.array(ztf.data['mag'][ztf.data['filtercode'] == 'zg'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_rsamp = np.array(ztf.data['mag'][ztf.data['filtercode'] == 'zr'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_isamp = np.array(ztf.data['mag'][ztf.data['filtercode'] == 'zi'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    
    # same for magnitude errors
    ztf_gerr = np.array(ztf.data['magerr'][ztf.data['filtercode'] == 'zg'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_rerr = np.array(ztf.data['magerr'][ztf.data['filtercode'] == 'zr'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_ierr = np.array(ztf.data['magerr'][ztf.data['filtercode'] == 'zi'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    
    # ...and dates
    ztf_gmjd = np.array(ztf.data['mjd'][ztf.data['filtercode'] == 'zg'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_rmjd = np.array(ztf.data['mjd'][ztf.data['filtercode'] == 'zr'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    ztf_imjd = np.array(ztf.data['mjd'][ztf.data['filtercode'] == 'zi'][ztf.data['mjd'] < obj['objectData']['discMjd']][ztf.data['magerr'] > 0])
    
    # calculating 5th - 95th percentile ranges for historic data 
    if len(ztf_gsamp) > 0:
        variables['ghist'] = ztf_gsamp[(obj['objectData']['discMjd'] - ztf_gmjd) > 100] # restricting data to 100 days before discovery date
        variables['ghisterr'] = ztf_gerr[(obj['objectData']['discMjd'] - ztf_gmjd) > 100]
        variables['g_p5'] = np.percentile(variables['ghist'], 5)
        variables['g_p95'] = np.percentile(variables['ghist'], 95)
        variables['g_5to95'] = variables['ghist'][(variables['ghist'] > variables['g_p5']) & (variables['ghist'] < variables['g_p95'])]
        variables['gerr_5to95'] = variables['ghisterr'][(variables['ghist'] > variables['g_p5']) & (variables['ghist'] < variables['g_p95'])]
    
    if len(ztf_rsamp) > 0:
        variables['rhist'] = ztf_rsamp[(obj['objectData']['discMjd'] - ztf_rmjd) > 100] # restricting data to 100 days before discovery date
        variables['rhisterr'] = ztf_rerr[(obj['objectData']['discMjd'] - ztf_rmjd) > 100]
        variables['r_p5'] = np.percentile(variables['rhist'], 5)
        variables['r_p95'] = np.percentile(variables['rhist'], 95)
        variables['r_5to95'] = variables['rhist'][(variables['rhist'] > variables['r_p5']) & (variables['rhist'] < variables['r_p95'])]
        variables['rerr_5to95'] = variables['rhisterr'][(variables['rhist'] > variables['r_p5']) & (variables['rhist'] < variables['r_p95'])]
    
    if len(ztf_isamp) > 0:
        variables['ihist'] = ztf_isamp[(obj['objectData']['discMjd'] - ztf_imjd) > 100] # restricting data to 100 days before discovery date
        variables['ihisterr'] = ztf_ierr[(obj['objectData']['discMjd'] - ztf_imjd) > 100]
        variables['i_p5'] = np.percentile(variables['ihist'], 5)
        variables['i_p95'] = np.percentile(variables['ihist'], 95)
        variables['i_5to95'] = variables['ihist'][(variables['ihist'] > variables['i_p5']) & (variables['ihist'] < variables['i_p95'])]
        variables['ierr_5to95'] = variables['ihisterr'][(variables['ihist'] > variables['i_p5']) & (variables['ihist'] < variables['i_p95'])]
    
    # Creating empty lists for each ZTF filter-band to be filled with apparent magnitudes and errors from alert packet
    gmag = []
    rmag = []
    imag = []
    
    gerr = []
    rerr = []
    ierr = []
    
    gmjd = []
    rmjd = []
    imjd = []
    alertmjd = []
    
    for alert in lc[obj['objectId']]['candidates']:
        if 'isdiffpos' in alert:
            alertmjd.append(alert['mjd'])
            
            # appending magnitudes and errors in each filter to their respective lists
            if alert['fid']== 1:
                variables['app_gmag'] = diff_to_app_mag(np.median(variables['ghist']), alert['magpsf'], alert['isdiffpos']) # calculating apparent magnitude 
                variables['app_gerr'] = app_mag_err_prop(np.median(variables['ghist']), alert['magpsf'], alert['isdiffpos'], np.median(variables['ghisterr']), alert['sigmapsf']) # propagating error      
                gmag.append(variables['app_gmag'])
                gerr.append(variables['app_gerr'])
                gmjd.append(alert['mjd'])
            
            if alert['fid']== 2:
                variables['app_rmag'] = diff_to_app_mag(np.median(variables['rhist']), alert['magpsf'], alert['isdiffpos']) # calculating apparent magnitude 
                variables['app_rerr'] = app_mag_err_prop(np.median(variables['rhist']), alert['magpsf'], alert['isdiffpos'], np.median(variables['rhisterr']), alert['sigmapsf']) # propagating error
                rmag.append(variables['app_rmag'])
                rerr.append(variables['app_rerr'])
                rmjd.append(alert['mjd'])
            
            if alert['fid']== 3:
                variables['app_imag'] = diff_to_app_mag(np.median(variables['ihist']), alert['magpsf'], alert['isdiffpos']) # calculating apparent magnitude 
                variables['app_ierr'] = app_mag_err_prop(np.median(variables['ihist']), alert['magpsf'], alert['isdiffpos'], np.median(variables['ihisterr']), alert['sigmapsf']) # propagating error
                imag.append(variables['app_imag'])
                ierr.append(variables['app_ierr'])
                imjd.append(alert['mjd'])
    
    
    # Restricting alert epoch magnitudes to data observed at least 30 days after the discovery date in Lasair
    alert_gsamp = np.array(gmag)[(np.array(gmjd) - obj['objectData']['discMjd']) >= 30]
    alert_rsamp = np.array(rmag)[(np.array(rmjd) - obj['objectData']['discMjd']) >= 30]
    alert_isamp = np.array(imag)[(np.array(imjd) - obj['objectData']['discMjd']) >= 30]
    
    # Performing K-S test with 95% confidence interval and calculating chi square statistic for each filter
    if len(alert_gsamp) > 0 and len(variables['ghist']) > 0:
        variables['KS_g'] = st.ks_2samp(alert_gsamp, variables['ghist'])
        variables['g_pvalue'] = variables['KS_g'].pvalue
        variables['chisq_g'] = np.sum((((variables['ghist'] - np.median(variables['ghist'])) / variables['ghisterr'])**2 ) / float(len(variables['ghist'])))
        variables['chisq_g_5to95'] = np.sum((((variables['g_5to95'] - np.median(variables['g_5to95'])) / variables['gerr_5to95'])**2 ) / float(len(variables['g_5to95'])))
        variables['g_med'] = np.median(variables['ghist'])
        variables['g_depth'] = np.median(alert_gsamp) - np.median(variables['ghist'])
        
    if len(alert_rsamp) > 0 and len(variables['rhist']) > 0:  
        variables['KS_r'] = st.ks_2samp(alert_rsamp, variables['rhist'])
        variables['r_pvalue'] = variables['KS_r'].pvalue
        variables['chisq_r'] = np.sum((((variables['rhist'] - np.median(variables['rhist'])) / variables['rhisterr'])**2 ) / float(len(variables['rhist'])))
        variables['chisq_r_5to95'] = np.sum((((variables['r_5to95'] - np.median(variables['r_5to95'])) / variables['rerr_5to95'])**2 ) / float(len(variables['r_5to95'])))
        variables['r_med'] = np.median(variables['rhist'])
        variables['r_depth'] = np.median(alert_rsamp) - np.median(variables['rhist'])
        
    if len(alert_isamp) > 0 and len(variables['ihist']) > 0: 
        variables['KS_i'] = st.ks_2samp(alert_isamp, variables['ihist'])
        variables['i_pvalue'] = variables['KS_i'].pvalue
        variables['chisq_i'] = np.sum((((variables['ihist'] - np.median(variables['ihist'])) / variables['ihisterr'])**2 ) / float(len(variables['ihist'])))
        variables['chisq_i_5to95'] = np.sum((((variables['i_5to95'] - np.median(variables['i_5to95'])) / variables['ierr_5to95'])**2 ) / float(len(variables['i_5to95'])))
        variables['i_med'] = np.median(variables['ihist'])
        variables['i_depth'] = np.median(alert_isamp) - np.median(variables['ihist']) 
   
    # Imposing conditions that K-S test p-value is < 0.05 (95% confidence interval) and Chi square statistic is < 50
    g = True
    r = True
    i = True

    # if g true then there are alert samples, at least a hundred historic samples, and they meet the stats
    if len(gmag) > 0:
        g = False
        if len(alert_gsamp) > 0 and len(variables['ghist']) >= 100:
            g = (variables['KS_g'].pvalue < 0.05) and (np.median(alert_gsamp) > variables['g_p95']) and (variables['chisq_g'] < 15) and (variables['chisq_g_5to95'] < 5) 
    
    if len(rmag) > 0:
        r = False
        if len(alert_rsamp) > 0 and len(variables['rhist']) >= 100:
            r = (variables['KS_r'].pvalue < 0.05) and (np.median(alert_rsamp) > variables['r_p95']) and (variables['chisq_r'] < 15) and (variables['chisq_r_5to95'] < 5) 
            
    if len(imag) > 0:
        i = False
        if len(alert_isamp) > 0 and len(variables['ihist']) >= 100:
             i = (variables['KS_i'].pvalue < 0.05) and (np.median(alert_isamp) > variables['i_p95']) and (variables['chisq_i'] < 15) and (variables['chisq_i_5to95'] < 5) 
        
    if g and r and i:
        
        
        
        candidates.append(obj['objectId'])
        
        cand_coord = SkyCoord(ra=obj['objectData']['ramean'], dec=obj['objectData']['decmean'], unit=(u.degree, u.degree), frame='icrs')
        cand_search = Gaia.cone_search_async(cand_coord, radius=u.Quantity(1.0, u.arcsec))
        cand_match = cand_search.get_data()
        cand_plx = float(cand_match['parallax']) # in milliarcseconds
        cand_bp_rp= float(cand_match['bp_rp']) # BP - RP color in magnitudes
        cand_g_appmag = float(cand_match['phot_g_mean_mag']) # apparent g-bang magnitude
        cand_GMAG = cand_g_appmag - 5*np.log10(1000/cand_plx) + 5 # calculating absolute magnitude
        
        cand_color.append(cand_bp_rp)
        cand_gmag.append(cand_GMAG)
        
        
        if latest_CSV:
            if (latest_stats['Object'] == obj['objectId']).any():
                variables['num'] = float(latest_stats['num'][latest_stats['Object'] == 'ZTF24aandanb']) + 1
            else: variables['num'] = 1
        else:
            variables['num'] = 1
        
        
        statistics = {
            "Object": obj['objectId'],
            "g_X2": variables.get('chisq_g', np.nan),
            "r_X2": variables.get('chisq_r', np.nan),
            "i_X2": variables.get('chisq_i', np.nan),
            "g_5-95_X2": variables.get('chisq_g_5to95', np.nan),
            "r_5-95_X2": variables.get('chisq_r_5to95', np.nan),
            "i_5-95_X2": variables.get('chisq_i_5to95', np.nan),
            "g_KSpvalue": variables.get('g_pvalue', np.nan),
            "r_KSpvalue": variables.get('r_pvalue', np.nan),
            "i_KSpvalue": variables.get('i_pvalue', np.nan),
            "g_depth": variables.get('g_depth', np.nan), 
            "r_depth": variables.get('r_depth', np.nan),
            "i_depth": variables.get('i_depth', np.nan),
            "nhist": len(ztf.data['mag']),
            "g_nhist": len(ztf_gsamp),
            "r_nhist": len(ztf_rsamp),
            "i_nhist": len(ztf_isamp),
            "nalert": len(alertmjd),
            "g_nalert": len(gmag),
            "r_nalert": len(rmag),
            "i_nalert": len(imag),
            "ramean": obj['objectData']['ramean'],
            "decmean": obj['objectData']['decmean'],
            "gaiara": float(cand_match['ra']),
            "gaiadec": float(cand_match['dec']),
            "g_med": variables.get('g_med', np.nan),
            "r_med": variables.get('r_med', np.nan),
            "i_med": variables.get('i_med', np.nan),
            "disc_mjd": obj['objectData']['discMjd'],
            "latest_mjd": np.max(alertmjd),
            "gaia_sourceid": int(cand_match['SOURCE_ID']),
            "gaia_app_gmag": cand_g_appmag,
            "gaia_plx": cand_plx, 
            "gaia_abs_gmag": cand_GMAG,
            "gaia_BP-RP": cand_bp_rp,
            "num": variables['num']
            
            
        }

        with open(broker, 'a', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=headers)
            writer.writerow(statistics)            
             
      
        
        # Plotting
        
        figname = path+ '/lightcurves/'+str(obj['objectId'])+ '_' + mjdnow +'.pdf'
        
        plt.figure(figsize = (15, 5))
        plt.subplot(111)
        
        
        
        ga = plt.errorbar(gmjd, gmag, yerr = gerr, fmt = '.', color= 'mediumblue', label = 'g (alert)')
        ra = plt.errorbar(rmjd, rmag, yerr = rerr, fmt = '.', color= 'crimson', label = 'r (alert)')
        ia = plt.errorbar(imjd, imag, yerr = ierr, fmt = '.', color= 'mediumorchid', label = 'i (alert)')

        gh = plt.errorbar(ztf_gmjd, ztf_gsamp, yerr = ztf_gerr, fmt = '.', color= 'skyblue', alpha = 0.4, label = 'g (historic)')
        rh = plt.errorbar(ztf_rmjd, ztf_rsamp, yerr = ztf_rerr, fmt = '.', color= 'indianred', alpha = 0.4, label = 'r (historic)')
        ih = plt.errorbar(ztf_imjd, ztf_isamp, yerr = ztf_ierr, fmt = '.', color= 'plum', alpha = 0.4, label = 'i (historic)')
        
        ymin, ymax = plt.gca().get_ylim()
        xmin, xmax = plt.gca().get_xlim()
        
        mjdmin = plt.vlines(obj['objectData']['jdmin']-2400000.5, ymin, ymax, color='silver', linestyles=':', label='mjd min')  
        disc = plt.vlines(obj['objectData']['discMjd'], ymin, ymax, color='silver', linestyles=':', label='disc date')
        
        if 'ghist' in variables:
            gmed = plt.hlines(np.median(variables['ghist']), xmin, xmax, color='mediumblue', alpha= 0.4, linestyles=':', label='historic g-band median magnitude')
        
        if 'rhist' in variables:
            rmed = plt.hlines(np.median(variables['rhist']), xmin, xmax, color='crimson', alpha= 0.4, linestyles=':', label='historic g-band median magnitude')
        
        if 'ihist' in variables:
            imed = plt.hlines(np.median(variables['ihist']), xmin, xmax, color='mediumorchid', alpha= 0.4, linestyles=':', label='historic g-band median magnitude')
        
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.minorticks_on()
        plt.gca().invert_yaxis()
        plt.ylabel('Apparent Magnitude')
        plt.xlabel('MJD')
        plt.title(obj['objectId'])
        
        
        plt.savefig(figname, format='pdf', dpi=300)
        print ('Light curve saved to \"'+figname+'\".')
        
        plt.show()
        
        # Printing K-S test results
        if 'KS_g' in variables:
            print(obj['objectId']+' g-band KS result: '+ str(variables['KS_g']))
            

        if 'KS_r' in variables:  
            print(obj['objectId']+' r-band KS result: '+ str(variables['KS_r']))
            
            
        if 'KS_i' in variables: 
            print(obj['objectId']+' i-band KS result: '+ str(variables['KS_i']))
            
            
        print()

        # Printing Chi Square test results
        if 'chisq_g' in variables:
            print(obj['objectId']+' g-band Chi Square statistic: '+ str(variables['chisq_g']))
            print(obj['objectId']+' g-band 5th - 95th percentile Chi Square statistic: '+ str(variables['chisq_g_5to95']))

        if 'chisq_r' in variables:  
            print(obj['objectId']+' r-band Chi Square statistic: '+ str(variables['chisq_r']))
            print(obj['objectId']+' r-band 5th - 95th percentile Chi Square statistic: '+ str(variables['chisq_r_5to95']))

        if 'chisq_i' in variables: 
            print(obj['objectId']+' i-band Chi Square statistic: '+ str(variables['chisq_i']))
            print(obj['objectId']+' i-band 5th - 95th percentile Chi Square statistic: '+ str(variables['chisq_i_5to95']))

        print() 
        
        print(obj['objectId']+' has '+ str(len(gmag) + len(rmag) + len(imag))+ ' alert packet data points')
        print(obj['objectId']+' has '+ str(len(ztf.data['mag']))+ ' historic data points')
        
        print ()
        
    else:
        trash.append(obj['objectId'])
        trash_coord = SkyCoord(ra=obj['objectData']['ramean'], dec=obj['objectData']['decmean'], unit=(u.degree, u.degree), frame='icrs')
        trash_search = Gaia.cone_search_async(trash_coord, radius=u.Quantity(1.0, u.arcsec))
        trash_match = trash_search.get_data()
        trash_plx = float(trash_match['parallax']) # in milliarcseconds
        trash_bp_rp= float(trash_match['bp_rp']) # BP - RP color in magnitudes
        trash_g_appmag = float(trash_match['phot_g_mean_mag']) # apparent g-bang magnitude
        trash_GMAG = trash_g_appmag - 5*np.log10(1000/trash_plx) + 5 # calculating absolute magnitude
        
        trash_color.append(trash_bp_rp)
        trash_gmag.append(trash_GMAG)
        
    variables.clear()
   
print(str(len(trash)) + ' candidates rejected')
print(str(len(candidates)) + ' candidates remain') 

CMD = path + '/CMDs/CMD_' + mjdnow + '.pdf'

plt.figure(figsize=(10, 8))
plt.hist2d(color, gaia_GMAG, bins=[300,300], cmap='binary', alpha = 0.7, range=[[-1, 5.5], [-2, 20]])

plt.plot(trash_color, trash_gmag, 'd', color = 'white', markeredgecolor = 'indigo', zorder= 10, markersize= 7, label = 'Stars that passed Lasair filter')
plt.plot(cand_color, cand_gmag, '*', color = 'white', markeredgecolor = 'crimson', zorder= 11, markersize = 15, label = 'Stars that passed filter + statistical tests')
plt.legend(loc='upper right')
plt.colorbar(label='Density')
plt.xlabel('BP-RP Color (Magnitudes)')
plt.ylabel('Absolute G Magnitude')
plt.autoscale()
plt.gca().invert_yaxis()
plt.title('Gaia Stars Color-Magnitude Diagram')


plt.savefig(CMD, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0)
print ('Plot saved to \"'+CMD+'\".')
            
plt.show()


# Current criteria
# 
# Lasair query criteria:
# - objects.sgscore1 > 0.9 (closer to 1 means the object is more likely to be a star than a galaxy)
# - Object is not classified as a supernova, nuclear transient, nor active galactic nucleus by the sherlock classification scheme in Lasair
# - Object has at least 10 data points in its alert packet
# - Object has been crossmatched with the Gaia catalog
# - Object is within 0.5 arcseconds of the best source match from sherlock
# - Closest source match object is brighter than 16th magnitude in either g or r band
# - Alert magnitudes are all fainter than the reference magnitude
# - Alert occurred within the past 1 day
# - Object's first alert happened in 2024
# 
# Further criteria:
# - If object has alert data from a given filter, it must have at least 100 data points in its historic light curve from ZTF for that filter
# - Object passes two-sample K-S test with a p-value less than 0.05 in every filter
# - Object's historic data has a chi-square statistic less than 15
# - Object's historic data that falls within the 5th to 95th percentile range has chi-square statistic less than 5
# - Object's median alert magnitude is greater/fainter than the 95th percentile magnitude of the historic data in every filter
# 
# Note: K-S and chi square tests are calculated on historic data from over 100 days before the alert epoch and alert data from at least 30 days after the discovery date in Lasair

