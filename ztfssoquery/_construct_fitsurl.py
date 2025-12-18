import pandas as pd
from tqdm import tqdm
from pathlib import Path
import time

#
# Parameters
#
PRIMARY_DES           = '2024 E1'                       # primary designation
STOBSDATE, ENDOBSDATE = '2025-01-01', '2025-10-30'      # start and end days
INTERVAL_DAYS         = 5                               # query time step in days
rh_max                = 10                              # maximum heliocentric distance to be quried in au
vmag_max              = 21                              # maximum v-magnitude to be quired
IS_CUTOUT             = True                            # retrieve cutout files
CUTOUT_SIZE           = "10arcmin"                      # cutout size

# Directory to save output file
outdir_name           = "".join(PRIMARY_DES.split())
OUTDIR                = Path(f"./{outdir_name}") 
OUTDIR.mkdir(exist_ok=True)

#
# Functions
#
def Extract_LastRecNum(error_msg):
    """
    Extracts the last record number from an error message returned by Horizons.

    Parameters:
        error_msg (str): The error message string from Horizons query.

    Returns:
        str or None: The last found record number or None if not found.
    """
    import re
    lines = error_msg.split('\n')
    data_lines = [line.strip() for line in lines if re.match(r'^\d{8}', line.strip())]

    if data_lines:
        last_record_line = data_lines[-1]
        last_record_number = last_record_line.split()[0]
        return last_record_number
    return None


def Query_Ephemerides(id, epochs, ephem_quantities):
    """
    Queries JPL Horizons for ephemerides of a given object.

    Parameters:
        id (str): Object identifier (e.g., '81P' for comet 81P/Wild 2).
        epochs (dict): Dictionary specifying observation time {'start': 'YYYY-MM-DD', 'stop': 'YYYY-MM-DD', 'step': 'Xd'}.
        ephem_quantities (str): Comma-separated list of ephemerides quantities.

    Returns:
        pd.DataFrame: Ephemerides data including datetime_jd, RA, and DEC.
    """
    from astroquery.jplhorizons import Horizons

    obj = Horizons(id=id, location='I41', epochs=epochs) # ZTF observatory code: I41

    try:
        eph = obj.ephemerides(quantities=ephem_quantities)
        
    except ValueError as e:
        e = str(e)
        id_record = Extract_LastRecNum(e)  # Extract record number for periodic comets

        if id_record:
            obj = Horizons(id=id_record, location='I41', epochs=epochs)
            try:
                eph = obj.ephemerides(quantities=ephem_quantities)
            except Exception as e:
                print(f"Failed to get ephemerides for record id: {id_record}. Error: {e}")
        else:
            print(f"No valid record id found for: {id}.")
            
    return eph.to_pandas()


def Query_ZTF(ra, dec, width, height, jd_start, jd_end, intersect="overlaps"):
    """
    Queries the ZTF Science Archive to retrieve metadata of matching images.

    Documentation: https://irsa.ipac.caltech.edu/docs/program_interface/ztf_api.html#accesscontrol

    Parameters:
        ra (float): Right Ascension in degrees.
        dec (float): Declination in degrees.
        jd_start (float): Start Julian Date.
        jd_end (float): End Julian Date.
        intersect (str): Intersection method (default: "overlaps").

    Returns:
        pd.DataFrame or None: DataFrame containing ZTF metadata or None if no results.
    """
    from astropy.io.votable import parse
    from io import BytesIO
    import requests

    ZTF_API_URL = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"

    params = {
        "POS": f"{ra},{dec}",
        "SIZE":f"{width},{height}",
        "INTERSECT": intersect,
        "WHERE": f"obsjd BETWEEN {jd_start} AND {jd_end}",
        "RESPONSEFORMAT": "VOTABLE"
    }

    response = requests.get(ZTF_API_URL, params=params)

    if response.status_code == 200:
        try:
            votable = parse(BytesIO(response.content))
            table = votable.get_first_table().to_table(use_names_over_ids=True)
            df = table.to_pandas()
            df["filefracday"] = df["filefracday"].astype('str')
            return df
        except Exception as e:
            print("Error parsing VOTable:", str(e))
            return None
    else:
        print(f"Error: API returned status code {response.status_code}")
        return None


def Construct_fitsurl(ztf_row, is_cutout, CUTOUT_SIZE):
    """
    Constructs a FITS file URL for a given ZTF observation metadata row.

    Documentation: https://irsa.ipac.caltech.edu/docs/program_interface/ztf_metadata.html
    
    <Filepath pattern>

    'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/'
    +year+'/'+month+day+'/'+fracday+'/ztf_'+filefracday+'_'+paddedfield+'_'+filtercode+'_c'+paddedccdid+'_'+imgtypecode+'_q'+qid+'_'+suffix

    <Example>
    https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/2018/0411/467847/ztf_20180411467847_000535_zr_c11_o_q3_sciimg.fits
    
    <Example cutout>
    https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/2018/0411/467847/ztf_20180411467847_000535_zr_c11_o_q3_sciimg.fits?center=255.57691,12.28378&size=50arcsec&gzip=false

    Parameters:
        ztf_row (pd.Series): A row of ZTF metadata.
        is_cutout (boolean): download cutout data
        CUTOUT_SIZE (str): cutout size (e.g., 1arcmin). Only available if is_cutout=True.  

    Returns:
        str: The URL for the FITS file.
    """
    ZTF_FITS_BASE_URL = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"

    FILEFRACDAY = str(ztf_row["filefracday"])
    year, monthday, fracday = FILEFRACDAY[:4], FILEFRACDAY[4:8], FILEFRACDAY[8:]

    FIELD       = f"{int(ztf_row['field']):06d}"
    FILTERCODE  = ztf_row["filtercode"]
    CCDID       = f"{int(ztf_row['ccdid']):02d}"
    IMGTYPECODE = ztf_row["imgtypecode"]
    QID         = ztf_row["qid"]

    fits_url = f"{ZTF_FITS_BASE_URL}/{year}/{monthday}/{fracday}/ztf_{FILEFRACDAY}_{FIELD}_{FILTERCODE}_c{CCDID}_{IMGTYPECODE}_q{QID}_sciimg.fits"
    
    if is_cutout:      
        RA  = ztf_row["RA"]
        DEC = ztf_row["DEC"]
        
        added_url = f"?center={RA},{DEC}&size={CUTOUT_SIZE}&gzip=false"
        fits_url  += added_url
        
    return fits_url


def Save_fitsurl(df_ztf, fpath_out, is_cutout, CUTOUT_SIZE):
    """
    Saves all FITS URLs from a ZTF DataFrame to a text file.

    Parameters:
        df_ztf (pd.DataFrame): DataFrame containing ZTF metadata.
        fpath_out (str or Path): Path to save the FITS URLs.

    Returns:
        None
    """
    if df_ztf is None or df_ztf.empty:
        print("❌ No valid ZTF results found.")
        return

    urls = []
    for _, row in df_ztf.iterrows():
        try:
            url = Construct_fitsurl(row, is_cutout=is_cutout, CUTOUT_SIZE=CUTOUT_SIZE)
            urls.append(url)
        except Exception as e:
            print(f"Error constructing URL for row {row}: {e}")

    if urls:
        with open(fpath_out, "w") as f:
            for url in urls:
                f.write(url + "\n")
        print(f"✅ Saved {len(urls)} FITS URLs to {fpath_out}")
    else:
        print("❌ No valid FITS URLs generated.")


# Query target ephemerides (with low resolution)
eph = Query_Ephemerides(
    id=PRIMARY_DES,
    epochs={'start':STOBSDATE, 'stop':ENDOBSDATE, 'step':f"{INTERVAL_DAYS}d"},
    ephem_quantities="1,3,9,19"
    )
# masking ephemerides
eph = eph[eph.r<rh_max].reset_index(drop=True)
eph = eph[eph.Tmag<vmag_max].reset_index(drop=True)
eph.to_csv(OUTDIR/"eph.csv") # save the epehemrides info.    

'''
quantities:
https://ssd.jpl.nasa.gov/horizons/manual.html#observer-table
        1. Astrometric RA & DEC
      * 2. Apparent RA & DEC
        3.   Rates; RA & DEC
      * 4. Apparent AZ & EL
        5.   Rates; AZ & EL
        6. Satellite X & Y, position angle
        7. Local apparent sidereal time
        8. Airmass and Visual Magnitude Extinction
        9. Visual magnitude & surface Brightness
       10. Illuminated fraction
       11. Defect of illumination
       12. Satellite angle of separation/visibility code
       13. Target angular diameter
       14. Observer sub-longitude & sub-latitude
       15. Sun sub-longitude & sub-latitude
       16. Sub-Sun position angle & distance from disc center
       17. North pole position angle & sistance from disc center
       18. Heliocentric ecliptic longitude & latitude
       19. Heliocentric range & range-rate
       20. Observer range & range-rate
       21. One-way down-leg light-time
       22. Speed of target with respect to Sun & observer
       23. Sun-Observer-Targ ELONGATION angle
       24. Sun-Target-Observer ~PHASE angle
       25. Target-Observer-Moon/Illumination%
       26. Observer-Primary-Target angle
       27. Position Angles; radius & -velocity
       28. Orbit plane angle
       29. Constellation Name
       30. Delta-T (TDB - UT)
     * 31. Observer-centered Earth ecliptic longitude & latitude
       32. North pole RA & DEC
       33. Galactic longitude and latitude
       34. Local apparent SOLAR time
       35. Earth->Site light-time
     > 36. RA & DEC uncertainty
     > 37. Plane-of-sky (POS) error ellipse
     > 38. Plane-of-sky (POS) uncertainty (RSS)
     > 39. Range & range-rate sigma
     > 40. Doppler/delay sigmas
       41. True anomaly angle
     * 42. Local apparent hour angle
       43. PHASE angle & bisector
       44. Apparent target-centered longitude of Sun (L_s)
     * 45. Inertial frame apparent RA & DEC
       46.   Rates: Inertial RA & DEC
     * 47. Sky motion: angular rate & angles
       48. Lunar sky brightness & target visual SNR
'''      

#
# Query ZTF metadata
#

df_ztf = pd.DataFrame() # query result

for _, eph_row in tqdm(eph.iterrows(), total=len(eph), desc="Querying ZTF"):

    try:
        datetime_jd, ra, dec, ra_dot, dec_dot = eph_row[["datetime_jd", "RA", "DEC", "RA_rate", "DEC_rate"]] # RA_rate, DEC_rate in arcsec/hour
        
        # Query region based on target velocity
        ra_dot  = ra_dot  / 3600 * 24 # arcsec/hr > deg/day
        dec_dot = dec_dot / 3600 * 24
        width   = abs(INTERVAL_DAYS * ra_dot)  # width of query size  (deg)
        height  = abs(INTERVAL_DAYS * dec_dot) # height of query size (deg)
        
        # ZTF query
        
        ztf = Query_ZTF(
            ra=ra,
            dec=dec,
            width=width,
            height=height,
            jd_start=datetime_jd - 0.5*INTERVAL_DAYS,
            jd_end=datetime_jd + 0.5*INTERVAL_DAYS
            )
        # print("Query done.")
        ztf = ztf.sort_values('obsjd').reset_index(drop=True) # make sure that query dataframe is ordered in time.
        
        if not ztf.empty:
            
            # Query ephemerides of given ZTF FOV
            eph_ztf = Query_Ephemerides(
                id=PRIMARY_DES,
                epochs=list(ztf.obsjd),
                ephem_quantities="1,3,4,8,9,16,18,19,20,24,27"
                )
            ztf = pd.concat([ztf, eph_ztf], axis=1)
        
            # check if target position is inside the FOV
            ztf["ra_min"]  = ztf[["ra1" , "ra2" , "ra3" , "ra4" ]].min(axis=1)
            ztf["ra_max"]  = ztf[["ra1" , "ra2" , "ra3" , "ra4" ]].max(axis=1)
            ztf["dec_min"] = ztf[["dec1", "dec2", "dec3", "dec4"]].min(axis=1)
            ztf["dec_max"] = ztf[["dec1", "dec2", "dec3", "dec4"]].max(axis=1)

            ztf_filtered      = ztf[
                (ztf.RA  > ztf.ra_min ) & 
                (ztf.RA  < ztf.ra_max ) & 
                (ztf.DEC > ztf.dec_min) & 
                (ztf.DEC < ztf.dec_max)
                ]
            
            # Attach to result dataframe
            if not ztf_filtered.empty:
                
                df_ztf = pd.concat([df_ztf, ztf_filtered], ignore_index=True)

        # time.sleep(3)
    
    except Exception as e:
        print(f"Skipping ephemeris row; {e}") 

df_ztf.to_csv(OUTDIR/"ztf.csv") # save the ZTF metadata

# Save all FITS URLs to a text file
Save_fitsurl(df_ztf, fpath_out=OUTDIR/"fits_urls.txt", is_cutout=IS_CUTOUT, CUTOUT_SIZE=CUTOUT_SIZE)
# the resultant fits url can be downloaded based on the following command line
# xargs -n 1 curl -O < fits_urls.txt
