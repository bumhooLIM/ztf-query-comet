# ztfssoquery/construct_fitsurl.py

import pandas as pd
from tqdm import tqdm
from pathlib import Path
import time
import re
import requests
from io import BytesIO
from astropy.io.votable import parse
from astroquery.jplhorizons import Horizons

def extract_lastrecnum(error_msg):
    """
    Extracts the last record number from an error message returned by Horizons.
    """
    lines = error_msg.split('\n')
    data_lines = [line.strip() for line in lines if re.match(r'^\d{8}', line.strip())]

    if data_lines:
        last_record_line = data_lines[-1]
        last_record_number = last_record_line.split()[0]
        return last_record_number
    return None

def query_sso_ephemeris(id, epochs, ephem_quantities):
    """
    Queries JPL Horizons for ephemerides of a given object.
    """
    obj = Horizons(id=id, location='I41', epochs=epochs) # ZTF observatory code: I41

    try:
        eph = obj.ephemerides(quantities=ephem_quantities)
        
    except ValueError as e:
        e = str(e)
        id_record = extract_lastrecnum(e)

        if id_record:
            obj = Horizons(id=id_record, location='I41', epochs=epochs)
            try:
                eph = obj.ephemerides(quantities=ephem_quantities)
            except Exception as e:
                print(f"Failed to get ephemerides for record id: {id_record}. Error: {e}")
        else:
            print(f"No valid record id found for: {id}.")
            return pd.DataFrame() # Return empty DF if failure
            
    return eph.to_pandas()

def query_ztf_metadata(ra, dec, width, height, jd_start, jd_end, intersect="overlaps"):
    """
    Queries the ZTF Science Archive to retrieve metadata of matching images.
    """
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
            # It is common to get no results, which can cause parsing errors or empty tables
            return None
    else:
        print(f"Error: API returned status code {response.status_code}")
        return None

def construct_fitsurl(ztf_row, is_cutout, cutout_size):
    """
    Constructs a FITS file URL for a given ZTF observation metadata row.
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
        
        added_url = f"?center={RA},{DEC}&size={cutout_size}&gzip=false"
        fits_url  += added_url
        
    return fits_url

def save_fitsurl(df_ztf, fpath_out, is_cutout, cutout_size):
    """
    Saves all FITS URLs from a ZTF DataFrame to a text file.
    """
    if df_ztf is None or df_ztf.empty:
        print("❌ No valid ZTF results found.")
        return

    urls = []
    for _, row in df_ztf.iterrows():
        try:
            url = construct_fitsurl(row, is_cutout=is_cutout, cutout_size=cutout_size)
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

def generate_fits_urls(
    primary_des, 
    start_date, 
    end_date, 
    interval_days=5, 
    rh_max=10, 
    vmag_max=21, 
    is_cutout=True, 
    cutout_size="10arcmin", 
    output_dir=None
):
    """
    Main execution logic to query Ephemerides, find matching ZTF images, 
    and save the list of FITS URLs.
    """
    
    # Determine output directory
    if output_dir is None:
        outdir_name = "".join(primary_des.split())
        OUTDIR = Path(f"./{outdir_name}") 
    else:
        OUTDIR = Path(output_dir)
    
    OUTDIR.mkdir(exist_ok=True, parents=True)

    # Query target ephemerides (with low resolution)
    eph = query_sso_ephemeris(
        id=primary_des,
        epochs={'start': start_date, 'stop': end_date, 'step': f"{interval_days}d"},
        ephem_quantities="1,3,9,19"
    )
    
    if eph.empty:
        print("No ephemerides found.")
        return

    # Masking ephemerides
    eph = eph[eph.r < rh_max].reset_index(drop=True)
    if 'Tmag' in eph.columns:
        eph = eph[eph.Tmag < vmag_max].reset_index(drop=True)
    
    eph.to_csv(OUTDIR / "eph.csv") 

    # Query ZTF metadata
    df_ztf = pd.DataFrame() 

    for _, eph_row in tqdm(eph.iterrows(), total=len(eph), desc="Querying ZTF"):
        try:
            datetime_jd = eph_row["datetime_jd"]
            ra = eph_row["RA"]
            dec = eph_row["DEC"]
            ra_dot = eph_row["RA_rate"]
            dec_dot = eph_row["DEC_rate"]

            # Query region based on target velocity (arcsec/hr > deg/day)
            ra_dot  = ra_dot  / 3600 * 24 
            dec_dot = dec_dot / 3600 * 24
            width   = abs(interval_days * ra_dot)  
            height  = abs(interval_days * dec_dot) 
            
            # ZTF query
            ztf = query_ztf_metadata(
                ra=ra,
                dec=dec,
                width=width,
                height=height,
                jd_start=datetime_jd - 0.5 * interval_days,
                jd_end=datetime_jd + 0.5 * interval_days
            )
            
            if ztf is not None and not ztf.empty:
                ztf = ztf.sort_values('obsjd').reset_index(drop=True)
                
                # Query ephemerides of given ZTF FOV
                eph_ztf = query_sso_ephemeris(
                    id=primary_des,
                    epochs=list(ztf.obsjd),
                    ephem_quantities="1,3,4,8,9,16,18,19,20,24,27"
                )
                ztf = pd.concat([ztf, eph_ztf], axis=1)
            
                # Check if target position is inside the FOV
                ztf["ra_min"]  = ztf[["ra1" , "ra2" , "ra3" , "ra4" ]].min(axis=1)
                ztf["ra_max"]  = ztf[["ra1" , "ra2" , "ra3" , "ra4" ]].max(axis=1)
                ztf["dec_min"] = ztf[["dec1", "dec2", "dec3", "dec4"]].min(axis=1)
                ztf["dec_max"] = ztf[["dec1", "dec2", "dec3", "dec4"]].max(axis=1)

                ztf_filtered = ztf[
                    (ztf.RA  > ztf.ra_min ) & 
                    (ztf.RA  < ztf.ra_max ) & 
                    (ztf.DEC > ztf.dec_min) & 
                    (ztf.DEC < ztf.dec_max)
                ]
                
                if not ztf_filtered.empty:
                    df_ztf = pd.concat([df_ztf, ztf_filtered], ignore_index=True)
        
        except Exception as e:
            print(f"Skipping ephemeris row; {e}") 

    df_ztf.to_csv(OUTDIR / "ztf.csv")

    # Save all FITS URLs to a text file
    save_fitsurl(df_ztf, fpath_out=OUTDIR / "fits_urls.txt", is_cutout=is_cutout, cutout_size=cutout_size)