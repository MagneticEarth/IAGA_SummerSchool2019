# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 12:27:07 2018

@author: djk
"""

import pandas as pd
import numpy as np

r2d = np.rad2deg
d2r = np.deg2rad

def IAGA2002_Header_Reader(IAGA2002_file):
    """
    This function counts the header and comment rows in an IAGA 2002 format
    file. It is designed to cope with the number of header lines being either
    12 or 13, and an arbitrary number of comment lines (including none).
    
    (The IAGA2002 format was last revised in June 2015 to allow an optional
    thirteenth header line 'Publication date'.   
    Ref: https://www.ngdc.noaa.gov/IAGA/vdat/IAGA2002/iaga2002format.html)
    
    The rows of data are preceded by a row of column headers starting with
    "DATE" in columns 0:3. This string cannot occur earlier in the file, so
    detecting the first occurence of this string may be used to count the total 
    number of header and comment lines.
          
    This function may be useful to define the number of rows to skip 
    (n_header + n_comment) in another function designed to read in the data. 
    While it is rather cumbersome, when reading in a long sequence of IAGA2002
    files, the 'safety first' approach would be to call this function for each
    file in case the number of header lines changes within the sequence of 
    files.
      
    Input parameter
    ---------------
    IAGA2002_file: string
    the full path and file name for the IAGA2002 data file
    
    Output
    ------
    A tuple: 
    with integer number of header rows (n_header),  integer number of comment 
    rows (n_comment), and headers, a dictionary containing the information in 
    the headers.
    
    Dependencies
    ------------
    pandas
      
    BGS Dependencies
    ----------------
    None
    
    Revision date
    -------------
    5 Feb 2018
      
    """
    COMMENT_STR   = '#'
    DATE_STR      = 'DATE'
    head          = '    '
    n_header = 0
    n_lines  = 0
    headers  = {}
    
    with open(IAGA2002_file) as ofile:
        while head[0:4] != DATE_STR:
            head = next(ofile)
            if head[1] != COMMENT_STR:
                key  = head[0:24].strip()
                val  = head[24:69].strip()
                headers[key] = val
                n_header    += 1
            n_lines += 1

    headers.pop(key)  # Remove the data column header line from the dictionary
    n_comment  = n_lines-n_header     # The number of comment lines
    n_header  -= 1                    # The number of header lines
    return (n_header, n_comment, headers)
    
###############################################################################  

def IAGA2002_Data_Reader(IAGA2002_file):
    """
    This function reads the data in an IAGA 2002 format file into a pandas
    dataframe.
      
    Input parameter
    ---------------
    IAGA2002_file: string
    the full path and file name for the IAGA2002 data file
    
    Output
    ------
    A pandas dataframe: 
    vals - has the data with a datetime index and the column labels from the
    IAGA2002 file
    
    Dependencies
    ------------
    pandas
      
    BGS Dependencies
    ----------------
    IAGA2002_Header_Reader
    
    Revision date
    -------------
    5 Feb 2018
      
    """
# Read the header and comment lines at the top of the file to get the number
# of rows to skip before reading the data
    header = IAGA2002_Header_Reader(IAGA2002_file)
    nskip  = header[0]+header[1]

# Read the data into a pandas dataframe (an IAGA2002 file has 'DATE' and 'TIME'
# as the first two column labels.) There's a trailing '|' on the column header
# line which is interpreted as the header for a column of nans and this 
# property is used to delete it. 

    DT_INDEX = 'DATE_TIME'
    vals = pd.read_csv(IAGA2002_file, 
                       delim_whitespace=True,
                       skiprows=nskip,
                       parse_dates=[DT_INDEX.split('_')],
                       index_col=DT_INDEX)
    vals.dropna(inplace=True, axis=1)
    
    return(vals)
    
###############################################################################  
"""
Created on Tue Jan 29 20:43:35 2019
    
    This function reads the hourly mean value data in yearly IAGA2002 format
    files into a pandas dataframe. (Note: The data may be reported in different
    ways in different years (e.g. DFHZ, FXYZ).)
      
    Input parameters
    ---------------
    obscode: the IAGA observatory code: string (3 or 4 characters)
    year_st: the start year for the data request
    year_fn: the final year for the data request
    folder : the location of the yearly hmv files
    
    Output
    ------
    A pandas dataframe: datareq
    This has columns of X, Y and Z data (only) and keeps the datetime index
    from the IAGA2002 files
    
    Dependencies
    ------------
    pandas
      
    Local Dependencies
    ----------------
    none
    
    Revision date
    -------------
    30 Jan 2019
    
Function to read in observatory annual mean files in IAGA2002 format

@author: djk
"""

def read_obs_hmv(obscode, year_st, year_fn, folder):
    OBSY   = obscode.upper()
    obsy   = obscode.lower()
# Read in the observatory data one year file at a time and construct filenames
    datareq = pd.DataFrame()        
    for year in range(year_st, year_fn+1):
        ystr    = str(year)
        file    = obsy + ystr + 'dhor.hor'
        fpf     =  folder + file
        tmp     = IAGA2002_Data_Reader(fpf)
        tmp.columns = [col.strip(OBSY) for col in tmp.columns]
        if('D' in tmp.columns):
            xvals, yvals  = dh2xy(tmp['D'], tmp['H'])
            tmp['X'] = xvals.round(decimals=1)
            tmp['Y'] = yvals.round(decimals=1)
        datareq = datareq.append(tmp[['X','Y', 'Z']])
    return(datareq)
    
###############################################################################  
"""
Created on Tue Jan 29 20:43:35 2019
    
    Function to read in the annual mean values for a single observatory
      
    Input parameters
    ---------------
    obscode:  the IAGA observatory code: string (3 or 4 characters)
    filename: the file of observatory annual mean values
    
    Output
    ------
    A pandas dataframe with the observatory annual mean values for all years
    available (indexed by year)
    
    Dependencies
    ------------
    pandas
         
    Revision date
    -------------
    30 Jan 2019

@author: djk

""" 
def read_obs_ann_mean(obscode, filename):
    count = 0
    with open(filename) as ofile:
        for line in ofile:
            if count==0:
                df = pd.DataFrame(columns=line.split())
                count += 1
            if str(obscode) in line:
                df.loc[count] = line.split()
                count +=1
    return(df)

###############################################################################  
"""
Created on Tue Jan 29 20:43:35 2019
    
    Function to read in the annual mean values for all observatories in a
    given year
      
    Input parameters
    ---------------
    obscode: the IAGA observatory code: string (3 or 4 characters)
    year:    the years are integers, (almost all) in the format yyyy.5
    
    Output
    ------
    A pandas dataframe with the observatory annual mean values for the year
    requested. 
    
    Dependencies
    ------------
    pandas
         
    Revision date
    -------------
    30 Jan 2019

@author: djk

"""     
def obs_ann_means_one_year(year, filename):
    
    count = 0
    with open(filename) as ofile:
        for line in ofile:
            if count==0:
                df = pd.DataFrame(columns=line.split())
                count += 1
            if str(year) in line:
                df.loc[count] = line.split()
                count +=1
    return(df)
    
###############################################################################
"""

Created on Tue Jan 29 20:43:35 2019
    
    Calculate X and Y from D (in 0.1 minutes) and H as reported in IAGA2002
    format files
      
    Input parameters
    ---------------
    D: declination (0.1 minutes)
    H: horizintal intensity (nT)
    
    Output
    ------
    A tuple: (X, Y), the North and East geomagnetic field components in nT
    
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    30 Jan 2019

@author: djk

""" 
   
def dh2xy(d, h):
    dec = d*np.pi/(180.*60.)
    return((h*np.cos(dec), h*np.sin(dec)))
    
###############################################################################
"""

Created on Tue Jan 29 20:43:35 2019
    
    Calculate D, H, I and F from (X, Y, Z)
      
    Input parameters
    ---------------
    X: north component (nT) 
    Y: east component (nT)
    Z: vertical component (nT)
    
    Output
    ------
    A tuple: (D, H, I, F)
    D: declination (degrees)
    H: horizontal intensity (nT)
    I: inclination (degrees)
    F: total intensity (nT)
    
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    28 April 2019
    
@author: djk

"""
def xyz2dhif(x, y, z):
    hsq = x*x + y*y
    hoz  = np.sqrt(hsq)
    eff = np.sqrt(hsq + z*z)
    dec = np.arctan2(y,x)
    inc = np.arctan2(z,hoz)
    return((r2d(dec), hoz, r2d(inc), eff))

###############################################################################
"""

Created on Tue Jan 29 20:43:35 2019
    
    Calculate secular variation in D, H, I and F from (X, Y, Z) and
    (Xdot, Ydot, Zdot)
      
    Input parameters
    ---------------
    X: north component (nT), and Xdot=dX/dt 
    Y: east component (nT), and Xdot=dX/dt 
    Z: vertical component (nT), and Xdot=dX/dt 
    
    Output
    ------
    A tuple: (Ddot, Hdot, Idot, Fdot)
    Ddot: rate of change of declination (degrees/year)
    Hdot: rate of change of horizontal intensity (nT/year)
    Idot: rate of change of inclination (degrees/year)
    Fdot: rate of change of total intensity (nT/year)
    
    Dependencies
    ------------
    numpy
         
    Revision date
    -------------
    28 April 2019
    
@author: djk

"""
def xyz2dhif_sv(x, y, z, xdot, ydot, zdot):
    h2  = x*x + y*y
    h   = np.sqrt(h2)
    f2  = h2 + z*z
    hdot = (x*xdot + y*ydot)/h
    fdot = (x*xdot + y*ydot + z*zdot)/np.sqrt(f2)
    ddot = r2d((xdot*y - ydot*x)/h2)*60
    idot = r2d((hdot*z - h*zdot)/f2)*60
    return((ddot, hdot, idot, fdot))