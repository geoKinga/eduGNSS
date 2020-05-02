#!/usr/bin/python
# -*- coding: UTF-8 -*- 

__copyright__   = "Copyright 2018, Faculty of Geodesy and Cartography, Warsaw University of technology"
__author__      = "Kinga Wezka"
__email__       = "kinga.wezka@pw.edu.pl"
__license__     = "GPL"
__version__     = "1.0.1"
__status__  = "under development"

"""
Written by Kinga Wezka, 2014 (2017 update to python 3.7)
Python  ver.: 2.7 and  3.7
Validated   : autopep8 -i --max-line-length 100 main.py
"""


import math
import numpy as np


class ConvertDate:
    """
    CLASS: The class 'ConvertDate' is use to calculate time
    DESCRIPTION:
        Available methods: \n       
        jd2mjd      : Convert Julian Day to Modified Julian Day \n
        mjd2jd      : Convert Modified Julian Day to Julian Day
        date2jd     : Convert a date to Julian Day. [day as fractional day!!!]
        jd2date     : Convert Julian Day to date
        days2hms    : Convert fractional days to hours, minutes, seconds
        hms2days    : Convert hours, minutes, seconds  to fractional days.
        jd2mjd      : Convert Julian Day to Modified Julian Day
        timeofday2sod- Convert time of the day to  Secound of the Day (SOD)
        jd2GPSdate  - Convert Julian Day to GPS date (GPS Week Number,Secound of the Day (SOD), Secound of the Week (SOW))
    REFERENCES:
        Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 4th ed., Duffet-Smith and Zwart, 2011.
    
    REMARKS:
        Note: The Python module assumes an infinitely valid Gregorian calendar.
        The Gregorian calendar took effect after 10-15-1582 
        and the dates 10-05 through  10-14-1582 never occurred. 
        Python 'ConvertDate' class will produce incorrect time if one date is from before 10-15-1582. 
    HTU:
        - import modul(this file) to another program:  
            from convertDate import *
        - define object of the class:
            date = ConvertDate()
        - invoke selected method:
            result = date.jd2date(2444244.5)  
    """
    def __init__(self): # metoda wykonuje się od razu po utworzeniu obiektu klasy
        self.sec_in_day     = 86400.0    # total number of the second in day
        self.half_week      = 302400.0   # number of the second in the half a week
        self.sec_in_week    = 604800.0   # total number of secound in the week
        self.JD80           = 2444244.5  # Julian Day : GPS  started at January 6, 1980 at midnight (0 h 0 m 0 s ) 
        self.MJD80          = 44244.0    # Modified Julian Day : GPS  started at January 6, 1980 at midnight (0 h 0 m 0 s ) 
   
    def mjd2jd(self, mjd):   # self - przyporzadkowanie do obiektu klasy
        '''
        Convert Modified Julian Day to Julian Day.
        -------
        Parameters:
        mjd : float  :Modified Julian Day
        -------
        Returns:
        jd : float   :Julian Day
        '''
        return mjd + 2400000.5
    
    def jd2mjd(self, jd):
        
        '''
        Convert Julian Day to Modified Julian Day
        -------
        Parameters:
        jd : float   :Julian Day
        -------
        Returns:
        mjd : float  :Modified Julian Day
        '''
        return jd - 2400000.5


    def jd2GPSdate(self, jd):
        '''
        Convert Julian Day to GPS date (GPS Week Number,Secound of the Day (SOD), Secound of the Week (SOW))
        -------
        Parameters:
        jd : float   :Julian Day
        -------
        Returns:
        mjd : float  :Modified Julian Day
        '''
        jd_diff  = self.jd2mjd(jd) - self.MJD80
        GPS_week = np.floor(jd_diff/7)
        DOW      =  int(np.floor(jd_diff - (GPS_week*7)))
        
        year, month, day, h, m, s = self.jd2datetime(jd)
        D0W = np.round(np.fmod(jd_diff, 7))
        
        SOD = self.timeofday2sod(hour=h, minute=m, second = s)
        if DOW == 0.0:
            SOW = SOD
        if DOW != 0.0:
            SOW = ((DOW) * self.sec_in_day) + SOD  
        return GPS_week, SOW, DOW, SOD
    
    
    def hhmmss2sod(self, hh, mm, ss):
        
        SOD = hh*3600 + mm*60 + ss
        return sod
    
    def timeofday2sod(self, hour=0, minute=0, second =0):
        
        sod = (hour*3600)+(minute*60) +second
        return sod

    def datetime2jd(self, year=2012, month=7, day=5, hour=0, minute=0, second=0):
        day_frac =self.hms2days(hour=hour, minute=minute , second= second)
        days = day + day_frac 
        jd = self.date2jd(year,month,days)
        return jd
    
    def jd2datetime(self,jd):
        year, month, day_frac = self.jd2date(jd)
        hd, m, s  = self.days2hms(day_frac)
        day = int(np.floor(day_frac))
        x = day_frac  - day
        h = int(np.floor(x*24))
        return year, month, day, h, m, s
    
    #
    def date2jd(self, year=2012, month=7, day=5):
        '''
        Convert a date to Julian Day.

        Parameters:
        year    : int    : Year as integer. Years preceding 1 A.D. should be 0 or negative.
                The year before 1 A.D. is 0, 10 B.C. is year -9.
        month   : int    : Month as integer, Jan = 1, Feb. = 2, etc.
        day     : float    : Day, may contain fractional part. eg. 12:00 > 0.5  (0.5x24)
        -------
        Returns
        jd : float    : Julian Day
        --------
        Examples
        Convert 6 a.m., February 17, 1985 to Julian Day
            >>> date_to_jd(1985,2,17.25)  
            >>> 2446113.75
        
        '''
        if month == 1 or month == 2:
            yearp = year - 1
            monthp = month + 12
        else:
            yearp = year
            monthp = month
        
        
        # --- this checks where we are in relation to October 15, 1582, the beginning of the Gregorian calendar.
        if ((year < 1582) or (year == 1582 and month < 10) or (year == 1582 and month == 10 and day < 15)):
            # --- before start of Gregorian calendar
            B = 0
        else:
            # --- after start of Gregorian calendar
            A = math.trunc(yearp / 100.)
            B = 2 - A + math.trunc(A / 4.)
            
        if yearp < 0:
            #C = math.trunc((365.25 * yearp) - 0.75)    #math.trunc = np.sign
            C = np.sign((365.25 * yearp) - 0.75)
        else:
            C = math.trunc(365.25 * yearp)
            
        D = math.trunc(30.6001 * (monthp + 1))
        jd = B + C + D + day + 1720994.5
        return jd


    def jd2date(self, jd):
        '''
        Convert Julian Day to date.
        
        -------
        Parameters:
            jd : float    :Julian Day
        -------
        Returns:
            year : int    : Year as integer. Years preceding 1 A.D. should be 0 or negative.
                The year before 1 A.D. is 0, 10 B.C. is year -9.
            month : int    : Month as integer, Jan = 1, Feb. = 2, etc.
            day : float    : Day, may contain fractional part. eg. 12:00 > 0.5  (0.5x24)
        --------
        Examples:
            Convert Julian Day 2446113.75 to year, month, and day.
        
        >>> jd_to_date(2446113.75) > (1985, 2, 17.25)
        
        '''
        jd = jd + 0.5
        F, I = math.modf(jd)  # The method modf() returns the fractional and integer parts of x in a two-item tuple. 
        I = int(I)
        
        A = math.trunc((I - 1867216.25)/36524.25) # return the truncated value of the input, element-wise.
        if I > 2299160:
            B = I + 1 + A - math.trunc(A / 4.)
        else:
            B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
        
        if G < 13.5:
            month = G - 1
        else:
            month = G - 13
            
        if month > 2.5:
            year = D - 4716
        else:
            year = D - 4715
            
        return year, month, day
   
  
    def hms2days(self, hour=0, minute=0 , second=0):
        '''
        Convert hours, minutes, seconds, and microseconds to fractional days.
        
        --------
        Parameters:
            hour : int       : Hour number.
            min : int        : Minute number.
            sec : int        : Second number.
        --------
        Returns:
            days : float    : A fractional number of days. Must be less than 1.
        -------          
        Examples
            >>> hmsm2days(hour=6)
            0.25
        
        '''
        days = second 
        days = minute + (days / 60.)
        days = hour + (days / 60.)
        return days / 24.


    def days2hms(self, days):
        '''
        Convert fractional days to hours, minutes, seconds
        ------- 
        Parameters:
            days : float    : A fractional number of days. Must be less than 1.
        -------        
        Returns:
            hour : int    : Hour integer number.
            min : int    : Minute integer number.
            sec : int    : Second integer number.
            ValueError  : if `days` is >= 1.
        --------    
        Example:
            >>> days_to_hmsm(0.1)
            >>> (2, 24, 0, 0)
        
        '''
        hours = days * 24.
        hours, hour = math.modf(hours)
        mins = hours * 60.
        mins, min = math.modf(mins)
        secs = mins * 60.
        return int(hour), int(min), np.round(secs)
    
    
    def UTCdate2GPS(self, year=2012, month=7, day=5, hour=0, minute=0, second=0):
        '''
        Convert year and day in year to GPS week and day in week
        ------- 
        Input:
            days : float    : A fractional number of days. Must be less than 1.
        -------        
        Output:
            GPS_week    : int    : GPS week
            day_in_week : int    : day in week (at midnight 00:00:00)

        --------    
        Example:
            >>> GPS_week, day_in_week =YERDOY2GPSWEEKDIW(year = 2012, doy = 187)
            >>> 1695,
        '''
        # --- for 1 January
        jd0 = self.datetime2jd(year=year, month=1, day=1, hour=0, minute=0, second=0)
        # --- for current date
        jd = self.datetime2jd(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
        DOY = jd - (jd0 - 1)
        GPS_week, SOW, DOW, SOD = self.jd2GPSdate(jd)
        return(int(GPS_week), int(SOW), int(DOW), int(SOD), int(DOY))
    
    def sow2ds(self, SOW):
        '''
        DIW, SOD = SOW2SODiDIW(SOW)
        '''
        DIW = int(SOW/86400)
        SOD = SOW - DIW*86400
        return(DIW, SOD)
        
    def sow2sod(self, sow):
        '''
        sow2sod(SOW) - converting sow (secound of the week) into 
        diw (day in the week) and sod (secound of the week)
        
        INPUT:
            sow - secound of the week
        OUTPUT:
            
        '''
        diw = int(sow/86400)  # day in week
        sod = sow - diw*86400 # secound in day
        return(diw, sod)

    def sod2hhmmss(self, sod):
        """
        sow2sod(SOW) - converting sow (secound of the week) into 
        diw (day in the week) and sod (secound of the week)
        INPUT:
            sow - secound of the week
        OUTPUT:
            
        """           
        hh  = math.floor(sod/3600)
        m   = sod - hh*3600.
        mm  = math.floor(m/60)
        ss  = m - mm*60
        # ----- time in format hhmmss
        hh_day  =  hh - (hh/24)
        time = "{:2.0f}:{:2.0f}:{:4.2f}".format(hh_day,mm,ss)
        return(time, hh, mm, ss)
    
    def one2Digits(self,one):
        """
        conver number of single digits into two digits:
        example: 2 > 02
        """
        d = str(int(one))
        if len(d) == 1:
            out = '0' + d
        else:
            out = d
        return str(out)
    
    def one3Digits(self,one):
        '''
        conver number of single digits into two digits:
        example: 2 > 002 or 23 > 023 
        '''
        d = str(one)
        if len(d) == 1.0:
            out = '00'+d
        elif len(d) == 2.0:
            out = '0'+d
        else:
            out = d
        return out
    
    def date_my_format(self, year, month, day, hour, minute, sec):
        '''
        converting date to my format as string: 2019-03-14|00:02:56.00
        '''
        y  = float(year)
        m  = float(month)
        d  = float(day)
        h  = float(hour)
        mi = float(minute)
        s  = float(sec)
#        print('{:4.0f}-{:02.0f}-{:02.0f}|{:02.0f}:{:2.0f}:{:05.2f}'.format(y, m, d, h, mi, s))
        return('{:4.0f}-{:02.0f}-{:02.0f}|{:02.0f}:{:2.0f}:{:05.2f}'.format(y, m, d, h, mi, s))
        
    def my_format2gps(my_format):
        '''
        converting date in my forma to GPS time: GPSweek, sow
        '''
        
        return gps_week, sow





if __name__=='__main__':
    convertdate = ConvertDate() # object of the class
    a = convertdate.date2jd(2019, 5, 1)
    print('aaa', a)
    b = convertdate.jd2datetime(a)
    print(b)

