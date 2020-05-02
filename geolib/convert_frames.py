#!/usr/bin/python
# -*- coding: UTF-8 -*-

__copyright__   = "Copyright 2018, Faculty of Geodesy and Cartography, Warsaw University of Technology"
__author__      = "Kinga Wezka"
__email__       = "kinga.wezka@pw.edu.pl"
__credits__      = "Kinga Wezka"
__copyright__   = "Copyright 2018 "
__license__     = "GPL"
__version__     = "1.0.1"
__status__      = "development"


from math import sin, cos, acos, asin, pi, atan, sqrt, atan2, radians, degrees
import numpy as np


#import sys
#assert ('linux' in sys.platform), "This code runs on Linux only." 


class ConvertFrames:
    """
    CLASS: The class 'ConvertFrames' is use to transform coorbinate between selected referance frames
    DESCRIPTION:
        Available methods: \n  
        deg2dms : converting angle given in decimal degree to deg, min, sec 
        dms2deg : convert angle given in deg, min, sec to decimal degree
        llh2xyz : convert geodetic frame (latitude, longitude, height) to geocentric reference frame x, y, z
        xyz2llh : lat, lon, alt = wgsxyz2lla(Xwgs,Ywgs,Zwgs)
        xyz2llaV:
        topocentrENU: topocentric coordinates East, North, Up = topocentrENU(Xwgs, Ywgs, Zwgs, ref_lat, ref_lon, ref_alt)
        
        AzElZeDist:
    
    REMARKS:
        - ...
    INFO:
        Written by Kinga Wezka, 2014 (2018 update to python 3.7)
        Python  ver.: 2.7 and  3.7
        Validated   : autopep8 -i --max-line-length 100 main.py
        ...
    HTU:
        - import modul(this file) to another program:  
            from convertDate import *
        - define object of the class:
            frame = ConvertFrames()
        - invoke selected method:
            result = date.llh2xyz(52, 21, 130)  
    """
    def __init__(self):
        self.OMEGAE84       = 7.292115*10**-5       # [rad/s]  - WGS84- Earth's  rotation rate as specified in ICD-GPS-200. # WGS84 [rads/s] - Earth's rotation rate [ICD-GPS-200C]   
        self.GM84           = 3.986005 * 10**14     # [m^3/s^2]- WGS84- universal  gravitation constant of gravitation [ICD-GPS-200C]
        self.GM             = 3.986004418 * 10E+14  # full number 
        self.G              = 6.67408 * 10**-11     # 1Newton [m^3 / kg s^2] - Gravitational cosntant
        self.f_WGS84        = 1/298.257223563;      #   WGS-84 Flattening.
        self.RE_WGS84       = 6378137.0             # [m] WGS84 equatorial radius 
        self.RP_WGS84       =  self.RE_WGS84 *(1-self.f_WGS84)   # [m] WGS84 polar radius
        self.REE            = 6378137.0             # [m] equatorial radius WGS84
        self.REM            = 6371000.0             # [m] mean radius
        self.REP            = 6356752.31413         # [m] polar radius
        self.flat           = 1/298.257223563       # WGS84 - Flattening
        self.REPclc         = self.REE*(1 - self.flat) # Polar radius (m)  = 6356752.31413: distance from perigee: q = a(1 − e) 
        self.ME             = 5.976 * 10**24        # [kg] mass of the earth    
        self.ECC            = 0.081819190842621    # WGS84 - Eccentricity 
    
    def deg2dms(self, decimalDeg):
        '''
        This method convert decimal degree to deg ,min, sec
        '''
        d = int(decimalDeg)
        m = int((decimalDeg - d) * 60)
        s = (decimalDeg - d - m/60.) * 3600.
        return(d, m, s)
    
    def dms2deg(self, dms):
        '''
        This method convert deg ,min, sec to decimal degree to 
        '''
        d = dms[0]
        m = dms[1]
        s = dms[2]
        decimal_degree = d + m/60 + s/3600
        return(decimal_degree)
    
    def llh2xyz(self, lat, lon, h):
        '''
        This method returns the position coordinate [X,Y,Z] given in the WGS-84 Earth Centered Earth Fixed 
        (ECEF) coordinate  for a user located at the goedetic coordinates lat,lon and h. The units 
        of the output position vector, p_e, are meters. latitude, longitude, altitude (reference 
        elipsoid WGS-84)
        INPUT:
            reference elipsoid WGS-84
            lat : latitude [degree] 
            lon : longitude [degree] 
            h   : altitude [meter]
            
        OUTPUT: 
            X,Y,Z  - geocentric coordinates
        '''
        #   Compute East-West Radius of curvature at current position
        R_E = self.RE_WGS84/(sqrt(1. - (self.ECC * sin(radians(lat)))**2))
        #  Compute ECEF coordinates
        X = (R_E + h)*cos(radians(lat)) * cos(radians(lon))
        Y = (R_E + h)*cos(radians(lat)) * sin(radians(lon))
        Z = ((1 - self.ECC**2)*R_E + h)*sin(radians(lat))
        return(X, Y, Z)

    def xyz2llaV(self, Xwgs, Ywgs, Zwgs):
        '''
        This method returns the position coordinate [lat, lon, h] given in the WGS-84
        INPUT:
            Xwgs,Ywgs,Zwgs  - geocentric coordinates
        OUTPUT: 
            reference elipsoid WGS-84
            lat - latitude [degree] 
            lon - longitude [degree] 
            alt - altitude/height [meter]
        '''
        #  reference elipsoid WGS-84: self.ECC, self.RE_WGS84, self.RP_WGS84
        # Calculate longitude
        lon = degrees(np.arctan2(Ywgs,Xwgs))

        #  Start computing intermediate variables needed to compute altitude/height
        normXY = np.array([[Xwgs,Ywgs]])
        p      = np.linalg.norm(normXY)
        E      = sqrt(self.RE_WGS84**2 - self.RP_WGS84**2)
        F      = 54*(self.RP_WGS84 * Zwgs)**2;
        G      = p**2 + (1 - self.ECC**2)*Zwgs**2 - (self.ECC * E)**2;
        c      = self.ECC**4 * F * p**2/G**3
        s      = (1 + c + sqrt(c**2 + 2*c))**(1/3)
        P      = (F/(3*G**2))/((s + (1/s) + 1)**2);
        Q      = sqrt(1 + 2 * self.ECC**4 * P)
        k_1    = -P * self.ECC**2 * p/(1 + Q)
        k_2    = 0.5 * self.RE_WGS84**2 *(1 + 1/Q)
        k_3    = -P * (1 - self.ECC**2) * Zwgs**2/(Q*(1 + Q))
        k_4    = -0.5 * P * p**2
        r_0    = k_1 + sqrt(k_2 + k_3 + k_4)
        k_5    = (p - self.ECC**2 * r_0)
        U      = sqrt(k_5**2 + Zwgs**2)
        V      = sqrt(k_5**2 + (1 - self.ECC**2) * Zwgs**2)
        alt = U*(1 - (self.RP_WGS84**2/(self.RE_WGS84 * V)))
        #  Compute additional values required for computing latitude
        z_0 = (self.RP_WGS84**2 * Zwgs)/(self.RE_WGS84*V)
        e_p = (self.RE_WGS84/self.RP_WGS84)*self.ECC
        lat = degrees(np.arctan((Zwgs + z_0 * (e_p)**2)/p) )
        return (lat, lon, alt)

    def xyz2llh(self, X, Y, Z, a = 6378137., finv = 298.2572221):
        '''
        Method to calculate geodetic coordinates latitude, longitude, height 
        INPUT: 
            - X,Y,Z     : Cartesian coordinates
            - a         : reference ellipsoid values semi-major axis
            - finv      : the inverse of flattening
        GRS 80: a= 6378137.0 m;    b= 6356752.314140 m;    finv=298.257222100882711
        WGS 84: a= 6378137.0 m;    b= 6356752.314245 m;    finv=298.257223563
        OUTPUT: 
            latitude1   [decimal degree]
            longitude   [decimal degree]
            hel         [meters]
            
        REMARKS:
            The units of linear parameters X,Y,Z,a must all agree (m,km,mi,ft,..etc) 
            The output units of angular quantities will be in decimal degrees (15.5 degrees not  15 deg 30 min).  
            The output units of h will be the  same as the units of X,Y,Z,a.
            latitude, longitude, height = XYZ2geod(self, X, Y, Z, a = 6378137., finv = 298.2572221):
        '''
        height  = 0
        tolsq   = 1.e-10
        r2d     = 180/np.pi #  radians-to-degree factor
        
        #  compute square of ECCentricity
        if finv < 1.e-20:
            esq = 0
        else:  
            esq = (2-1/finv)/finv
        oneesq = 1-esq
        P = sqrt(X**2 + Y**2) # first guess : P is distance from spin axix
        # direct calculation of longitude
        if P > 1.e-20:
            longitude = atan2(Y,X)*r2d
        else:
            longitude = 0 
        if longitude < 0:
            longitude = longitude + 360
        r = sqrt(P**2 + Z**2) # r is distance from origin (0,0,0)
        if r > 1.e-20:
            sinphi = Z/r
        else:
            sinphi = 0 
        latitude = asin(sinphi)   
        # initial value of height  =  distance from origin minus
        # approximate distance from origin to surface of ellipsoid
        if r < 1.e-20:
            height = 0
        height  = r - a * (1 - sinphi * sinphi/finv)
        height  = height
        # iterations
        maxit  = 10
        i=0
        while i < maxit:
            sinphi = sin(latitude);
            cosphi = cos(latitude);
            #  compute radius of curvature in prime vertical direction
            N_phi = a/sqrt(1-esq*sinphi*sinphi)
            #  compute residuals in P and Z
            dP = P - (N_phi + height) * cosphi;
            dZ = Z - (N_phi*oneesq + height) * sinphi;
            #  update height and latitude
            height  = height + (sinphi*dZ+cosphi*dP)
            latitude= latitude + (cosphi*dZ-sinphi*dP)/(N_phi + height)
            #  test for convergence
            if (dP*dP + dZ*dZ < tolsq):
                break
            if i == maxit: # if not converged--Warn user!!!
                print( ' Problem in toGEOD, did not converge')
            i+=1
        latitude = degrees(latitude)
        return(latitude, longitude, height)

    def topocentrENU(self, Xuser, Yuser, Zuser, lat_ref, lon_ref, h_ref):
        '''
        This function returns the position coordinates of a user at the WGS-84 ECEF coordiantes in east-north-up coordinates
        relative to the reference position located at lat_ref (latitude in degrees),     
        lon_ref (longitude in degrees) and h_ref (in meters).
        
        INPUT:
            Xuser, Xuser, Xuser : WGS-84 ECEF [meters]
            lat_ref             : [degrees] latitude,     
            lon_ref             : [degrees] longitude, 
            h_ref               : [meters] altitude
        OUTPUT: 
            lat,lon             : [degree], 
            U[ in [meters]
        EXAMPLE: 
            enu = xyz2enu(p_e,lat_ref,lon_ref,h_ref)           
        '''
        #  convert lat, lon, h  to XYZ WGS84 
        XYZ_ref = self.llh2xyz(lat_ref, lon_ref, h_ref)
        delta_X = Xuser - XYZ_ref[0]
        delta_Y = Yuser - XYZ_ref[1]
        delta_Z = Zuser - XYZ_ref[2] 
        #  Calculate ENU coordinates
        Up      =  cos(radians(lat_ref)) * cos(radians(lon_ref)) * delta_X + cos(radians(lat_ref)) * sin(radians(lon_ref)) * delta_Y + sin(radians(lat_ref)) * delta_Z
        East    = -sin(radians(lon_ref)) * delta_X + cos(radians(lon_ref)) * delta_Y
        North   = -sin(radians(lat_ref)) * cos(radians(lon_ref)) * delta_X - sin(radians(lat_ref)) * sin(radians(lon_ref)) * delta_Y + cos(radians(lat_ref)) * delta_Z
        return(East, North, Up)

    def AzElZeDist(self, Xrec, Yrec, Zrec, Xsat, Ysat, Zsat):
        '''
        Transformation of vector dx(sat-rec) into topocentric coordinate system with origin at X (receiver coordinate).
        This modul calculates the elevation and azimuth from a reference position specified in ECEF coordinates (e.g. receiver) 
        to another position specified in ECEF coordinates (e.g. satellite ):
        INPUT :
            Xrec, Yrec, Zrec - ECEF receiver coordinates  
            Xsat, Ysat, Zsat - ECEF satellite coordinates
        
        OUTPUT: 
            D    - vector length in units like the input
            Az   - azimuth from north positive clockwise [degrees]
            El   - elevation angle, [degrees]
        '''
        r2d = pi/180
        # Transformation of vector dx into topocentric coordinate system with origin at X.    
        phi, lam, h = self.xyz2llh(Xrec, Yrec, Zrec) # lat, lon, height 
        phi_rad = phi * r2d
        lam_rad = lam * r2d

        # define transformation matrix F: b-phi, l- lam
        #|-sl, -sb*cl, cb*cl |
        #| cl, -sb*sl, cb*sl |
        #| 0,      cb,  sb   |
        F  = np.matrix([[-sin(lam_rad), -sin(phi_rad)*cos(lam_rad), cos(phi_rad)*cos(lam_rad)],
                        [ cos(lam_rad), -sin(phi_rad)*sin(lam_rad), cos(phi_rad)*sin(lam_rad)],
                        [       0,               cos(phi_rad),                 sin(phi_rad)  ]])

        #  compute unit vector from observation station to satellite position
        dx    = ( np.matrix([Xsat-Xrec, Ysat-Yrec, Zsat-Zrec]) ).T  # było

        local_vector = F.T * dx
        E = local_vector[0]
        N = local_vector[1]
        U = local_vector[2]
     
        distHz = sqrt(E**2 + N**2)
        dist3d   =  sqrt((Xsat-Xrec)**2 + (Ysat-Yrec)**2 + (Zsat-Zrec)**2)
        if distHz < 1.e-20:
            print('distant is too short, probable the same points is used')
            Az = 0
            El = 90
            Ze = 0
        else:
            Az = (atan2(E,N)) / r2d 
            El = (atan2(U, distHz)) / r2d  #elevation angle
            Ze = (acos(U / dist3d)) / r2d #zenith angle
        if Az < 0:
            Az = Az + 360
        D = sqrt(dx[0,0]**2 + dx[1,0]**2 + dx[2,0]**2)
        return Az, El, Ze, D  # El [degree]; Az [degree]: D [units like the input]
    
    
    def AzZe(self, Xrec, Yrec, Zrec, Xsat, Ysat, Zsat):
        '''
        Transformation of vector dx(sat-rec) into topocentric coordinate system with origin at X (receiver coordinate).
        This modul calculates the elevation and azimuth from a reference position specified in ECEF coordinates 
        (e.g. receiver) to another position specified in ECEF coordinates (e.g. satellite ):
        INPUT :
            Xrec, Yrec, Zrec - ECEF receiver coordinates  
            Xsat, Ysat, Zsat - ECEF satellite coordinates
        OUTPUT: 
            D    - vector length in units like the input
            Az   - azimuth from north positive clockwise [degrees]
            El   - elevation angle, [degrees]
        '''
        dtr = pi/180
        # Transformation of vector dx into topocentric coordinate system with origin at X.    
        [phi, lam, h] = self.xyz2lla(Xrec, Yrec, Zrec) # lat, lon, alt 
        phi_rad = (phi)*dtr
        lam_rad = (lam)*dtr
        
        # define transformation matrix F: b-phi, l- lam
        #|-sl, -sb*cl, cb*cl |
        #| cl, -sb*sl, cb*sl |
        #| 0,      cb,  sb   |
        F  = np.matrix([[-sin(lam_rad), -sin(phi_rad)*cos(lam_rad), cos(phi_rad)*cos(lam_rad)],
                        [ cos(lam_rad), -sin(phi_rad)*sin(lam_rad), cos(phi_rad)*sin(lam_rad)],
                        [       0,               cos(phi_rad),                 sin(phi_rad)  ]])
        
        U  = np.matrix([[ cos(phi_rad)*cos(lam_rad)],
                        [ cos(phi_rad)*sin(lam_rad)],
                        [      sin(phi_rad)        ]])
    
        N = np.matrix([[-sin(phi_rad)*cos(lam_rad) ],
                       [-sin(phi_rad)*sin(lam_rad) ],
                       [      cos(phi_rad)         ]])    
        
        E  = np.matrix([[-sin(lam_rad)],
                        [ cos(lam_rad)],
                        [       0     ]])
    
#        verification: E == N X U: 

    
        #  compute unit vector from observation station to satellite position
        dx    = ( np.matrix([Xsat-Xrec, Ysat-Yrec, Zsat-Zrec]) ).T  # było
        #        dx    = ( np.matrix([Xrec-Xsat, Yrec- Ysat, Zrec- Zsat]) ).T 
        
        local_vector = F.T * dx
        E = local_vector[0]
        N = local_vector[1]
        U = local_vector[2]
        
        HZdist = sqrt(E**2 + N**2);
        dist   =  sqrt((Xsat-Xrec)**2 + (Ysat-Yrec)**2 + (Zsat-Zrec)**2)
        if HZdist < 1.e-20:
            print('distant is too short, probably the same points is used')
            Az = 0
            El = 0
            Z  = 0
        else:
            Az = (atan2(E,N))/dtr
            Ze = (acos(U/dist)) /dtr    #zenith angle
            El = (atan2(U,HZdist))/dtr  #elevation angle
            print('Az, Ze, El', Az, Ze, El)
            
        if Az < 0:
            Az = Az + 360
        D = sqrt(dx[0,0]**2 + dx[1,0]**2 + dx[2,0]**2);
        return El, Az, Ze, Ze+El    # El [degree]; Az [degree]: D [units like the input]

    def ll2gemagnetic(self, MJD, lon, lat):
        """
        Written by Kinga Wezka, June 2015
        DESCRIPTION:
            Calculating
        INPUT:
            MJD: Modified Julian Day
            lon: geographic coordinate: longitude
            lat: geographic coordinate: latitude
        OUTPUT:
            iono_2ho: geomagnetic longitude [deg]
        
        REFERENCES:
            [1] Gizawy, M. L. E. (2003) Development of an Ionosphere Monitoring Technique Using GPS 
            Measurements for High Latitude GPS Users The University of Calgary, The University of Calgary, 2003
            Appendix B: Geomagnetic coordinates calculations
       REMARKS:
           Hs to be validated
        """
        lon = radians(lon)
        lat = radians(lon)
        
        # gemagnetic pole coordinates
        Lon0_pole = radians( 78.8 + 4.283*10**(-2) * ( (MJD -46066.)/365.25) ) # gemagnetic pole longitide
        Lat0_pole = radians( 289.1 + 1.413*10**(-2) * ( (MJD -46066.)/365.25) )# gemagnetic pole latitude
        
        # dipolar coordinates
        Lambda_dip  = sin(lat) * sin(Lat0_pole)  + cos(lat) * cos(Lat0_pole) * cos(lon -Lon0_pole) 
        Phi_dip     =  (cos(lat) * sin(lon -Lon0_pole)) / cos(Lambda_dip)
        
        print( "geomagnetic", np.arcsin(Phi_dip), degrees(np.arcsin(Phi_dip)))
        return Lambda_dip, Phi_dip
    
    def groundtrack_latlon(self, xs, ys, zs):
        """
        function which converts ECEF position vectors to latitude and longitude
        Based on rvtolatlong.m in Richard Rieber's orbital library on mathwork.com
        
        Note: that this is only suitable for groundtrack visualization, not rigorous 
        calculations.
        """
        R = [xs, ys, zs]
        r_delta = np.linalg.norm(R[0:1]);
        sinA = R[1]/r_delta;
        cosA = R[0]/r_delta;
    
        Lon = atan2(sinA,cosA);
    
        if Lon < - pi:
            Lon = Lon + 2 * pi;
        Lat = asin(R[2]/np.linalg.norm(R));
        return degrees(Lat), degrees(Lon)
    
    def xyz2gemagnetic(self, X, Y, Z):
       '''
       INPUT:
            X, Y, Z : geodetic cartesian coordinate (ECEF)
            
            #lon = 7.035   # E > Geomagnetic 90.55E
            #lat = 48.116  # N > Geomagnetic 48.76N
            lon = 13.  # E > Geomagnetic 53.98E
            lat = 52.  # N > Geomagnetic 85.27N
       REMARKS:
           Hs to be finished
       '''      
       res = 2+5
       return res
   
    def satellite_body_fixed(time, XS):
        ''':
            time: GPS time
            XS
        
        INPUT
        
        OUTPUT:
            i - unit vector that completes the right-handed system
            j - resulting unit vector of the cross product of k vector with the unit vector from the satellite to Sun
            k - unit vector pointing from the Satellite Mass Centre (MC) to the Earth's centre
        
        #'''
        #t_sun = SP3.t_sun
        #X_sun = SP3.X_sun

        ##supposing t_sun regularly sampled
        ## Sun posirioning vector from sat_moon_pos
        #q = round((time - t_sun(1)) / p_rate) + 1
        #X_sun = X_sun(q)
        #e = (X_sun - X_sat) / norm(X_sun - X_sat) # es - sun
        #k = -X_sat/np.linalg.norm(X_sat) #ez
        ##j = np.cross(ez,es)
        #j = [k[2].*e[3]-k[3]*e[2]
             #k[3].*e[1]-k[1]*e[3]
             #k[1].*e[2]-k[2]*e[1]] # ey
        ##i = np.cross(j,k)
        #i = [j[2]*k[3]-j[3]*k[2]    #ex
             #j[3]*k[1]-j[1]*k[3]
             #j[1]*k[2]-j[2]*k[1]]

        #j = j / norm[j]
        #i = i / norm[i]
        return 1#i,j,k #Es = [ex,ey, ez] - transoformation matrix
        
 
if __name__=='__main__':
    frames = ConvertFrames() # class object
    X1 = 4075580.6852; Y1 = 931853.6596; Z1 = 4801568.0542 #IGS: WTZR 
    X2 = 4075535.3000; Y2 = 931822.18506; Z2 = 4801608.9150 #IGS: WTZ2 

    xw, yw, zw = frames.llh2xyz(52,21,130)
    xx, yx, zx = frames.llh2xyz(53,21,130)
 
    el, az, ze, d = frames.AzElZeDist(xw, yw, zw, xx, yx, zx)
    print( f'\n  az, el, ze, d: {az:.2f}, {el:.2f}, {ze:.2f}, {d:.2f}{az:.2f}, {el+ze:.2f}')
    
    
    e, n, u = frames.topocentrENU(xw, yw, zw, 53, 21, 130)
    print( f'\n e, n, u: {e:.2f}, {n:.2f}, {u:.2f}')



    

