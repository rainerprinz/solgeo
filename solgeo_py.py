def solgeo(lat, lon, time, lsm):
    # Solar constant, W/m²
    Io = 1362  
    
    # Extract day of year and minute from datetime
    doy = time.dayofyear
    #mn = time.minute
    
    # Calculate day angle, Iqbal 1983 (1.2.2) day angle
    gamma = 2 * np.pi * (doy - 1) / 365   
    
    # Excentricity
    Eo = 1.000110 + 0.034221 * np.cos(gamma) + 0.001280 * np.sin(gamma) + \
          0.000719 * np.cos(2 * gamma) + 0.000077 * np.sin(2 * gamma)  
    
    # Declination
    delta = (0.006918 - 0.399912 * np.cos(gamma) + 0.070257 * np.sin(gamma) - \
              0.006758 * np.cos(2 * gamma) + 0.000907 * np.sin(2 * gamma) - \
              0.002697 * np.cos(3 * gamma) + 0.00148 * np.sin(3 * gamma))  
    
    # Equation of time
    Et = (0.000075 + 0.001868 * np.cos(gamma) - 0.032077 * np.sin(gamma) -
          0.014615 * np.cos(2 * gamma) - 0.04089 * np.sin(2 * gamma)) * 229.18  
    
    # Difference to local apparent time, Iqbal 1983 (1.4.2) difference to local apparent time
    LATd = Et + 4 * (lon - lsm)  
    
    # Local apparent time
    time_LT = time + pd.to_timedelta(LATd, unit='minutes')  
    
    # Decimal hours of local apparent time
    hro = time_LT.hour + time_LT.minute / 60  
    
    # Hour angle
    omega = (hro - 12) * np.deg2rad(15)  
    
    # Latitude in radians
    latrad = np.deg2rad(lat)  
    
    # Sun elevation
    sunel = np.arcsin(np.sin(latrad) * np.sin(delta) + np.cos(latrad) * np.cos(delta) * np.cos(omega))  
    
    # Top of atmosphere radiation
    TOA = Io * Eo * np.sin(sunel)  
    
    # Length of day
    omegas = np.arccos(np.tan(latrad) * -1 * np.tan(delta)) #nan for values outside [-1,1]
    LOD = 2 / 15 * np.rad2deg(omegas)  
