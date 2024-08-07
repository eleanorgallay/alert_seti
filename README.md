### Schema
- Object: Object ID from Lasair Broker.
   - Takes the form ZTFyyxxxxxxx.
- g_X2: Chi-Square statistic for historical g-band data.
   - This data is pulled from the ZTF archive and must be at least 100 days older than the object's alert epoch.
    We define the start of an alert epoch by the object's discovery date on Lasair.
    This date marks the first point in the ZTF alert packet, 30 days before the first data point that flagged an alert.
- r_X2: Chi-Square statistic for historical r-band data.
   - This data is pulled from the ZTF archive and must be at least 100 days older than the object's alert epoch.
    We define the start of an alert epoch by the object's discovery date on Lasair.
    This date marks the first point in the ZTF alert packet, 30 days before the first data point that flagged an alert.
- i_X2: Chi-Square statistic for historical i-band data.
   - This data is pulled from the ZTF archive and must be at least 100 days older than the object's alert epoch.
    We define the start of an alert epoch by the object's discovery date on Lasair.
    This date marks the first point in the ZTF alert packet, 30 days before the first data point that flagged an alert.
- g_5-95_X2
- r_5-95_X2
- i_5-95_X2
- g_KSpvalue
- r_KSpvalue
- i_KSpvalue 
- g_depth
- r_depth
- i_depth
- nhist
- g_nhist
- r_nhist
- i_nhist
- nalert
- g_nalert
- r_nalert
- i_nalert 
- ramean
- decmean
- gaiara
- gaiadec
- g_med
- r_med
- i_med
- disc_mjd
- latest_mjd
- gaia_sourceid 
- gaia_app_gmag
- gaia_plx
- gaia_abs_gmag
- gaia_BP-RP
- num
