#compute the geocentric rectangular coordinates
function ECEF(Lat,Long,Depth)
    ## Earth's Data
    R_Ecuador=6378.137; #Equatorial radius in km
    R_Polar=6356.752; #Polar radius in km
    N_LAT=(R_Ecuador*R_Ecuador)./((√).(R_Ecuador*R_Ecuador.*cosd.(Lat).*cosd.(Lat) + R_Polar*R_Polar.*sind.(Lat).*sind.(Lat)));

    # geocentric rectangular coordinates of the point as follow: https://en.wikipedia.org/wiki/Reference_ellipsoid
    X=(N_LAT-Depth).*cosd.(Lat).*cosd.(Long);
    Y=(N_LAT-Depth).*cosd.(Lat).*sind.(Long);
    #Z=(((R_Polar*R_Polar)/(R_Ecuador*R_Ecuador))*(N_LAT-Depth)).*sind.(Lat);
    Z=-(((R_Polar*R_Polar)/(R_Ecuador*R_Ecuador))*(Depth)).*sind.(Lat);
return X,Y,Z

end
