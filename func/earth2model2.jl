function earth2model2(TTDATA0,dx,dy,dz,dray,AK135,Iterations)

## Earth's Data
    # Transform Earth to Model properly
    X_S,Y_S,Z_S=ECEF(TTDATA0[:,1],TTDATA0[:,2],TTDATA0[:,3]);
    X_E,Y_E,Z_E=ECEF(TTDATA0[:,4],TTDATA0[:,5],TTDATA0[:,6]);

    #Get some Values for the transformation
    LONG_0=minimum([X_S; X_E]);
    LAT_0=minimum([Y_S; Y_E]);
    #Larger window than data
    LONG_min=floor(minimum([X_S; X_E]))-dx;
    LAT_min=floor(minimum([Y_S; Y_E]))-dy;
    LONG_max=ceil(maximum([X_S; X_E]))+dx;
    LAT_max=ceil(maximum([Y_S; Y_E]))+dy;
    #Depth Min and Max
    DEPTH_min=floor(minimum(Z_S));
    DEPTH_max=ceil(maximum(Z_E));
    # Distances in the Model, from degrees to km
    LONG_Distance=round((LONG_max-LONG_min);sigdigits=2);
    LAT_Distance=round((LAT_max-LAT_min);sigdigits=2);

    baselo=LONG_0-LONG_min; #%shift between (0,0) and data
    basela=LAT_0-LAT_min; #shift between (0,0) and data

    #DATA IN MODEL COORDNATES
    TTDATA=[ X_S.-LONG_0.+baselo Y_S.-LAT_0.+basela Z_S.-0.001 X_E.-LONG_0.+baselo Y_E.-LAT_0.+basela Z_E.+0.001 TTDATA0[:,7] ];


    ## THE MODEL!
    # You can use const (minor improvements)
    xmin = 0; # min longitud in km
    xmax = LONG_Distance; # max longitud in km
    ymin = 0; # min latitude in km
    ymax = LAT_Distance; # max latitude in km
    zmin = DEPTH_min; #min depth in the model (we should include topography)
    movedz = dz*dray;
    # Create nodes for raytracing! This is a fine mesh in the Z direction!
    vx=xmin:dx:xmax; # Vector of X coordinates
    vy=ymin:dy:ymax; # Vector of Y coordinates
    vz2=collect([zmin:dz:0;0:dz:zmax]); #zmin:dz:zmax; # Vector of Z coordinates for MODEL

    iVp=LinearInterpolation(AK135[:,1],  AK135[:,2]);

    ModelSp=zeros(size(vx,1),size(vy,1),size(vz2,1));
    #Fill Slowness Model with AK135 Values
    for k=1:size(vz2,1)
        for j=1:size(vy,1)
            for i=1:size(vx,1)
                ModelSp[i,j,k]=1/8.0#iVp(vz2[k]); # Matrix Space of P wave Slowness value
            end
        end
    end

    #Create a Velocity Function
    #knots = ([x for x = xmin:dx:xmax], [y for y = ymin:dy:ymax], [z for z = zmin:1:zmax]);
    knots2 = ([x for x = xmin:dx:xmax], [y for y = ymin:dy:ymax], vz2);
    F_Sp = interpolate(knots2, ModelSp, Gridded(Linear()));

    #Evaluate the Slowness in the nodes
    Tomography_Sp=zeros(size(vx,1),size(vy,1),size(vz2,1),Iterations);
    Tomography_Sp[:,:,:,1]=F_Sp(vx,vy,vz2);

return TTDATA,xmin,xmax,ymin,ymax,zmin,movedz,vx,vy,vz2,F_Sp,Tomography_Sp,knots2
end
