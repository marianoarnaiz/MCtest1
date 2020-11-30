function earth2model(TTDATA0,dx,dy,dz,dray,AK135,Iterations)

    ## Transform real coordinates to Model Coordinates!
    LONG_0=minimum([TTDATA0[:,1]; TTDATA0[:,4]]);
    LAT_0=minimum([TTDATA0[:,2]; TTDATA0[:,5]]);
    #Larger window than data
    LONG_min=floor(minimum([TTDATA0[:,1]; TTDATA0[:,4]]))-1;
    LAT_min=ceil(minimum([TTDATA0[:,2]; TTDATA0[:,5]]))-1;
    LONG_max=floor(maximum([TTDATA0[:,1]; TTDATA0[:,4]]))+1;
    LAT_max=ceil(maximum([TTDATA0[:,2]; TTDATA0[:,5]]))+1;
    #Depth Min and Max
    DEPTH_min=floor(minimum(TTDATA0[:,3]));
    DEPTH_max=ceil(maximum(TTDATA0[:,6]));
    # Distances in the Model, from degrees to km
    LONG_Distance=round((LONG_max-LONG_min)*111.3195;sigdigits=2);
    LAT_Distance=round((LAT_max-LAT_min)*111.3195;sigdigits=2);
    # New coordinate system Base
    baselo=LONG_0-LONG_min; #%shift between (0,0) and data
    basela=LAT_0-LAT_min; #shift between (0,0) and data

    #data in the XYZ cartesian coordinate system
    TTDATA=[ ((TTDATA0[:,1].-LONG_0).*112).+baselo ((TTDATA0[:,2].-LAT_0).*112).+basela TTDATA0[:,3].-0.001 ((TTDATA0[:,4].+(pi/1000).-LONG_0).*112).+baselo ((TTDATA0[:,5].+(pi/1000).-LAT_0).*112).+basela TTDATA0[:,6].+0.001 TTDATA0[:,7] ];

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
                ModelSp[i,j,k]=1/iVp(vz2[k]); # Matrix Space of P wave Slowness value
            end
        end
    end

#ANOMALY=[+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0
#+1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0
#+1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0
#+1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0
#+1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0
#+1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0
#+1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0
#+1.0 +1.0 -0.1 -0.1 -0.1 -0.1 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +0.2 +0.2 +0.2 +0.2 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
#+1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0];


#for k=1:size(vz2,1)
#            ModelSp[:,:,k]=1 ./((1 ./ModelSp[:,:,k])+(1 ./ModelSp[:,:,k]).*ANOMALY); # Matrix Space of P wave Slowness value
#end

    #Create a Velocity Function
    #knots = ([x for x = xmin:dx:xmax], [y for y = ymin:dy:ymax], [z for z = zmin:1:zmax]);
    knots2 = ([x for x = xmin:dx:xmax], [y for y = ymin:dy:ymax], vz2);
    F_Sp = interpolate(knots2, ModelSp, Gridded(Linear()));

    #Evaluate the Slowness in the nodes
    Tomography_Sp=zeros(size(vx,1),size(vy,1),size(vz2,1),Iterations);
    Tomography_Sp[:,:,:,1]=F_Sp(vx,vy,vz2);

return TTDATA,xmin,xmax,ymin,ymax,zmin,movedz,vx,vy,vz2,F_Sp,Tomography_Sp,knots2
end
