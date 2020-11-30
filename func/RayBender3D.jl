function RayBender3D(ray0::Array{Float64, 2}, dr::Int64, zmin::Float64, zmax::Int64, F_Sp, movedz::Float64)
       T0  = 0.0::Float64
       T1  = 0.0::Float64
       T2  = 0.0::Float64
       NCC = 0::Int64;
    T0=GetTime(ray0, F_Sp);
#NCC=0;
while NCC<dr-2
    NCC = 0::Int64;
    for p=2:dr-1

        #Bend a point of the ray up and down
        #Also get the travel time of each test

        #ray down

        rtestz = ray0[p,3];
        ray0[p, 3] = rtestz + movedz;
        if ray0[p,3]<zmax
            T1=GetTime(ray0,F_Sp);
        else
            ray0[p, 3] = rtestz;
            T1=T0;
        end

        #ray up

        ray0[p, 3] = rtestz - movedz;
        if ray0[p,3]>zmin
            T2=GetTime(ray0,F_Sp);
        else
            ray0[p, 3] = rtestz;
            T2=T0;
        end

        #ray right
        rtesty = ray0[p,2];
        ray0[p, 2] = rtesty + movedz;
        if ray0[p,2]<ymax
        T3=GetTime(ray0,F_Sp);
        else
            ray0[p, 2] = rtesty;
            T3=T0;
        end

        #ray left
        rtesty = ray0[p,2];
        ray0[p, 2] = rtesty - movedz;
        if ray0[p,2]>ymin
        T4=GetTime(ray0,F_Sp);
        else
            ray0[p, 2] = rtesty;
            T4=T0;
        end


        if T1<T0 && T1<T2 && T1<T3 && T1<T4 # If RAY 1 is the solution
            ray0[p, 3] = rtestz + movedz;
            T0=T1;
            #println("T1")

        elseif T2<T0 && T2<T1 && T2<T3 && T2<T4  # If RAY 2 is the solution
            ray0[p, 3] = rtestz - movedz;
            T0=T2;
            #println("T2")

        elseif T3<T0 && T3<T1 && T3<T2 && T3<T4  # If RAY 2 is the solution
            ray0[p, 2] = rtesty + movedz;
            T0=T3;
            #println("T2")
        elseif T4<T0 && T4<T1 && T4<T2 && T4<T3  # If RAY 2 is the solution
            ray0[p, 2] = rtesty - movedz;
            T0=T4;
        else
            #NCC=NCC+1; # Keep Count of NO CHANGES
            NCC += 1;
        end

    end
end
return T0,ray0

end
