function RayBenderG(ray0::Array{Float64, 2}, dr::Int64, zmin::Float64, zmax::Int64, F_Sp, movedz::Float64)
       T0  = 0.0::Float64
       T1  = 0.0::Float64
       T2  = 0.0::Float64
       NCC = 0::Int64;
       g=[0.135335,  0.606531,  1.0,  0.606531,  0.135335]*movedz;
    T0=GetTime(ray0, F_Sp);
#NCC=0;
while NCC<dr-7
    NCC = 0::Int64;
    for p=2:dr-5

        #Bend a point of the ray up and down
        #Also get the travel time of each test

        rtest = ray0[p:p+4,3];
        ray0[p:p+4, 3] = rtest.+g;
        if ray0[p,3]<zmax
            T1=GetTime(ray0,F_Sp);
        else
            ray0[p:p+4, 3] = rtest
            T1=T0;
        end

        ray0[p:p+4, 3] = rtest.-g;
        if ray0[p,3]<zmax
            T2=GetTime(ray0,F_Sp);
        else
            ray0[p:p+4, 3] = rtest
            T2=T0;
        end



        if T1<T0 && T1<T2 # If RAY 1 is the solution
            ray0[p:p+4, 3] = rtest.+g;
            T0=T1;
            #println("T1")

        elseif T2<T0 && T2<T1 # If RAY 2 is the solution
            ray0[p:p+4, 3] = rtest.-g;
            T0=T2;
            #println("T2")

        else
            #NCC=NCC+1; # Keep Count of NO CHANGES
            NCC += 1;
        end

    end
end
return T0,ray0

end
