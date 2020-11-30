# It is a good practice specialize your function
function GetTime(ray0::Array{Float64, 2}, F_Sp)
    time = 0.0::Float64;
    #Loop to evaluate slowness
    for ri=2:size(ray0,1)
        #B = (F_Sp(ray0[ri-1,1],ray0[ri-1,2],ray0[ri-1,3]) + F_Sp(ray0[ri,1],ray0[ri,2],ray0[ri,3]))*0.5;
        #C = √((ray0[ri,1] - ray0[ri-1,1])^2 + (ray0[ri,2] - ray0[ri-1,2])^2 + (ray0[ri,3] - ray0[ri-1,3])^2);
        #time += B * C;
        time += ((F_Sp(ray0[ri-1,1],ray0[ri-1,2],ray0[ri-1,3]) + F_Sp(ray0[ri,1],ray0[ri,2],ray0[ri,3]))*0.5)*(√((ray0[ri,1] - ray0[ri-1,1])^2 + (ray0[ri,2] - ray0[ri-1,2])^2 + (ray0[ri,3] - ray0[ri-1,3])^2));
    end

    return time
end
