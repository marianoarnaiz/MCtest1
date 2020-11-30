function perturbations(Iterations,vx,vy,vz2,αmin,αmax)
    #Location of the perturbations in the mode
    Perturbations_XYZ=[rand(1:size(vx,1),Iterations,1) rand(1:size(vy,1),Iterations,1) rand(1:size(vz2,1),Iterations,1)];
    #Initialize Velocity perturbation vector
    Perturbations_α=zeros(Iterations);
    # Fill the Velocity perturbation vector with geology
    for i=1:Iterations
        if vz2[Perturbations_XYZ[i,3]] < 10
            Perturbations_α[i]=rand(αmin:0.1:7);
        elseif 10 <= vz2[Perturbations_XYZ[i,3]] <= 70
            Perturbations_α[i]=rand(6:0.1:8.5);
        elseif 70 < vz2[Perturbations_XYZ[i,3]] <= zmax
            Perturbations_α[i]=rand(7.7:0.1:9.5);
        end
    end
    #Create the distribution of acceptance for the Markov Chain
    MCrandom=rand(0:0.001:1,Iterations);
return Perturbations_XYZ, Perturbations_α, MCrandom
end
