##  Markov Chain Monte Carlo Inversion for Travel Time TOMOGRAPHY
# Still Under Development: Test 1
println("****************************")
println("* Markov Chain Monte Carlo *")
println("*  Travel Time Tomography  *")
println("*         Test 1           *")
println("****************************")

clearconsole()
## Before using
# import Pkg; Pkg.add("Plots")
#using Pkg; Pkg.add("Interpolations")
#Pkg.add("VectorizedRoutines")
# ] DelimitedFiles
# import Pkg; Pkg.add("Statistics")

## USIN AND & INCLUDE
using DelimitedFiles, Interpolations, Plots, VectorizedRoutines, StatsPlots, Statistics

## FUNTIONS
include("func/earth2model.jl") #Creates the model in XYZ coordinates and the initial velocity models
include("func/GetTime.jl") #Computes Traveltimes in a slowness model
include("func/RayBender.jl") # Bends rays and computes traveltimes in 2.5 D
include("func/RayBenderG.jl") # Bends rays and computes traveltimes in 2.5 D
include("func/RayBender3D.jl") # Bends rays and computes traveltimes in 3D
include("func/initialrays.jl") #traces initial linear rays and creates other variables
include("func/RootMeanSquareError.jl") #compute RMSE error for 2 vectors
include("func/perturbations.jl") #Creates the perturbation matrixes for the inversion
include("func/likelyhood.jl") #Computes likelyhood between two vectors
#include("func/ECEF.jl")

## LOAD DATA SECTION
#Load AK135 Data
AK135=readdlm("data/jAK135.txt",Float64);
#Read input travel time data
TTDATA0=readdlm("data/TestData.txt",Float64);

println("Data Have Been Read!")

## INPUT PARAMETERS
const zmax = 300; #max depth in the model;
const dx = 100; #horizontal x cell size
const dy = 100; #horizontal y cell size
const dz = 20;  #vertical z cell size
const dray = 0.05; # Porcentage to change ray
const dr = 14; #Nodes in ray. THis inclures the source and receiver
const αmin=2.5 #min allowed velocity km/2
const αmax=9.5 #min allowed velocity km/2
Iterations=100; # Number of Iterations
σp=100; #Standard Deviation for the inversion (Sigma square)

## Setting things up

TTDATA,xmin,xmax,ymin,ymax,zmin,movedz,vx,vy,vz2,F_Sp,Tomography_Sp,knots2=earth2model(TTDATA0,dx,dy,dz,dray,AK135,Iterations);

println("Data In Model Coordinates")
println("Initial Model Created")

# Inicialize ray tracing  and initial time!
ALL_T, ALL_RAYS, ALL_T_Perturb, ALL_RAYS_Perturb =initialrays(TTDATA);
for i=1:size(ALL_T,1)
    ALL_T[i], ALL_RAYS[:,:,i] = RayBender3D(ALL_RAYS[:,:,i],   dr, zmin,zmax, F_Sp, movedz);
    #println("$i")
end
println("Initial Rays Traced and Initial Times Computed")

## Set Iterations
RMSE_Perturb=zeros(Iterations); # Initialize RMSE vector
L=zeros(Iterations); # Initialize likelyhood vector
Prob=zeros(Iterations); # Initialize Probability vector

RMSE_Perturb[1]=RootMeanSquareError(ALL_T[:,1],TTDATA[:,7]); # FIRST RMSE
L[1]=likelyhood(ALL_T[:,1],TTDATA[:,7],σp);
Prob[1]=1;

println("Initial RMSE, Likelyhood and Probability Computed")

## Set Perturbations
Perturbations_XYZ, Perturbations_α, MCrandom= perturbations(Iterations,vx,vy,vz2,αmin,αmax);

## INITIAL TOMOGRAPHY MODEL IS THE SAMPLED RAY VELOCItY MODEL


## ITERATIONS BEGIN AT 2
println("Begining Iterations!")
println(" ")

ACC=0; #accepted candidates counter
for IT=2:Iterations
    #println("VOY POR LA ITERACION: $IT")
    ALL_RAYS_BEFORE=ALL_RAYS; #Save the previous rays to perturbation in a variable
    #Perturbate Model
    Origincal_Cell_Value=Tomography_Sp[Perturbations_XYZ[IT,1],Perturbations_XYZ[IT,2],Perturbations_XYZ[IT,3],IT-1];
    Tomography_Sp[Perturbations_XYZ[IT,1],Perturbations_XYZ[IT,2],Perturbations_XYZ[IT,3],IT-1]=1/Perturbations_α[IT];
    #Candidate slowness model
    F_Tomo_Sp = interpolate(knots2, Tomography_Sp[:,:,:,IT-1], Gridded(Linear()));
    #Trace the rays and compute the times in the candidate model
    for i=1:size(TTDATA,1);
        # Trace ray by Pseudo-Bending Algorithm
        ALL_T_Perturb[i],   ALL_RAYS_Perturb[:,:,i]   = RayBender(ALL_RAYS[:,:,i],   dr, zmin,zmax, F_Tomo_Sp, movedz);
    end
    #Compute Parameters for the candidate model
        RMSE_Perturb[IT]=RootMeanSquareError(ALL_T_Perturb,TTDATA[:,7]);
        L[IT]=likelyhood(ALL_T_Perturb,TTDATA[:,7],σp);
        ProbCandidate=min(1,L[IT]/L[IT-1]);

        #if MCrandom[IT] <= ProbCandidate
        #    println("Accepted Perturbation in iteration: $IT")
        #elseif
        #    println("Rejected Perturbation in iteration: $IT")
        #end

     #IF error is reduced
    if RMSE_Perturb[IT] < RMSE_Perturb[IT-1]
        #Print this, keep count and save the rays
        println("Accepted Perturbation in iteration: $IT")
        #Save the Model in the Markov Chain
        Tomography_Sp[:,:,:,IT]=F_Tomo_Sp(vx,vy,vz2);
        #Save Probability
        Prob[IT]=ProbCandidate;
        #Save Rays
        global ALL_RAYS=ALL_RAYS_Perturb;
        #Keep count of accepted candidates
        global ACC+=1;

    else
        # Keep the same RMSE,
        RMSE_Perturb[IT] = RMSE_Perturb[IT-1];
        #Delete the perturbation from the model, and forget the traced rays
        Tomography_Sp[Perturbations_XYZ[IT,1],Perturbations_XYZ[IT,2],Perturbations_XYZ[IT,3],IT-1]=Origincal_Cell_Value;
        #Repeat the same model
        Tomography_Sp[:,:,:,IT]=Tomography_Sp[:,:,:,IT-1];
        Prob[IT]=Prob[IT-1];
        L[IT]=L[IT-1];
    end
end

## Print Report

pacc=ACC*100/Iterations;
rmse_reduction=((RMSE_Perturb[2]-RMSE_Perturb[Iterations])*100)/RMSE_Perturb[2];
meanresidual=mean(ALL_T_Perturb-TTDATA[:,7]);
stdresidual=std(ALL_T_Perturb-TTDATA[:,7]);

println(" ")
println("********************************** INVERSION REPORT ***************************")
println("*  Naive Walk Minimization Minimization Succesfully Run!        ")
println("*  Numbre of Iterations: $Iterations                            ")
println("*  Numbre of Accepted Candidates: $ACC                          ")
println("*  Percentage of admitance: $pacc%                              ")
println("*  RMSE reduction: $rmse_reduction%                             ")
println("*  Trave Time Residual -> Mean: $meanresidual Std: $stdresidual ")
println("*******************************************************************************")
println(" ")

## FIGURES! LET'S MAKE THEM PRETTY!

println("Plotting Figures! Hold your horses!")

#Figure 1. Velocity model and Stations and Events
p1=Plots.scatter(TTDATA[:,1],TTDATA[:,2],markershape = :dtriangle, color = [:orange],label="Stations")
p2=Plots.scatter(TTDATA[:,4],TTDATA[:,5],markershape = :circle, color = [:red],label="Events")
#Plots.scatter!(ALL_RAYS[:,1,:],ALL_RAYS[:,2,:],markershape = :circle, color = [:blue],label="")

p3=plot(2:Iterations,RMSE_Perturb[2:Iterations], color=:green, label="RMSE", xlabel = "Iterations", ylabel = "RMSE (s)")
p4=plot(2:Iterations,log10.(L[2:Iterations]), color=:red4, label="Log(L)", xlabel = "Iterations", ylabel = "Log10(Likelyhood)")
p5=plot(2:Iterations,Prob[2:Iterations], color=:blue4, label="RMSE", xlabel = "Iterations", ylabel = "Probability")


p6=histogram(ALL_T-TTDATA[:,7],color=:blue, fillalpha=0.3,label="Initial Residuals",xlabel = "Residual Time (s)",ylabel = "Frecuency (Counts)")
histogram!(ALL_T_Perturb-TTDATA[:,7],color=:red, fillalpha=0.3,label="Final Residuals",xlabel = "Residual Time (s)",ylabel = "Frecuency (Counts)")


plot(p1,p2,p3,p4,p5,p6)
savefig("figs/Figure1.pdf")

#= Figure 2
scatter(TTDATA[:,1],TTDATA[:,2],TTDATA[:,3], markershape = :dtriangle, color = [:orange],xlim=(xmin,xmax),ylim=(ymin,ymax),zlim=(zmin,zmax),label="Stations")
for ii=1:size(ALL_T,1)-1
    plot!(ALL_RAYS_Perturb[:,1,ii],ALL_RAYS_Perturb[:,2,ii],ALL_RAYS_Perturb[:,3,ii], color = [:blue],label="")
end
plot!(ALL_RAYS_Perturb[:,1,size(ALL_T,1)],ALL_RAYS_Perturb[:,2,size(ALL_T,1)],ALL_RAYS_Perturb[:,3,size(ALL_T,1)], color = [:blue],label="Rays")
scatter!(TTDATA[:,4],TTDATA[:,5],TTDATA[:,6] , markershape = :circle, color = [:red], label="Events")
savefig("figs/Rays.pdf")


gr()

p5= contourf(vy,vz2,1 ./ Tomography_Sp[5,:,:]',yflip=:true, c=:viridis)
p6= contourf(vy,vz2,1 ./ Tomography_Sp[10,:,:]',yflip=:true, c=:viridis)
p7= contourf(vy,vz2,1 ./ Tomography_Sp[12,:,:]',yflip=:true, c=:viridis)
p8= contourf(vy,vz2,1 ./ Tomography_Sp[18,:,:]',yflip=:true, c=:viridis)

plot(p5,p6,p7,p8)
savefig("figs/Profiles.pdf")
=#
