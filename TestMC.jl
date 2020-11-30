clearconsole()
##  Markov Chain Monte Carlo Inversion for Travel Time TOMOGRAPHY
# Still Under Development: Test 1
println("****************************")
println("* Markov Chain Monte Carlo *")
println("*  Travel Time Tomography  *")
println("*         Test 1           *")
println("****************************")

## Before using
# import Pkg; Pkg.add("Plots")
#using Pkg; Pkg.add("Interpolations")
#Pkg.add("VectorizedRoutines")
# ] DelimitedFiles
# import Pkg; Pkg.add("Statistics")
# import Pkg; Pkg.add("StatsBase")


## USIN AND & INCLUDE
using DelimitedFiles, Interpolations, Plots, VectorizedRoutines, StatsPlots, Statistics

## FUNTIONS
include("func/earth2model.jl") #Creates the model in XYZ coordinates and the initial velocity models
include("func/earth2model2.jl") #Creates the model in XYZ coordinates and the initial velocity models
include("func/GetTime.jl") #Computes Traveltimes in a slowness model
include("func/RayBender.jl") # Bends rays and computes traveltimes in 2.5 D
include("func/RayBenderG.jl") # Bends rays and computes traveltimes in 2.5 D
include("func/RayBender3D.jl") # Bends rays and computes traveltimes in 3D
include("func/initialrays.jl") #traces initial linear rays and creates other variables
include("func/RootMeanSquareError.jl") #compute RMSE error for 2 vectors
include("func/perturbations.jl") #Creates the perturbation matrixes for the inversion
include("func/likelyhood.jl") #Computes likelyhood between two vectors
include("func/asknumber.jl") #Asks for the cut in the burn in
include("func/ECEF.jl") ##compute the geocentric rectangular coordinates
include("func/Sm.jl") ##compute the geocentric rectangular coordinates


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
Iterations=25000; # Number of Iterations
σp=50; #Standard Deviation for the inversion (Sigma square)

## Setting things up

TTDATA,xmin,xmax,ymin,ymax,zmin,movedz,vx,vy,vz2,F_Sp,Tomography_Sp,knots2=earth2model2(TTDATA0,dx,dy,dz,dray,AK135,Iterations);

println("Data in Model Coordinates")
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
        ALL_T_Perturb[i],   ALL_RAYS_Perturb[:,:,i]   = RayBender3D(ALL_RAYS[:,:,i],   dr, zmin,zmax, F_Tomo_Sp, movedz);
    end
    #Compute Parameters for the candidate model
        RMSE_Perturb[IT]=RootMeanSquareError(ALL_T_Perturb,TTDATA[:,7]);
        L[IT]=likelyhood(ALL_T_Perturb,TTDATA[:,7],σp);
        ProbCandidate=min(1,L[IT]/L[IT-1]);
        println("Iteration: $IT and P:$ProbCandidate ")
    #If the probability is greater than a random number we accept the candidates
    #this means usually that the perturbation improves the data fit or that we accept
    # a bad model randomly
    if MCrandom[IT] <= ProbCandidate #RMSE_Perturb[IT] < RMSE_Perturb[IT-1] #
        #Print this, keep count and save the rays
        rmseforshow=RMSE_Perturb[IT];
        showl=L[IT];
        println("Accepted Perturbation in iteration: $IT with RMSE: $rmseforshow and P:$ProbCandidate ")
        #Save the Model in the Markov Chain
        Tomography_Sp[:,:,:,IT]=F_Tomo_Sp(vx,vy,vz2);
        #Delete the perturbation from the model in the previous iteration
        Tomography_Sp[Perturbations_XYZ[IT,1],Perturbations_XYZ[IT,2],Perturbations_XYZ[IT,3],IT-1]=Origincal_Cell_Value;
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

## Time to select the Burn in

p01=plot(2:Iterations,RMSE_Perturb[2:Iterations], color=:green, xscale=:log10, label="RMSE", xlabel = "Iterations", ylabel = "RMSE (s)")
p02=plot(2:Iterations,RMSE_Perturb[2:Iterations], color=:green, label="RMSE", xlabel = "Iterations", ylabel = "RMSE (s)");
plot(p01,p02)
savefig("figs/RMSE.pdf")

println(" ")
println("Please! Check Figure RMSE.pdf in the figs directory...")
println("...I need this information to continue! :s ")
println(" ")
CUT=asknumber()
println(" ")
println("Thank you! :)")
println(" ")
println("Please wait a bit longer...")
println(" ")

## Some Models for Output and Figures
Mean_Model= mean(1 ./Tomography_Sp[:,:,:,CUT:Iterations],dims=4);
Median_Model= median(1 ./Tomography_Sp[:,:,:,CUT:Iterations],dims=4);
Std_Model= std(1 ./Tomography_Sp[:,:,:,CUT:Iterations],dims=4);
Var_Model= var(1 ./Tomography_Sp[:,:,:,CUT:Iterations],dims=4);

mean1D=reshape(mean(mean(Mean_Model,dims=2),dims=1),1,:)';
median1D=reshape(median(median(Median_Model,dims=2),dims=1),1,:)'
initial1D=reshape(mean(mean(1 ./F_Sp,dims=2),dims=1),1,:)'

## Rays in the Mean MODEL


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
println("*  Cut of Burn in: $CUT                            ")
println("*  RMSE reduction: $rmse_reduction%                             ")
println("*  Trave Time Residual -> Mean:~ $meanresidual Std:~ $stdresidual ")
println("*******************************************************************************")
println(" ")

## FIGURES! LET'S MAKE THEM PRETTY!

println("Plotting Figures! Hold your horses!")

#Stations and Events
p1=Plots.scatter(TTDATA0[:,1],TTDATA0[:,2],markershape = :dtriangle, color = [:orange],label="Stations",aspect_ratio=:equal)
p2=Plots.scatter(TTDATA0[:,4],TTDATA0[:,5],marker_z=TTDATA0[:,6], markershape = :circle,label="Events",aspect_ratio=:equal)
plot(p1, p2, layout = (2, 1))
savefig("figs/Stations&Events.pdf")

# Markov Chain
p3=plot(2:Iterations,RMSE_Perturb[2:Iterations], color=:green, label="RMSE", xlabel = "Iterations", ylabel = "RMSE (s)")
vline!([CUT],label="Burn In/Sampling")
p4=plot(2:Iterations,log10.(L[2:Iterations]), color=:red4, label="Log(L)", xlabel = "Iterations", ylabel = "Log10(Likelyhood)")
vline!([CUT],label="Burn In/Sampling")
p5=plot(2:Iterations,Prob[2:Iterations], color=:blue4, label="RMSE", xlabel = "Iterations", ylabel = "Probability")
vline!([CUT],label="Burn In/Sampling")
plot(p3, p4, p5, layout = (3, 1))
savefig("figs/Iterations.pdf")

#Residuals
p6=histogram(ALL_T-TTDATA[:,7],color=:blue, fillalpha=0.3,label="Initial Residuals",xlabel = "Residual Time (s)",ylabel = "Frecuency (Counts)")
histogram!(ALL_T_Perturb-TTDATA[:,7],color=:red, fillalpha=0.3,label="Final Residuals",xlabel = "Residual Time (s)",ylabel = "Frecuency (Counts)")
plot(p6)
savefig("figs/Residuals.pdf")

#1D Models
p7=plot(initial1D,vz2,color=:black,label="Initial Model", xlabel = "Velocity (km/s)", ylabel = "Depth (km)",yflip=:true)
plot!(mean1D,vz2,color=:red4,label="Mean Model", xlabel = "Velocity (km/s)", ylabel = "Depth (km)",yflip=:true)
plot!(median1D,vz2,color=:blue4,label="Median Model", xlabel = "Velocity (km/s)", ylabel = "Depth (km)",yflip=:true)
plot!(AK135[:,2],AK135[:,1],color=:green,label="AK135 Model", xlabel = "Velocity (km/s)", ylabel = "Depth (km)",yflip=:true)
ylims!((zmin,zmax))
plot(p7)
savefig("figs/1D_Models.pdf")


#Rays
scatter(TTDATA[:,1],TTDATA[:,2],TTDATA[:,3], markershape = :dtriangle, color = [:orange],xlim=(xmin,xmax),ylim=(ymin,ymax),zlim=(zmin,zmax),label="Stations")
for ii=1:size(ALL_T,1)-1
    plot!(ALL_RAYS_Perturb[:,1,ii],ALL_RAYS_Perturb[:,2,ii],ALL_RAYS_Perturb[:,3,ii], color = [:blue],label="",linealpha=0.3)
end
plot!(ALL_RAYS_Perturb[:,1,size(ALL_T,1)],ALL_RAYS_Perturb[:,2,size(ALL_T,1)],ALL_RAYS_Perturb[:,3,size(ALL_T,1)], color = [:blue],label="Rays",linealpha=0.3)
scatter!(TTDATA[:,4],TTDATA[:,5],TTDATA[:,6] , markershape = :circle, color = [:red], label="Events")
savefig("figs/Rays.pdf")


#Profiles
p8=contourf(vx,vz2,Mean_Model[:,4,:]',yflip=:true, c=:viridis)
p9=contourf(vx,vz2,Median_Model[:,4,:]',yflip=:true, c=:viridis)
p10=contourf(vx,vz2,Std_Model[:,4,:]',yflip=:true, c=:plasma)
p11=contourf(vx,vz2,Var_Model[:,4,:]',yflip=:true, c=:inferno)
p12=contourf(vy,vz2,Mean_Model[10,:,:]',yflip=:true, c=:viridis)
p13=contourf(vy,vz2,Median_Model[10,:,:]',yflip=:true, c=:viridis)
p14=contourf(vy,vz2,Std_Model[10,:,:]',yflip=:true, c=:plasma)
p15=contourf(vy,vz2,Var_Model[10,:,:]',yflip=:true, c=:inferno)
plot(p8,p9,p10,p11, layout = (2, 2))
savefig("figs/Profile_A.pdf")
plot(p12,p13,p14,p15, layout = (2, 2))
savefig("figs/Profile_B.pdf")

# Cell Histogram
CELLS=[rand(1:size(vx,1),8,1) rand(1:size(vy,1),8,1) rand(1:size(vz2,1),8,1)];
p16=histogram(1 ./Tomography_Sp[CELLS[1,1],CELLS[1,2],CELLS[1,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p17=histogram(1 ./Tomography_Sp[CELLS[2,1],CELLS[2,2],CELLS[2,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p18=histogram(1 ./Tomography_Sp[CELLS[3,1],CELLS[3,2],CELLS[3,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p19=histogram(1 ./Tomography_Sp[CELLS[4,1],CELLS[4,2],CELLS[4,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p20=histogram(1 ./Tomography_Sp[CELLS[5,1],CELLS[5,2],CELLS[5,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p21=histogram(1 ./Tomography_Sp[CELLS[6,1],CELLS[6,2],CELLS[6,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p22=histogram(1 ./Tomography_Sp[CELLS[7,1],CELLS[7,2],CELLS[7,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
p23=histogram(1 ./Tomography_Sp[CELLS[8,1],CELLS[8,2],CELLS[8,3],CUT:Iterations], bar_width=0.1, normalize=:probability, color=:blue, fillalpha=0.3,label="",xlabel = "Velocity (km/s)",ylabel = "P")
plot(p16,p17,p18,p19,p20,p21,p22,p23, layout = (4, 2))
savefig("figs/Cells_Histograms.pdf")
