function initialrays(TTDATA)
    ALL_RAYS=zeros(dr,3,size(TTDATA,1));
    ALL_RAYS_Perturb=zeros(dr,3,size(TTDATA,1));
    ALL_T=zeros(size(TTDATA,1));
    ALL_T_Perturb=zeros(size(TTDATA,1));


for i=1:size(TTDATA,1)
    RECEIVER = [TTDATA[i,1], TTDATA[i,2],TTDATA[i,3]]; #Source Position
    SOURCE = [TTDATA[i,4],TTDATA[i,5],TTDATA[i,6]]; #Receiver Position
    #Initial RAY
    ALL_RAYS[:,:,i]=[collect(range(SOURCE[1],RECEIVER[1],length=dr))  collect(range(SOURCE[2],RECEIVER[2],length=dr)) collect(range(SOURCE[3],RECEIVER[3],length=dr))];
end
return ALL_T, ALL_RAYS, ALL_T_Perturb, ALL_RAYS_Perturb
end
