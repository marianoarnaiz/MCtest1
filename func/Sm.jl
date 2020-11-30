#According to Mosegaard and Tarantola 1995 eq 15
function Sm(A,B)
    RMSE=(0.5)*(sum((A-B).^2))
return RMSE
end
