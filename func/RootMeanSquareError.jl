function RootMeanSquareError(A,B)
    RMSE=(√).((sum((A-B).^2)/size(A,1)))
return RMSE

end
