#According to Mosegaard and Tarantola 1995 eq 15
function likelyhood(A,B,σp)
    l1=sum((A-B).^2);
    l2=0.5*l1;
    L=exp(-l2/(σp^2))
    #l1=sum(abs(A-B)./σp);
    #L=exp(-l1)
    return L
end
