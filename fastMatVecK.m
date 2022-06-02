function w2 = fastMatVecK(Av,z,j,n)
    w=zeros(n,1);
    for i=1:j
        w=w+Av{i}*z((i-1)*n+1:i*n);
    end
    if j==1
        w2=w; 
        return
    else
        w2=[w;z(n+1:end)];
        return
    end
end