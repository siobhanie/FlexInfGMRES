function vec=fastMatVecM(v,j,n)
    %Nonscaling version
    vec=[zeros(n,1);v(1:end)];
    
%     v = reshape(v,[n,j]); 
%     v2(:,1:j) = bsxfun (@rdivide, v, 1:j);
%     v2 = v2(:); 
%     vec=[zeros(n,1);v2(1:end)];
    
end