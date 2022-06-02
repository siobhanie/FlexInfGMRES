function [vec]=fastMatSolveK_alt_bicg(A0,Av,z,j,tol)
%With random perturbations 
    n=length(Av{1});
    j
    if j==1
        vec=bicgstab(A0,z(1:n),tol,n);
        return
    else
        vec=zeros(n,1);
        for i=1:j
            if i==1
                vec=vec-z(1:n);
            else
                vec=vec+Av{i}*z((i-1)*n+1:i*n);
            end
        end
        vec1=bicgstab(A0,vec,tol,n);
        vec=[-vec1;z(n+1:j*n)];
    end
end