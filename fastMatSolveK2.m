function [vec]=fastMatSolveK2(A0,Av,z,j,tol,A_of_zero_inv)
%Direct solve 
    n=length(Av{1});

    if j==1
        %vec= A0\z(1:n);
        vec=A_of_zero_inv(z(1:n)); 
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
        %vec1=A0\vec; 
        vec1=A_of_zero_inv(vec); 
        vec=[-vec1;z(n+1:j*n)];
    end
end
