function [Av,b,A0,A1,A2,A3,A4]=two_D_inf_gmres_setup(max_it,s)
    load('FEM.mat'); 
    n=length(b); 
    Av=cell(max_it,1);  
    Av{1} = A0;
    Av{2} = s*A1 + s*A4;
    coeff = 1; 
    for i=4:2:max_it 
       coeff = coeff*(-1);
       Av{i} = (1/factorial(i-1))*coeff*(s^(i-1))*A4; 
    end
    Av{4} = Av{4} + (1/factorial(4-1))*(6*(s^3))*A3;
    for i=3:2:max_it
        Av{i} = sparse(n,n); 
    end
    Av{3} = (1/factorial(3-1))*(4*(s^2))*A2;
    A_of_zero = A0;
    b=b'; 
end