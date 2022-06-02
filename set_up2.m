function [Av1] = set_up2(A0,A1,In,m)
%This corresponds to the time-delay problem 

n=length(A0);
Av1=cell(m+1,1);
Av1{1} = A0+A1; %A(0)
Av1{2} = -In-A1; %A'(0)
j=1; 
for k=3:m+1 %Higher order derivatives
    Av1{k}=j*(1/factorial(k-1))*A1;
    j=j*(-1);
end


end