clear all; clc; clf; close all; 

m = 40;
epsilon = 1.e-12;
mu=2;
s=4;
mu2=mu/s;
 
[Av,b,A0,A1,A2,A3,A4]=two_D_inf_gmres_setup(m+1,s); 
n = length(b);
b2=[b;zeros(m*(n),1)];

%Set up for direct solve
[L,U,P,Q] = lu(A0); 
A_of_zero_inv = @(x) Q*(U\(L\(P*x)));

which=1;
[x,Z,Q,H,e_inner,tol_v,rtilde_vec,pvec,thist]=flex_gmres(A0,Av,m,b2,n,mu2,which,epsilon,A_of_zero_inv); 

which=2; 
[x2,Z2,Q2,H2,e_inner2,tol_v2,rtilde_vec2,pvec2,thist2]=flex_gmres(A0,Av,m,b2,n,mu2,which,epsilon,A_of_zero_inv); 

A_of_mu = (A0 + mu*A1 + (2*mu^2)*A2 + (mu^3)*A3 + sin(mu)*A4);
res(1) = norm(b)/norm(b);  %Initial 
res2(1) = norm(b)/norm(b); 
for i=1:m    
    xa=x(mu2,i); xb=x2(mu2,i); 
    res(i+1) = norm(A_of_mu*xa(1:n)-b)/norm(b); 
    res2(i+1) = norm(A_of_mu*xb(1:n)-b)/norm(b); 
end

figure(1)
semilogy(res,'*')
hold on
semilogy(res2,'ro')
semilogy(e_inner)
semilogy(tol_v,'*')
legend('res','LU','e_inner','tol of it meth')


function [x,Z,V,H,e_inner,tol_v,rtilde_vec,pvec,thist]=flex_gmres(A0,Av,m,b,n,mu,which,epsilon,A_of_zero_inv)
    
    [Z,V,H,e_inner,tol_v,rtilde_vec,pvec,thist] = pc_arnoldi_fac(A0,Av,n,m,b,mu,which,epsilon,A_of_zero_inv); 
    
    e = zeros(m+1,1); e(1) = 1;
    
    %GMRES
    x = @(mu,l) Z(:,1:l)*((-mu*H(1:l+1,1:l) + eye(l + 1,l))\(e(1:l+1)*norm(b))); 
    
    %FOM
    %x = @(mu,l) Z(:,1:l)*((-mu*H(1:l,1:l) + eye(l,l))\(e(1:l)*norm(b))); 
end

function [Z,V,H,e_inner,tol_v,rtilde_vec,pvec,thist]=pc_arnoldi_fac(A0,Av,n,m,r0,mu,which,epsilon,A_of_zero_inv)

    V(:,1)=r0/norm(r0);
    thist=zeros(m+1,1);
    tt=tic();
    e_inner=epsilon; 

    for j=1:m
        
        if which==1
            if j==1
                rtilde_n=norm(r0); 
            else 
                rtilde_n = norm(rtilde_vec(j-1)); 
            end
        end
         
        tol=e_inner(end);
        if which==1
            z=fastMatSolveK_alt_bicg(A0,Av,V(:,j),j,tol);
            P=V(1:length(z),j)-fastMatVecK(Av,z,j,n);
            while norm(P)*rtilde_n>epsilon
                tol=tol/10; 
                z=fastMatSolveK_alt_bicg(A0,Av,V(:,j),j,tol); 
                P=V(1:length(z),j)-fastMatVecK(Av,z,j,n);
            end
        elseif which==2
            z=fastMatSolveK2(A0,Av,V(:,j),j,tol,A_of_zero_inv); 
            P=0; 
        end

        pvec(j)=norm(P);
        
        if which == 1
            e_inner(j)=epsilon/rtilde_n;
            tol_v(j)=tol; 
        else 
            e_inner(j)=0; 
            tol_v(j)=tol; 
        end
        w=fastMatVecM(z,j,n);
        %w grows only by a nonzero block of dimension n
        
        %Double G-S
        Vactive=V(1:(j+1)*n,1:j);
        [h,w] = orthogonalize(Vactive,w);

        H(1:j,j)=h; H(j+1,j) = norm(w);
        V(1:length(w),j+1) = w/H(j+1,j);
        z=[z;zeros(n,1)]; Z(1:length(z),j)=z; %Increased memory

        %To compute rtilde 
        if which==1
            Im=eye(j + 1, j); 
            e = zeros(j+1,1); e(1) = 1;
            ym=((-mu*H(1:j+1,1:j) + eye(j + 1, j))\(e*norm(r0)));
            rtilde = V*(e*norm(r0)-(Im - mu*H)*ym);   
            rtilde_vec(j) = norm(rtilde); 
        else
            rtilde_vec(j) = 0; 
        end 

        thist(j+1)=toc(tt);
    end
end