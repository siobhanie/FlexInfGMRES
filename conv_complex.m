clear all; clc; close all; 
randn('seed',0); rand('seed',0);

epsilon1=1.e-12; 
n=100; 
npts=200; 
m=150; %Krylov subspace dim


x=.5; y=.5; 
mu1=x + y*j; %run the algorithm on this (determines the stopping criteria) 

xx=linspace(-1.2,1.2,npts); yy=linspace(-1.2,1.2,npts);
[X,Y] = meshgrid(xx,yy); mus=X+Y*j; 
relres = zeros(size(mus));

randn('seed',0);
Av1=cell(m+1,1); In=speye(n,n); 
A0 = spdiags(rand(n,5),-2:2,n,n);
A1 = spdiags(rand(n,5),-2:2,n,n);
b=randn(n,1);

[Av1] = set_up2(A0,A1,In,m); 
b2=[b;zeros(m*(n),1)];
which=1; 

%Run the algorithm 1 time
[x1,Z,Q,H,e_inner1,tol_v1,rtilde_vec1]=flex_gmres(Av1,m,b2,n,mu1,which,epsilon1); 

figure(1)
for i=1:m
	xeval = x1(mu1,i); A_of_mu = A(A0,A1,In,mu1);
    res(i) = norm(A_of_mu*xeval(1:n)-b)/norm(b); 
end
semilogy(res)

for i=1:length(mus(1,:))
    for j=1:length(mus(:,1))

        %experiment for mus(i);
        mu = mus(i,j); 
        xeval = x1(mu,m); A_of_mu = A(A0,A1,In,mu);
        relres(i,j) = norm(A_of_mu*xeval(1:n)-b)/norm(b); 
    end
end


figure(2) 
pl = imagesc(xx,yy,(relres))
c=colorbar
caxis([1.e-8 1.e-1])
set(gca,'YDir','normal') 
xlabel('real')
ylabel('imag')
%title('Norm rel res at convergence, inner tol set for mu=.5+.5j')
%set(gca, 'Visible', 'off')
set(gca,'ColorScale','log')
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
%saveas(gca,'conv.pdf', 'pdf')







% data1 = [[0:1:m]',res'];
% data2 = [[0:1:m]',res2'];
% data3 = [[1:1:m]',e_inner'];
% data4 = [[1:1:m]',tol_v'];
% 
% header = {};
% header{1} = 'column 1'; header{2} = 'column 2'; 
% header = strjoin(header, ',');
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res1b.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res1b.csv',data1,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res2b.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res2b.csv',data2,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res3b.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res3b.csv',data3,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res4b.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/res4b.csv',data4,'-append');






function [x,Z,V,H,e_inner,tol_v,rtilde_vec]=flex_gmres(Av,m,b,n,mu,which,epsilon)
    [Z,V,H,e_inner,tol_v,rtilde_vec] = pc_arnoldi_fac(Av,n,m,b,mu,which,epsilon);
    %thist=1; 

    e = zeros(m+1,1); e(1) = 1;
    
    %GMRES
    x = @(mu,l) Z(:,1:l)*((-mu*H(1:l+1,1:l) + eye(l + 1,l))\(e(1:l+1)*norm(b))); 
    
    %FOM
    %x = @(mu,l) Z(:,1:l)*((-mu*H(1:l,1:l) + eye(l,l))\(e(1:l)*norm(b))); 
end

function [Z,V,H,e_inner,tol_v,rtilde_vec]=pc_arnoldi_fac(Av,n,m,r0,mu,which,epsilon)
%This one is NOT for timing purposes 

    V(:,1)=r0/norm(r0);
    
    for j=1:m
        
        if j==1
            rtilde_n=norm(r0); 
        else 
            rtilde_n = norm(rtilde_vec(j-1)); 
        end
        
        tol=1.e3; %Start with something very large 
        z=fastMatSolveK(Av,V(:,j),j,tol,which);
        P=V(1:length(z),j)-fastMatVecK(Av,z,j,n);
        
        while norm(P)*rtilde_n>epsilon
            tol=tol/10; 
            z=fastMatSolveK(Av,V(:,j),j,tol,which);
            P=V(1:length(z),j)-fastMatVecK(Av,z,j,n);
        end
        
        e_inner(j)=epsilon/rtilde_n;
        tol_v(j)=tol;        
        w=fastMatVecM(z,j,n);
        %w grows only by a nonzero block of dimension n
        
        %Double G-S
        Vactive=V(1:(j+1)*n,1:j);
        [h,w] = orthogonalize(Vactive,w);

        H(1:j,j)=h; H(j+1,j) = norm(w);
        V(1:length(w),j+1) = w/H(j+1,j);
        z=[z;zeros(n,1)]; Z(1:length(z),j)=z; %Increased memory

        %To compute rtilde 
        Im=eye(j + 1, j); 
        e = zeros(j+1,1); e(1) = 1;
        ym=((-mu*H(1:j+1,1:j) + eye(j + 1, j))\(e*norm(r0)));
        rtilde = V*(e*norm(r0)-(Im - mu*H)*ym);   
        rtilde_vec(j) = norm(rtilde); 
        
        %norm(A*Z-V(:,1:j+1)*H)
    end
end

function [h,y] = orthogonalize(Q,w)
    h=Q'*w;
    y = w - Q*h;
    g = Q'*y;
    y = y - Q*g;
    h = h + g;
end

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

function vec=fastMatVecM(v,j,n)
    vec=[zeros(n,1);v(1:end)];
end

function [vec]=fastMatSolveK(Av,z,j,tol,which)
    n=length(Av{1});

    if j==1
        
        if which==1
            vec= Av{1}\z(1:n) + randn(n,1)*tol;
        else
            vec=Av{1}\z(1:n);
        end
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
        
        if which==1
            vec1=Av{1}\vec + randn(n,1)*tol; 
        else
            vec1=Av{1}\vec;
        end
        vec=[-vec1;z(n+1:j*n)];
    end
end

function A_of_mu = A(A0,A1,In,mu)
    n=length(A0); 
    A_of_mu = -mu*In + A0 + A1*exp(-mu); 
end

