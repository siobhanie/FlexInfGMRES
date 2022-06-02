clear all; clc; close all;
randn('seed',0); rand('seed',0);

n=1000; 
mu=.2; m=70;
mu2=.1;
mu3=.05;
mu4=.025;

randn('seed',0);
Av1=cell(m+1,1); In=speye(n,n); 
A0 = spdiags(rand(n,5),-2:2,n,n);
A1 = spdiags(rand(n,5),-2:2,n,n);
b=randn(n,1);

[Av1] = set_up2(A0,A1,In,m); 

A_of_mu = A(A0,A1,In,mu);
A_of_mu2 = A(A0,A1,In,mu2);
A_of_mu3 = A(A0,A1,In,mu3);
A_of_mu4 = A(A0,A1,In,mu4);
b2=[b;zeros(m*(n),1)];

which=1; 
[x,Z,Q,H,e_inner,tol_v,rtilde_vec]=flex_gmres(Av1,m,b2,n,mu,which); 

for i=1:m
    xgmresi = x(mu,i);
    xgmresi2 = x(mu2,i);
    xgmresi3 = x(mu3,i);
    xgmresi4 = x(mu4,i);
    res(i) = norm(A_of_mu*xgmresi(1:n)-b)/norm(b);
    res2(i) = norm(A_of_mu2*xgmresi2(1:n)-b)/norm(b);
    res3(i) = norm(A_of_mu3*xgmresi3(1:n)-b)/norm(b);
    res4(i) = norm(A_of_mu4*xgmresi4(1:n)-b)/norm(b);
end

[x2,Z,Q,H,e_inner2,tol_v,rtilde_vec2]=flex_gmres(Av1,m,b2,n,mu2,which);
[x3,Z,Q,H,e_inner3,tol_v,rtilde_vec3]=flex_gmres(Av1,m,b2,n,mu3,which);
[x4,Z,Q,H,e_inner4,tol_v,rtilde_vec4]=flex_gmres(Av1,m,b2,n,mu4,which);

which=2;
[x1a,Z,Q,H,e_inner1,tol_v,rtilde_vec1]=flex_gmres(Av1,m,b2,n,mu,which);
[x2a,Z,Q,H,e_inner2,tol_v,rtilde_vec2]=flex_gmres(Av1,m,b2,n,mu2,which);
[x3a,Z,Q,H,e_inner3,tol_v,rtilde_vec3]=flex_gmres(Av1,m,b2,n,mu3,which);
[x4a,Z,Q,H,e_inner4,tol_v,rtilde_vec4]=flex_gmres(Av1,m,b2,n,mu4,which);

for i=1:m
    xgmresi = x1a(mu,i);
    xgmresi2 = x2a(mu2,i);
    xgmresi3 = x3a(mu3,i);
    xgmresi4 = x4a(mu4,i);
    resa(i) = norm(A_of_mu*xgmresi(1:n)-b)/norm(b);
    res2a(i) = norm(A_of_mu2*xgmresi2(1:n)-b)/norm(b);
    res3a(i) = norm(A_of_mu3*xgmresi3(1:n)-b)/norm(b);
    res4a(i) = norm(A_of_mu4*xgmresi4(1:n)-b)/norm(b);
end


data1 = [[1:1:m]',res'];
data2 = [[1:1:m]',res2'];
data3 = [[1:1:m]',res3'];
data4 = [[1:1:m]',res4'];
data5 = [[1:1:m]',resa'];
data6 = [[1:1:m]',res2a'];
data7 = [[1:1:m]',res3a'];
data8 = [[1:1:m]',res4a'];



% header = {};
% header{1} = 'column 1'; header{2} = 'column 2'; 
% header = strjoin(header, ',');
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res1.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res1.csv',data1,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res2.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res2.csv',data2,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res3.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res3.csv',data3,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res4.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res4.csv',data4,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res5.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res5.csv',data5,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res6.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res6.csv',data6,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res7.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res7.csv',data7,'-append');
% 
% fid = fopen('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res8.csv','w'); 
% fprintf(fid,'%s\n',header); 
% fclose(fid);
% dlmwrite('/Users/siobhanie/Dropbox/Flex Inf GMRES/doc/raw data/rp_res8.csv',data8,'-append');




clf;
figure(1)
semilogy([1:1:m],res,'*')
hold on
semilogy([1:1:m],res2,'*')
semilogy([1:1:m],res3,'*')
semilogy([1:1:m],res4,'*')
semilogy([1:1:m],resa,'o')
semilogy([1:1:m],res2a,'o')
semilogy([1:1:m],res3a,'o')
semilogy([1:1:m],res4a,'o')
legend('location','southwest')
legend('mu=.2 rp','mu=.1 rp','mu=.05 rp','mu=.025 rp','mu=.2 ds','mu=.1 ds','mu=.05 ds','mu=.025 ds')
xlabel('iters')
ylabel('rel res')
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf,'ran_pert_dir_solve.pdf', 'pdf')

figure(2)
semilogy([1:1:m],e_inner,'*')
hold on
semilogy([1:1:m],e_inner2,'*')
semilogy([1:1:m],e_inner3,'*')
semilogy([1:1:m],e_inner4,'*')
ylabel('norm E-inner')
xlabel('iters')
legend('location','northwest')
legend('mu=.2','mu=.1','mu=.05','mu=.025')
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf,'einner1.pdf', 'pdf')

figure(3)
semilogy([1:1:m],rtilde_vec,'*')
hold on
semilogy([1:1:m],rtilde_vec2,'*')
semilogy([1:1:m],rtilde_vec3,'*')
semilogy([1:1:m],rtilde_vec4,'*')
ylabel('norm r-tilde')
xlabel('iters')
legend('location','southwest')
legend('mu=.2','mu=.1','mu=.05','mu=.025')
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf,'rtilde1.pdf', 'pdf')


function [x,Z,V,H,e_inner,tol_v,rtilde_vec]=flex_gmres(Av,m,b,n,mu,which)
    [Z,V,H,e_inner,tol_v,rtilde_vec] = pc_arnoldi_fac(Av,n,m,b,mu,which);
    e = zeros(m+1,1); e(1) = 1;
    
    %GMRES
    x = @(mu,l) Z(:,1:l)*((-mu*H(1:l+1,1:l) + eye(l + 1,l))\(e(1:l+1)*norm(b))); 
    
    %FOM
    %x = @(mu,l) Z(:,1:l)*((-mu*H(1:l,1:l) + eye(l,l))\(e(1:l)*norm(b))); 
end

function [Z,V,H,e_inner,tol_v,rtilde_vec]=pc_arnoldi_fac(Av,n,m,r0,mu,which)
%This one is NOT for timing purposes 

    V(:,1)=r0/norm(r0);
    epsilon=1.e-10; %set to something very small 
    
    for j=1:m
        if j==1
            rtilde_n=norm(r0); 
        else 
            rtilde_n = norm(rtilde_vec(j-1)); 
        end
        
        tol=1.e3; %start with something very large
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
