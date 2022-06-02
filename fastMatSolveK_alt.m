function [vec]=fastMatSolveK_alt(A0,Av,z,j,tol,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB,switch1)
%With random perturbations 
    n=length(Av{1});
    j
    if j==1
        vec=agmg(A0,z(1:n),RESTART,TOL,MAXIT,VERBOSE,X0,IJOB); 
        %vec=gmres(A0,z(1:n),[],1.e-12,n);
        %vec=bicgstab(A0,z(1:n),tol,n);
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
        if switch1==0 %In the beginning you definitely need multigrid 
            vec1=agmg(A0,vec,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB);  
            %vec1=gmres(A0,vec,[],tol,n);
            %vec1=bicgstab(A0,vec,tol,n);
        else
            %vec1=agmg(A0,vec,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB);
            %vec1=gmres(A0,vec,[],tol,n);
            %vec1=bicgstab(A0,vec,tol,n);
            vec1=vec; 
            disp('nothing')
        end
        %vec1=A0\vec; 
        vec=[-vec1;z(n+1:j*n)];
    end
end