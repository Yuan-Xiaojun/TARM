function [Xhat,mse] = TARM(b,dim,A,At,r,params,errorfunction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TARM algorithm for Robust Affine Rank Minimization problem
%%% b: observation of low-rank matrix, of size m*1
%%% A: linear maping, n1*n2-->m
%%% At: adjunt linear maping of A, n-->n1*n2
%%% r: initial rank of X
%%% param.mu: step size type, 0,1,2..
%%% errorfunction: calculate psnr of each iteration
%%% lansvd package required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=dim.m;
n1=dim.n1;
n2=dim.n2;
n=n1*n2;
delta=m/n;
% params
mu=params.mu;
iter=params.iter; % max iteration time
tol=params.tol; % tol for stopping
divtype=params.divtype; %divergence type
epsi=0.0001; % divergence simulation parameter
% step size selection
times=1;
Xhat = zeros(n1,n2);
[U,Sig,V] = lansvd(At(b-A(Xhat)),r);
Xhat = U*Sig*V';
for ii=1:iter
    Grad=At(b-A(Xhat));
    % adoptive step size
    if mu == 0
        step=1/delta;
    end
    if mu == 1
        Pu=U*U';
        PsG=Pu*Grad;
        step=norm(PsG,'fro')^2/norm(A(PsG))^2;
    end
    if mu == 2
        Pv=V*V';
        PsG = Grad*Pv;
        step=norm(PsG,'fro')^2/norm(A(PsG))^2;
    end
    if mu==3
        Pu=U*U';Pv=V*V';
        PsG=Pu*Grad+Grad*Pv-Pu*Grad*Pv;
        step=norm(PsG,'fro')^2/norm(A(PsG))^2;
    end
    R = Xhat+step*(Grad);
    %---
    [U,sig,V] = lansvd(R+eps*randn(size(R)),r);
    Z=U*sig*V';
    zr = sum(sum(Z.*R));
    rr = norm(R,'fro')^2;
    zz = norm(Z,'fro')^2;
    % calculate alpha
    if divtype==0
        div=0;
        for ti=1:times
            noi=randn(size(R)); Rn=R+epsi*noi;
            [U1,sig1,V1] = lansvd(Rn,r);
            Zn=U1*sig1*V1';
            div=sum(sum((Zn-Z).*noi))/epsi/n+div;
        end
        alpha=div/times;
    else
        if divtype==1
            alpha=(abs(n1-n2)*r-r^2+2.6*(min(n1,n2))*r)/n;
        else
            if(ii<=2)
                div=(abs(n1-n2)*r-r^2+2.3*(min(n1,n2))*r)/n;
                divrange = div-0.01:0.001:div+0.01;
                for al = 1:length(divrange)
                    alpha = divrange(al);
                    c = (zr-alpha*rr)/(zz+alpha^2*rr-2*alpha*zr);
                    err = norm(b-c*(A(Z)-alpha*A(R)))^2;
                    cors(al) = err;
                end
                [~,indx]=  min(abs(cors));
                alpha = divrange(indx);
            end
        end
    end
    
    alpha=max(0.001,min(alpha,0.9));
    % calculation of extrinsic denoiser
    ext=Z-alpha*R;
    c = (zr-alpha*rr)/(zz+alpha^2*rr-2*alpha*zr);
    Xhat=c*ext;
    % next iteration
    % calculate psnr
    mse(ii)=errorfunction(Xhat);
    if mse(ii)<tol
        break;
    end
end
end