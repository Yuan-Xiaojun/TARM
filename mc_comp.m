% test matrix completion
clear;
warning off;
addpath('./LRGeomCG');
%--------
r=20;
n2=1000;
n1=n2*1;
n=n1*n2;
p=(n1+n2-r)*r;
rate=0.12;
m=fix(n*rate);
r*(n1+n2-r)/m
dof=(n1+n2-sqrt((n1+n2)^2-4*m))/2;
iter=40;
%--
t0=0;t1=0;t2=0;t3=0;t4=0;t5=0;t6=0;
l0=0;l1=0;l2=0;l3=0;l4=0;l5=0;l6=0;
s0=0;s1=0;s2=0;s3=0;s4=0;s5=0;s6=0;
num=1;

for ii =1:num
    M=randn(n1,r)*randn(r,n2);
    %[Um,~,Vm]=lansvd(M,r);% Sm=diag(diag(Sm)-min(diag(Sm))+0.001);
    %dg = exp(-0.5.*(1:r));
    %dg = randn(1,r);
    %Sm=diag(dg);
    %M=Um*Sm*Vm';
    M=sqrt(n)*M/norm(M,'fro');
    %------mc-----
    perm=randperm(n);
    indexs=perm(1:m);
    sigma=0; %
    %sigma=sqrt(10^(-50/10));
    w=sigma*randn(m,1);
    error_function = @(qval) 20*log10(norm(qval - M,'fro')/norm(M,'fro'));
    %----------
    dim.m=m;
    dim.n1=n1;
    dim.n2=n2;
    A=@(z) subsref(z(:),struct('type','()','subs',{{indexs}}));
    At=@(z) reshape(put_vector(n,indexs,z),[n1,n2]);
    b=A(M)+w;
    tol = -60;
    %------LRGeomCG----
    params.iter=iter; 
    params.tol=tol;
    tic;
    [~,mse6] = LRG(b,dim,A,At,r,params,error_function);
    if(mse6(length(mse6))<tol)
        t6=t6+toc;
        l6=l6+length(mse6);
        s6=s6+1;
    end
    %--------------ALPS-------------
    aplsparams.tol = 1e-10;
    aplsparams.tol2 = tol;
    aplsparams.xpath = 1;
    aplsparams.svdMode = 'propack';
    aplsparams.ALPSiters = iter;
    aplsparams.svdApprox = 0;
    aplsparams.cg_tol = 1e-10;
    aplsparams.cg_maxiter = 500;
    tic;
    [X_hat, numiter, mse3] = matrixALPSII(b, A, At, n1, n2, r, aplsparams, error_function);
    if(mse3(length(mse3))<tol)
        t3=t3+toc;
        l3=l3+length(mse3);
        s3=s3+1;
    end
    %-----------LMAFit-------------
    opts.tol2=tol;
    opts.maxit=iter;
    opts.A = A;
    opts.At =At;
    opts.Zfull = 0;
    tic;
    [X,Y,Out] = lmafit_mc_adp(n1,n2,r,indexs,b',opts,error_function);
    mse1 = Out.psnr;
    if(mse1(length(mse1))<tol)
        t1=t1+toc;
        l1=l1+Out.iter;
        s1=s1+1;
    end
    %---TARM parameters----
    params.mu=3; %0: 1/delta,1: auto tuning type 1,2: auto-tuning 2 3: auto-tuning 3
    params.iter=iter; % max iteration time
    params.tol=tol; % tol for stopping
    params.divtype=1; %0: simulation, 1: approximation
    params.sigma=sigma;
    params.ptype = 'MC';
    tic;
    [~,mse0] = TARM(b,dim,A,At,r,params,error_function);
    if(mse0(length(mse0))<tol)
        s0=s0+1;
        t0=t0+toc;
        l0=l0+length(mse0);
    end
    %------RCG------
    tic;
    [~,mse2]=RCG(b,dim,A,At,r,params,error_function);
    if(mse2(length(mse2))<tol)
        t2=t2+toc;
        s2=s2+1;
        l2=l2+length(mse2);
    end
    %------RGrad------
    tic;
    [~,mse5]=RGrad_IV(b,dim,A,At,r,params,error_function);
    if(mse5(length(mse5))<tol)
        t5=t5+toc;
        s5=s5+1;
        l5=l5+length(mse5);
    end
    %------NIHT------
    tic;
    [~,mse4,timedata] = NIHT(b,dim,A,At,r,params,error_function);
    if(mse4(length(mse4))<tol)
        t4=t4+toc;
        s4=s4+1;
        l4=l4+length(mse4);
    end
end

plot(1:length(mse0),mse0,'-p b',1:length(mse1),mse1,'-. r',1:length(mse2),mse2,'-> g',1:length(mse3),mse3,'-v k',1:length(mse4),mse4,'--o',1:length(mse5),mse5,'--p',1:length(mse6),mse6,'- <','LineWidth',1.5);
xlabel('Iteration');
ylabel('NMSE (dB)');
legend('TARM','LMAFit','RCG','ALPS','NIHT','RGrad','LRGeomCG');
set(gca,'FontSize',14,'FontName','Times');
grid on;

%---TARM---
if(s0>0)
    fprintf('TARM-iter:%f \n', l0/s0);
    fprintf('TARM:%f \n', s0);
    fprintf('TARM-time:%f \n', t0/s0);
end
%------LMAfit------
if(s1>0)
    fprintf('LMAFit-iter:%f \n', l1/s1);
    fprintf('LMAFit:%f \n', s1);
    fprintf('LMAFit-time:%f \n', t1/s1);
end
%------RCG------
if(s2>0)
    fprintf('RCG-iter:%f \n', l2/s2);
    fprintf('RCG:%f \n', s2);
    fprintf('RCG-time:%f \n', t2/s2);
end
%------ALPS------
if(s3>0)
    fprintf('ALPS-iter:%f \n', l3/s3);
    fprintf('ALPS:%f \n', s3);
    fprintf('ALPS-time:%f \n', t3/s3);
end
%-----NIHT----
if(s4>0)
    fprintf('NIHT-iter:%f \n', l4/s4);
    fprintf('NIHT:%f \n', s4);
    fprintf('NIHT-time:%f \n', t4/s4);
end
%-----RGrad---
if(s5>0)
    fprintf('RGrad-iter:%f \n', l5/s5);
    fprintf('RGrad:%f \n', s5);
    fprintf('RGrad-time:%f \n', t5/s5);
end
%-----LRGeomCG---
if(s6>0)
    fprintf('LRGeomCG-iter:%f \n', l6/s6);
    fprintf('LRGeomCG:%f \n', s6);
    fprintf('LRGeomCG-time:%f \n', t6/s6);
end
