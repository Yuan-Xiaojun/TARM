% with different singular value distribution
clear;
warning off;
%--------
r=20;
n2=1000;
n1=n2*1;
n=n1*n2;
p=(n1+n2-r)*r;
rate=0.08;
m=fix(n*rate);
r*(n1+n2-r)/m
ss=(n1+n2-sqrt((n1+n2)^2-4*m))/2;
iter=200;
%----------
t0=0;t1=0;t2=0;t3=0;t4=0;
l0=0;l1=0;l2=0;l3=0;l4=0;
s0=0;s1=0;s2=0;s3=0;s4=0;
num=3;
kk = [0.1,0.5,1];
for ii =1:num
    M=randn(n1,r)*randn(r,n2);
    [Um,~,Vm]=lansvd(M,r);
    k = kk(ii);
    dg = exp(-k*(1:r));
    Sm=diag(dg);
    M=Um*Sm*Vm';
    M=sqrt(n)*M/norm(M,'fro');
    [Um,Sm,Vm]=lansvd(M,r);
    %---low-rank matrix recovery 1:(partial orthogonal case)---
    perm=randperm(n);
    indexs=perm(1:m);
    sign1=2*(rand(n,1)>0.5)-1;
    %-----parameters of TARM-----
    dim.m=m;
    dim.n1=n1;
    dim.n2=n2;
    %----measurement type-----
    A=@(z) subsref(dct(sign1.*z(:)),struct('type','()','subs',{{indexs}}));
    At=@(z) reshape(sign1.*idct(put_vector(n,indexs,z)),size(M));
    %--------
    sigma=0;%sigma=sqrt(10^(-20/10));
    b=A(M)+sigma*randn(m,1);
    %---parameters----
    params.mu=0;
    params.iter=iter;
    params.tol=-100;
    params.divtype= 0;
    error_function = @(qval) 20*log10(norm(qval - M,'fro')/norm(M,'fro'));
    tic;
    [Mhat,psnr] = TARM(b,dim,A,At,r,params,error_function);
    plot(1:length(psnr),psnr,'-o b','LineWidth',1.5);
    hold on;
end
xlabel('Iteration');
ylabel('NMSE (dB)');
legend('k=0.1','k=0.5','k=1');
set(gca,'FontSize',14,'FontName','Times');
grid on;
hold off;