function [wirf , res] = WoldBoot(y,p,H,c)
if nargin<7
    opt=0;
end
n=size(y,2); 

%_____________________________________
% OLS estimates
%_____________________________________
[Y X] = VarStr(y,c,p);       % yy and XX all the sample

T=size(Y,1);
Bols=inv(X'*X)*X'*Y;
VecB=reshape(Bols,n*(n*p+1),1);
B=companion(VecB,p,n,c);
C=Bols(1,:)';
res=Y-X*Bols;

% impulse response functions
for h=1:H
    irf(:,:,h)=B^(h-1); %companion form
    wirf(:,:,h)=irf(1:n,1:n,h); %Reduce form IRF. Point estimate
end
