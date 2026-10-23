function z=ct_analytic(x)
%CT_ANALYTIC Analytic signal along rows; no toolbox required.
n=size(x,2); h=zeros(1,n); h(1)=1;
if rem(n,2)==0
    h(2:n/2)=2; h(n/2+1)=1;
else
    h(2:(n+1)/2)=2;
end
z=ifft(fft(x,[],2).*h,[],2);
end
