function m=ct_metrics(x,xhat,modes)
%CT_METRICS Reconstruction and global cross-energy orthogonality diagnostics.
e=x(:)-xhat(:); nx=norm(x(:)); ne=norm(e);
if nx==0, nre=double(ne~=0)*Inf; if ne==0, nre=0; end
else, nre=ne/nx; end
if ne==0, snr=Inf; elseif nx==0, snr=-Inf; else, snr=20*log10(nx/ne); end
M=real(modes); G=M*M.'; den=sum(x(:).^2);
off=G-diag(diag(G));
m=struct('NRE',nre,'RMSE',sqrt(mean(e.^2)),'SNRdB',snr, ...
    'OI',sum(abs(off(:)))/max(den,realmin));
if size(M,1)>0
    Mc=M-mean(M,2); V=Mc*Mc.'; d=sqrt(diag(V));
    m.correlation=V./(d*d.');
else
    m.correlation=zeros(0);
end
end
