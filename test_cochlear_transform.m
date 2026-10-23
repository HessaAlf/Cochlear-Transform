function test_cochlear_transform
% Numerical checks of the wavelet, synthesis and feature conventions.
theta=[0 15 200 500 990];
assert(max(abs(cochlear_map(cochlear_map(theta),'angle')-theta))<1e-10);
fs=8000; n=8000; t=(0:n-1)/fs;
o=struct('Modulate',false,'PaddingSeconds',0,'SmoothingHz',0, ...
    'MinPeakSeparationHz',20,'Normalization','peak','Sigma',0.1);
x=2*cos(2*pi*400*t)+0.7*cos(2*pi*1300*t)+3;
r=cochlear_transform(x,fs,o);
assert(numel(r.centerHz)==2);
assert(max(abs(sort(r.centerHz)-[400 1300]))<1e-8);
assert(r.metrics.NRE<1e-10);
assert(max(abs(ict(r)-r.xhat))<1e-12);
% Overlapping filters need not yield perfectly isolated multitone modes.
assert(abs(median(r.IF(1,:))-1300)<1);
single=cochlear_transform(2*cos(2*pi*400*t),fs,o);
assert(max(abs(single.IF(1,:)-400))<1e-7);
bin=find(abs(single.frequencyHz-400)<1e-10,1);
expectedAmplitude=2*abs(single.banks{1}.filters(1,bin));
assert(max(abs(single.IA(1,:)-expectedAmplitude))<1e-8);
% Mutating the stored input must not affect inverse synthesis.
r.x(:)=NaN; assert(all(isfinite(ict(r))));
f=r.frequencyHz; th=cochlear_map([400 1300],'angle');
[H,~]=ct_wavelet(f,th,struct('Order',4));
for j=1:2
    fc=cochlear_map(th(j)); lam=2*pi*24.7*(4.37*fc/1000+1);
    direct=1i*2*pi*f*factorial(3)./(lam+1i*2*pi*(f-fc)).^4;
    direct(f<=0)=0;
    assert(norm(H(j,:)-direct)/norm(direct)<1e-12);
end
assert(all(H(:,f<=0)==0,'all'));
assert(all(abs(sum(ifft(H,[],2),2))<1e-15));
% Finite-difference derivative of the generating function (Appendix 6.3).
tt=linspace(0.001,0.02,10000); fc=400; lam=2*pi*24.7*(4.37*fc/1000+1);
beta=-lam+1i*2*pi*fc; g=tt.^3.*exp(beta*tt);
psi=(3*tt.^2+beta*tt.^3).*exp(beta*tt);
dg=gradient(g,tt(2)-tt(1));
assert(norm(dg(2:end-1)-psi(2:end-1))/norm(psi(2:end-1))<1e-5);
% Automatic modulation must raise analysis fs before shifting low-rate data.
fs2=200; t2=(0:399)/fs2; x2=cos(2*pi*4*t2)+0.7*cos(2*pi*20*t2);
om=struct('PaddingSeconds',0,'SmoothingHz',0,'Normalization','peak', ...
    'Sigma',0.1,'MinPeakSeparationHz',1);
rm=cochlear_transform(x2,fs2,om);
assert(rm.analysisFs>2*max(cellfun(@(b)max(b.centerHz),rm.banks)));
assert(max(abs(sort(rm.centerHz)-[4 20]))<1e-8);
assert(rm.metrics.NRE<1e-8);
% Upper out-of-range content is translated downward and recovered.
fs3=48000; t3=(0:fs3-1)/fs3;
rh=cochlear_transform(cos(2*pi*18000*t3),fs3,om);
assert(rh.banks{1}.carrierHz<0 && rh.metrics.NRE<1e-7);
% Odd lengths, constant signals, empty selection, manual dense grid.
rz=cochlear_transform(7*ones(1,101),200,struct('PaddingSeconds',0));
assert(isempty(rz.modes) && max(abs(rz.xhat-7))==0 && rz.metrics.NRE==0);
ro=cochlear_transform(cos(2*pi*400*(0:1000)/8000),8000, ...
    struct('Theta',0:30:990,'Modulate',false,'ComputeTFR',false));
assert(numel(ro.xhat)==1001 && all(isfinite(ro.xhat)));
% Scaled-mother dual inversion and positive Jacobian quadrature weights.
os=o; os.Family='scaledMother'; os.Normalization='none'; os.Theta=0:15:990;
rs=cochlear_transform(x,fs,os); assert(rs.metrics.NRE<1e-8);
assert(all(rs.banks{1}.synthesisWeights>0));
assert(all(isfinite(ict(rs,'quadrature'))));
% Clustering does not destroy the unmerged analysis coefficients used by ICT.
oc=o; oc.ClusterLabels=[1 1]; rc=cochlear_transform(x,fs,oc);
assert(size(rc.modes,1)==1 && rc.metrics.NRE<1e-10);
[T,~,~]=ct_tfr([1 2;3 4],[100 100;100 NaN],1000,11);
assert(abs(sum(T(:))-6)<1e-12);
fprintf('All Cochlear Transform tests passed.\n');
fprintf('Two-tone NRE: %.3g; modulated NRE: %.3g; high-band NRE: %.3g\n', ...
    r.metrics.NRE,rm.metrics.NRE,rh.metrics.NRE);
end
