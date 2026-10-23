function [T,f,t]=ct_tfr(IA,IF,fs,nBins,rangeHz)
%CT_TFR Eq. (29), amplitude accumulated in nearest frequency bin.
% Values are amplitude per bin (not an energy density).
if nargin<4, nBins=512; end
if nargin<5, rangeHz=[0 fs/2]; end
validateattributes(nBins,{'numeric'},{'scalar','integer','>=',2});
assert(isequal(size(IA),size(IF)),'IA and IF sizes must agree.');
f=linspace(rangeHz(1),rangeHz(2),nBins).'; t=(0:size(IA,2)-1)/fs;
T=zeros(nBins,numel(t));
for k=1:size(IA,1)
    good=isfinite(IF(k,:)) & IF(k,:)>=f(1) & IF(k,:)<=f(end);
    cols=find(good); bins=1+round((IF(k,good)-f(1))/(f(end)-f(1))*(nBins-1));
    ix=sub2ind(size(T),bins,cols); T(ix)=T(ix)+IA(k,good);
end
end
