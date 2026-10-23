function [IA,IP,IF]=ct_features(z,fs,relativeFloor)
%CT_FEATURES Eqs. (24)-(27). IF is NaN where phase is unreliable.
if nargin<3, relativeFloor=1e-8; end
IA=abs(z); IP=unwrap(angle(z),[],2); n=size(z,2);
IF=nan(size(IA));
if n>=2
    IF(:,1)=(IP(:,2)-IP(:,1))*fs/(2*pi);
    IF(:,end)=(IP(:,end)-IP(:,end-1))*fs/(2*pi);
end
if n>=3, IF(:,2:end-1)=(IP(:,3:end)-IP(:,1:end-2))*fs/(4*pi); end
if ~isempty(z)
    IF(IA<=relativeFloor*max(IA,[],2))=NaN;
end
end
