function out=cochlear_transform(x,fs,options)
%COCHLEAR_TRANSFORM Paper-based FFT implementation of the Cochlear Transform.
% Alfalahi, Khandoker & Hadjileontiadis, Measurement 292 (2026) 123573.
% Rows of out.modes are complex modes at the ORIGINAL sampling rate.
% See README.md for equation mapping, conventions and limitations.
if nargin<3, options=struct; end
validateattributes(x,{'numeric'},{'real','finite','vector','nonempty'});
validateattributes(fs,{'numeric'},{'real','finite','scalar','positive'});
if numel(x)<4, error('CT:Length','At least four samples required.'); end
o=defaults(options); x=double(x(:).'); n0=numel(x); mu=mean(x); x0=x-mu;
pad=round(o.PaddingSeconds*fs);
if pad>0
    % Reflection with arbitrary padding length, excluding duplicated ends.
    idx=(1-pad):(n0+pad); q=mod(idx-1,2*(n0-1));
    q(q>n0-1)=2*(n0-1)-q(q>n0-1); xp=x0(q+1);
else
    xp=x0;
end
np=numel(xp); original=pad+(1:n0);
za=ct_analytic(xp); P=abs(fft(xp)).^2/np;
fp=(0:floor(np/2))*fs/np; pp=P(1:numel(fp));
occupied=find(pp>=o.OccupancyThreshold*max(pp) & fp>0);
if norm(x0)==0, occupied=[]; end
lo=cochlear_map(990); hi=cochlear_map(0);
splits={za}; carriers=0; bands=[0 fs/2];
if o.Modulate && ~isempty(occupied)
    band=[fp(occupied(1)) fp(occupied(end))]; cut=mean(band);
    if band(1)<lo || band(2)>hi
        f0=(0:np-1)*fs/np; Z=fft(za);
        masks={f0<=cut, f0>cut & f0<=fs/2};
        splits={}; carriers=[]; bands=[];
        for j=1:2
            active=occupied(masks{j}(occupied));
            if isempty(active), continue; end
            bd=[fp(active(1)) fp(active(end))]; fm=0;
            if bd(1)<lo || bd(1)>hi
                if diff(bd)>diff(o.TargetBandHz)
                    error('CT:Bandwidth','Branch exceeds target band; split it externally or disable modulation.');
                end
                fm=mean(o.TargetBandHz)-mean(bd);
            end
            splits{end+1}=ifft(Z.*masks{j}); %#ok<AGROW>
            carriers(end+1)=fm; bands(end+1,:)=bd; %#ok<AGROW>
        end
    end
end
% Preserve original Nyquist content even when only one branch was occupied.
% Spectral occupancy is a selection heuristic; omitted weak branches are
% reflected in the reported residual, not secretly added to reconstruction.
maxHz=fs/2;
if ~isempty(carriers), maxHz=max([maxHz; bands(:,2)+carriers(:)]); end
L=max(1,ceil(2*maxHz*1.1/fs));
if all(carriers==0), L=1; end
fa=fs*L; na=np*L; ta=(0:na-1)/fa;
f=(0:na-1)*fa/na; f(f>fa/2)=f(f>fa/2)-fa;
if rem(na,2)==0, f(na/2+1)=-fa/2; end
banks=cell(1,numel(splits)); zAll=complex(zeros(0,n0));
centers=[]; angles=[]; branches=[];
for b=1:numel(splits)
    zm=interpft(splits{b},na,2).*exp(1i*2*pi*carriers(b)*ta);
    [peaks,peakInfo]=select_peaks(zm,fa,o,[lo min(hi,fa/2)]);
    if ~isempty(o.Theta)
        theta=sort(unique(o.Theta(:).')); peaks=cochlear_map(theta);
        valid=peaks<fa/2; theta=theta(valid); peaks=peaks(valid);
    else
        theta=cochlear_map(peaks,'angle');
    end
    if isempty(theta)
        H=zeros(0,na); info=struct;
    else
        [H,info]=ct_wavelet(f,theta,o);
    end
    C=ifft(conj(H).*fft(zm),[],2);
    zd=C.*exp(-1i*2*pi*carriers(b)*ta);
    zd=interpft(zd,np,2); zd=zd(:,original);
    zAll=[zAll; zd]; %#ok<AGROW>
    centers=[centers peaks-carriers(b)]; angles=[angles theta]; %#ok<AGROW>
    branches=[branches b*ones(1,numel(theta))]; %#ok<AGROW>
    weights=ones(numel(theta),1); admissibility=NaN;
    if strcmpi(o.Family,'scaledMother') && ~isempty(theta)
        % Trapezoidal integration in theta with |a'|/a^2=|fc'|/fref.
        if numel(theta)>1
            d=diff(theta); quadrature=[d(1)/2 (d(1:end-1)+d(2:end))/2 d(end)/2];
        else
            quadrature=1;
        end
        df=165.4*3251^2.1*(2.1*1.149)*(theta+177.3).^(-2.1*1.149-1);
        weights=(quadrature.*df/o.ReferenceHz).';
        wc=2*pi*o.ReferenceHz; lam=2*pi*24.7*(4.37*o.ReferenceHz/1000+1);
        % Dimensionless integration keeps the peak near v=1 even for
        % a high-frequency reference and very small raw wavelet gains.
        integrand=@(v) v./((lam/wc)^2+(v-1).^2).^o.Order;
        admissibility=exp(2*gammaln(o.Order)+(2-2*o.Order)*log(wc))* ...
            (integral(integrand,0,1,'RelTol',1e-8)+ ...
             integral(integrand,1,Inf,'RelTol',1e-8));
    end
    banks{b}=struct('filters',H,'coefficients',C,'carrierHz',carriers(b), ...
        'theta',theta,'centerHz',peaks,'waveletInfo',info, ...
        'peakInfo',peakInfo,'synthesisWeights',weights,'admissibility',admissibility);
end
[~,ix]=sort(centers,'descend'); zAll=zAll(ix,:);
centers=centers(ix); angles=angles(ix); branches=branches(ix);
[modes,groups]=ct_cluster(zAll,centers,o.ClusterSeparationHz,o.ClusterLabels);
[IA,IP,IF]=ct_features(modes,fs,o.PhaseFloor);
out=struct('x',x,'mean',mu,'fs',fs,'time',(0:n0-1)/fs,'options',o, ...
    'analysisFs',fa,'analysisLength',na,'analysisTime',ta,'frequencyHz',f, ...
    'paddedLength',np,'originalIndices',original,'banks',{banks}, ...
    'channels',zAll,'centerHz',centers,'theta',angles,'branch',branches, ...
    'modes',modes,'groups',{groups},'IA',IA,'IP',IP,'IF',IF);
if o.ComputeTFR
    [out.TFR,out.tfrFrequencyHz,out.tfrTime]=ct_tfr(IA,IF,fs,o.FrequencyBins,[0 fs/2]);
end
[out.xhat,out.inverseDiagnostics]=ict(out);
out.residual=x-out.xhat;
out.metrics=ct_metrics(x,out.xhat,modes);
end

function o=defaults(user)
o=struct('Order',4,'Sigma',0.2,'SmoothingHz',1,'MinPeakSeparationHz',1, ...
    'MaxModes',128,'PaddingSeconds',0.1,'Modulate',true, ...
    'TargetBandHz',[380 3400],'OccupancyThreshold',1e-6,'Theta',[], ...
    'Family','localERB','ReferenceHz',cochlear_map(0),'Normalization','none', ...
    'ClusterSeparationHz',0,'ClusterLabels',[],'PhaseFloor',1e-8, ...
    'InverseTolerance',1e-12,'ComputeTFR',true,'FrequencyBins',512);
names=fieldnames(user);
for k=1:numel(names)
    if ~isfield(o,names{k}), error('CT:Option','Unknown option: %s',names{k}); end
    o.(names{k})=user.(names{k});
end
validateattributes(o.Sigma,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(o.PaddingSeconds,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(o.SmoothingHz,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(o.MinPeakSeparationHz,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(o.MaxModes,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.InverseTolerance,{'numeric'},{'scalar','real','>',0,'<',1});
validateattributes(o.OccupancyThreshold,{'numeric'},{'scalar','real','>',0,'<=',1});
validateattributes(o.TargetBandHz,{'numeric'},{'vector','numel',2,'real','finite','positive','increasing'});
if ~isempty(o.Theta), validateattributes(o.Theta,{'numeric'},{'vector','real','finite','>=',0,'<=',990}); end
end

function [peaks,info]=select_peaks(z,fs,o,limits)
n=numel(z); f=(0:floor(n/2))*fs/n; P=abs(fft(z)).^2/(fs*n);
P=P(1:numel(f)); span=max(1,round(o.SmoothingHz*n/fs));
P=conv(P,ones(1,span)/span,'same');
valid=f>=limits(1) & f<=limits(2); P(~valid)=0;
if max(P)>0, P=P/max(P); end
% Strict comparison on one side selects a unique point on a plateau.
idx=find(P>=[Inf P(1:end-1)] & P>[P(2:end) Inf] & P>=o.Sigma & P>0);
[~,order]=sort(P(idx),'descend'); selected=[];
for j=order
    if isempty(selected) || all(abs(f(idx(j))-f(selected))>=o.MinPeakSeparationHz)
        selected(end+1)=idx(j); %#ok<AGROW>
        if numel(selected)>=o.MaxModes, break; end
    end
end
peaks=sort(f(selected),'descend');
info=struct('frequencyHz',f,'normalizedPSD',P,'selectedHz',peaks,'smoothingBins',span);
end
