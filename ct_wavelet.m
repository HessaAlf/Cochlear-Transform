function [H, info] = ct_wavelet(f, theta, opts)
%CT_WAVELET Fourier response of the first derivative of a gammatone.
% H is K-by-numel(f); f is in Hz in any order. Eq. (12), Appendix 6.3.
% localERB uses lambda(theta)=2*pi*ERB(fc). scaledMother implements
% sqrt(a)*Phi_ref(a*omega) from Appendix 6.4 instead.
if nargin<3, opts=struct; end
if ~isfield(opts,'Order'), opts.Order=4; end
if ~isfield(opts,'Family'), opts.Family='localERB'; end
if ~isfield(opts,'Normalization'), opts.Normalization='none'; end
if ~isfield(opts,'ReferenceHz'), opts.ReferenceHz=cochlear_map(0); end
validateattributes(opts.Order,{'numeric'},{'scalar','integer','>=',2,'<=',16});
validateattributes(theta,{'numeric'},{'real','finite','vector','>=',0,'<=',990});
validateattributes(opts.ReferenceHz,{'numeric'},{'real','finite','scalar','positive'});
fc=cochlear_map(theta(:)); w=2*pi*f(:).';
erb=24.7*(4.37*fc/1000+1); a=opts.ReferenceHz./fc;
H=complex(zeros(numel(fc),numel(w))); gain=ones(numel(fc),1);
for k=1:numel(fc)
    if strcmpi(opts.Family,'localERB')
        wk=w; wc=2*pi*fc(k); lam=2*pi*erb(k); scale=1;
    elseif strcmpi(opts.Family,'scaledMother')
        wk=a(k)*w; wc=2*pi*opts.ReferenceHz;
        lam=2*pi*24.7*(4.37*opts.ReferenceHz/1000+1); scale=sqrt(a(k));
    else
        error('CT:Family','Family must be localERB or scaledMother.');
    end
    % Log form avoids factorial overflow and unnecessary underflow.
    keep=w>0;
    H(k,keep)=scale*exp(gammaln(opts.Order)+log(1i*wk(keep))- ...
        opts.Order*log(lam+1i*(wk(keep)-wc)));
    % Positive-frequency projection is explicit, not a claim that the
    % causal quadrature approximation has zero negative-frequency energy.
    switch lower(opts.Normalization)
        case 'none'
        case 'peak'
            gain(k)=max(abs(H(k,:)));
            if gain(k)>0, H(k,:)=H(k,:)/gain(k); end
        case 'energy'
            gain(k)=sqrt(sum(abs(H(k,:)).^2)*abs(f(2)-f(1)));
            if gain(k)>0, H(k,:)=H(k,:)/gain(k); end
        otherwise
            error('CT:Normalization','Use none, peak, or energy.');
    end
end
info=struct('centerHz',fc,'erbHz',erb,'scale',a,'removedGain',gain, ...
    'family',opts.Family,'normalization',opts.Normalization);
end
