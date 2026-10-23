function [xhat,diagnostics]=ict(out,method)
%ICT Synthesis using stored full coefficients, never the reference signal.
% 'dual' is a thresholded discrete canonical-dual inverse.
% 'quadrature' uses Jacobian-weighted adjoint synthesis / Cphi, only
% meaningful for an unnormalized scaledMother bank covering enough scales.
if nargin<2, method='dual'; end
n=out.analysisLength; z=complex(zeros(1,n)); diagnostics=cell(size(out.banks));
for b=1:numel(out.banks)
    bank=out.banks{b}; H=bank.filters; C=bank.coefficients;
    if isempty(H)
        diagnostics{b}=struct('coveredPositiveFraction',0); continue
    end
    weights=bank.synthesisWeights(:);
    numerator=sum(weights.*H.*fft(C,[],2),1);
    S=sum(weights.*abs(H).^2,1);
    valid=S>out.options.InverseTolerance*max(S);
    Z=zeros(1,n);
    switch lower(method)
        case 'dual'
            Z(valid)=numerator(valid)./S(valid);
        case 'quadrature'
            if ~strcmpi(out.options.Family,'scaledMother') || ...
                    ~strcmpi(out.options.Normalization,'none')
                error('CT:Quadrature','Use unnormalized scaledMother filters.');
            end
            Z=numerator/bank.admissibility;
        otherwise
            error('CT:Inverse','Use dual or quadrature.');
    end
    zb=ifft(Z).*exp(-1i*2*pi*bank.carrierHz*out.analysisTime);
    z=z+zb;
    pos=out.frequencyHz>0;
    diagnostics{b}=struct('coveredPositiveFraction',sum(valid & pos)/sum(pos), ...
        'frameResponse',S,'coveredBins',valid);
end
% Resample the entire padded interval before removing the padding.
xfull=real(interpft(z,out.paddedLength,2));
xhat=xfull(out.originalIndices)+out.mean;
end
