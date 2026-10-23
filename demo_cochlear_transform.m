function out=demo_cochlear_transform(makePlots)
%DEMO_COCHLEAR_TRANSFORM Run a reproducible example; no toolbox required.
if nargin<1, makePlots=true; end
fs=200; t=(0:399)/fs;
x=cos(2*pi*4*t)+0.5*cos(2*pi*20*t);
opts=struct('Sigma',0.05,'SmoothingHz',0,'PaddingSeconds',0, ...
    'Normalization','peak','FrequencyBins',401);
out=cochlear_transform(x,fs,opts);
fprintf('Analysis fs: %g Hz; %d channels; NRE: %.4g\n', ...
    out.analysisFs,size(out.channels,1),out.metrics.NRE);
disp(table(out.centerHz(:),out.theta(:),'VariableNames',{'OriginalHz','ThetaDegrees'}));
if ~makePlots, return; end
figure('Color','w','Name','Cochlear Transform','Position',[100 100 1000 800]);
tiledlayout(3,1,'Padding','loose','TileSpacing','loose');
nexttile; plot(t,x,'k',t,out.xhat,'--r');
xlabel('Time (s)'); ylabel('Amplitude'); legend('Input','Dual reconstruction');
nexttile; plot(t,real(out.modes).'); xlabel('Time (s)'); ylabel('Real modes');
nexttile; imagesc(out.tfrTime,out.tfrFrequencyHz,out.TFR); axis xy;
ylim([0 30]); xlabel('Time (s)'); ylabel('Frequency (Hz)'); colorbar;
title('CT-TFR: instantaneous amplitude per frequency bin');
end
