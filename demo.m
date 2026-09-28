fs = 2000;
t = 0:1/fs:2-1/fs;
x = sin(2*pi*50*t) + 0.6*sin(2*pi*120*t)+rand;

params = struct();
params.delta_theta = 20;
params.sigma = 0.1;
params.minPeakSeparationHz = 20;
params.clusterSeparationHz = 10;
params.waveletOrder = 4;
params.tfrFreqBins = 200;
params.verbose = true;
params.doModulation = false;

out = CT_full(x, fs, params);

figure;
subplot(3,1,1);
plot(t, x); title('Input signal'); xlabel('Time (s)');

subplot(3,1,2);
imagesc(out.tfrTimeAxis, out.tfrFreqAxis, out.TFR);
axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('CT-TFR'); colorbar;

subplot(3,1,3);
plot(t, real(out.modes_display.'));
xlabel('Time (s)'); title('Demodulated cochlear modes');