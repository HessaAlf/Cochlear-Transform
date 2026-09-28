function out = CT_full(x, fs, params)
% CT_FULL 

% Complete Cochlear Transform implementation code for H. Alfalahi, A.
% Khandoker and L. Hadjileontiadis, "Cochlear Transform: A Novel
% Bio-inspired Digital Stethoscope for Biomedical Signals", Measurement,
% 2026. 


% Inputs:

% x: input signal 
% fs: sampling frequency
% params: optional struct with fields:
%    .pad = true 
%    .delta_theta = 15 
%    .waveletOrder = 4 
%    .sigma = 0.20 
%    .minPeakSeparationHz = 15 
%    .smoothingSpan = 9
%    .TargetBandHz = [380 3400] 
%    .TargetCenterHz = 1700 
%    .doModulation = true 
%    .doClusterModes = true 
%    .clusterSeparationHz = 20 
%    .normalizeFilters = true 
%    .displayNormalizedModes = true 
%    .tfrFreqBins = 512 
%    .ifClipHz = [0 fs/2] 
%    .verbose = true 

% Output
%   out: struct with fields: 
%    .x -original zero mean signal 
%    .x_mod -modulated signal for analysis 
%    .fs 
%    .params 
%    .theta_all 
%    .centerFreq_all_Hz
%    .selectedPeakFreqs_Hz 
%     .peakFreqs_original_Hz 
%    .selectedTheta_deg 
%    .fm -modulation frequency 
%    .modes_mod -modes in modulated domain 
%    .modes -demodulated modes 
%    .modes_display -display-normalized demodulated modes 
%    .filters -selected cochlear filters in the FFT domain
%    .IA 
%    .IP 
%    .IF 
%    .TFR 
%    .tfrFreqAxis 
%    .tfrTimeAxis 
%    .xhat -approximate reconstructed
%    .residual 
%    .debug 

% Notes: 
% -this implementation used FFT-based cochlear filtering 
% -reconstruction is approximate 
% -For reviewers, keep raw modes for metrics and use modes_display only for
% figures 

if nargin < 3 
    params = struct(); 
end 

params = apply_defaults(params, fs); 

%   -----------prepare signal ------------ % 
x = x(:).'; 
x = x -mean(x); 
n0 = length(x); 

if params.pad 
    nfft = 2^nextpow2(n0); 

    if nfft < n0
        nfft = 2^(nextpow2(n0)+1); 
    end 
    xpad = [x zeros(1, nfft - n0)]; 
else 
    xpad = x; 
    nfft = length(xpad); 
end 

tpad = (0:nfft-1)/fs; 

%  ----------- cochlear map ------------% 

theta_all = 0:params.delta_theta:990; 
centerFreq_all_Hz = theta_to_freq(theta_all); 

cochleaMinHz = theta_to_freq(990); 
cochleaMaxHz = theta_to_freq(0); 

% ----------- cochlear modulation -------%

[x_mod, fm, modInfo] = cochlear_modulation(xpad, fs, params, cochleaMinHz, cochleaMaxHz);

% ---------- PSD peak selection ---------% 
[peakFreqs_Hz, peakInfo] = select_dominant_peaks(x_mod, fs, params);

if isempty(peakFreqs_Hz)
    if params.verbose
        warning('No peaks found with current sigma; falling back to dominant FFT peaks.'); 
    end 
    peakFreqs_Hz = fallback_peak(x_mod, fs); 
end 

% ------- map peaks to theta ---------% 
selectedTheta_deg = freq_to_theta(peakFreqs_Hz); 
selectedTheta_deg = max(0, min(990, selectedTheta_deg)); 

[selectedTheta_deg, ia] = unique(round(selectedTheta_deg, 6), 'stable'); 
peakFreqs_Hz = peakFreqs_Hz(ia); 
peakFreqs_original_Hz = peakFreqs_Hz - fm;

% sort high to low freqs 
[peakFreqs_Hz, idxSort] = sort(peakFreqs_Hz, 'descend'); 
selectedTheta_deg = selectedTheta_deg(idxSort); 


% ------- build selected cochlear wavelet filters ----% 
[filters, omega, fgrid, filterInfo] = build_selected_cochlear_filters(fs, nfft, selectedTheta_deg, params); 

% ------ extract modes (modulated domain)  ---------% 
modes_mod = extract_modes_fft(x_mod, filters); 

% ------ demodulate modes ---------% 
modes = demodulate_modes(modes_mod, fs, fm); 

modes = modes(:, 1:n0); 
modes_mod = modes_mod(:, 1:n0); 

 % -------------------- optional clustering / merging -----
clusterInfo = struct();
if params.doClusterModes && size(modes,1) > 1
   [modes, peakFreqs_Hz, selectedTheta_deg, clusterInfo] = ...
   cluster_close_modes(modes, peakFreqs_Hz, selectedTheta_deg, params.clusterSeparationHz);
end

% -------------------- display normalization -------------
modes_display = modes;
if params.displayNormalizeModes
    modes_display = normalize_modes_for_display(modes, selectedTheta_deg);
end

% -------------------- IA / IP / IF ----------------------
[IA, IP, IF] = analytic_mode_features(modes, fs, params.ifClipHz);

% -------------------- CT-TFR ----------------------------
[TFR, tfrFreqAxis, tfrTimeAxis] = build_ct_tfr(IF, IA, fs, params.tfrFreqBins);

% -------------------- reconstruction --------------------
%xhat = real(sum(modes, 1));
%residual = x - xhat(1:n0);
%xhat_mod = reconstruct_from_modes(x_mod, filters, modes_mod);
%xhat = real(demodulate_signal(xhat_mod, fs, fm));
xhat = reconstruct_modes_lsq(modes,x);
xhat = xhat(1:n0);
residual = x - xhat;
% -------------------- output ----------------------------
out = struct();
out.x = x;
out.x_mod = x_mod(1:n0);
out.fs = fs;
out.params = params;

out.theta_all = theta_all;
out.centerFreq_all_Hz = centerFreq_all_Hz;

out.selectedPeakFreqs_Hz = peakFreqs_Hz;
out.peakFreqs_original_Hz = peakFreqs_original_Hz; 
out.selectedTheta_deg = selectedTheta_deg;
out.fm = fm;

out.modes_mod = modes_mod;
out.modes = modes;
out.modes_display = modes_display;
out.filters = filters;

out.IA = IA;
out.IP = IP;
out.IF = IF;

out.TFR = TFR;
out.tfrFreqAxis = tfrFreqAxis;
out.tfrTimeAxis = tfrTimeAxis;

out.xhat = xhat(1:n0);
out.residual = residual;

out.debug = struct();
out.debug.nfft = nfft;
out.debug.omega = omega;
out.debug.fgrid = fgrid;
out.debug.modInfo = modInfo;
out.debug.peakInfo = peakInfo;
out.debug.filterInfo = filterInfo;
out.debug.clusterInfo = clusterInfo;
end

% ========================================================================
% DEFAULTS
% ========================================================================
function params = apply_defaults(params, fs)
    def.pad                   = true;
    def.delta_theta           = 15;
    def.waveletOrder          = 4;
    def.sigma                 = 0.20;
    def.minPeakSeparationHz   = 15;
    def.smoothingSpan         = 9;
    def.targetBandHz          = [380 3400];
    def.targetCenterHz        = 1700;
    def.doModulation          = true;
    def.doClusterModes        = true;
    def.clusterSeparationHz   = 20;
    def.normalizeFilters      = true;
    def.displayNormalizeModes = true;
    def.tfrFreqBins           = 512;
    def.ifClipHz              = [0 fs/2];
    def.verbose               = true;

    fn = fieldnames(def);
    for k = 1:numel(fn)
        if ~isfield(params, fn{k}) || isempty(params.(fn{k}))
            params.(fn{k}) = def.(fn{k});
        end
    end
end


% ========================================================================
% COCHLEAR THETA <-> FREQ MAP
% ========================================================================
function fHz = theta_to_freq(thetaDeg)
% Frequency-position map from manuscript
    fHz = 165.4 .* (3251^2.1 .* (thetaDeg + 177.3).^(-2.1 .* 1.149) - 0.88);
    fHz = max(fHz, 1e-6);
end

function thetaDeg = freq_to_theta(fHz)
% Inverse of the adopted frequency-position map
% f = 165.4 * (32512.1*(theta+177.3)^(-2.413) - 0.88)
% (theta+177.3)^(-2.413) = ((f/165.4)+0.88)/32512.1
% theta = (...)^(-1/2.413) - 177.3

    exponentVal = -2.1 * 1.149;
    tmp = ((fHz ./ 165.4) + 0.88) ./ 3251^2.1;
    tmp = max(tmp, 1e-16);
    thetaDeg = tmp.^(1/exponentVal) - 177.3;
end

% ========================================================================
% MODULATION
% ========================================================================
function [x_mod, fm, info] = cochlear_modulation(x, fs, params, cochleaMinHz, cochleaMaxHz)

    info = struct();
    x = x(:).';
    n = length(x);
    t = (0:n-1)/fs;

    % Estimate dominant occupied band
    X = fft(x);
    P2 = abs(X).^2 / n;
    P1 = P2(1:floor(n/2)+1);
    f = (0:floor(n/2)) * fs / n;

    % crude occupied-band estimate
    thr = 0.05 * max(P1);
    idx = find(P1 >= thr);

    if isempty(idx)
        bandMin = 0;
        bandMax = fs/2;
    else
        bandMin = f(min(idx));
        bandMax = f(max(idx));
    end

    bandCenter = 0.5*(bandMin + bandMax);
    targetCenter = params.targetCenterHz;

    fm = 0;
    x_mod = x;

    info.bandMin = bandMin;
    info.bandMax = bandMax;
    info.bandCenter = bandCenter;
    info.cochleaMinHz = cochleaMinHz;
    info.cochleaMaxHz = cochleaMaxHz;

    if ~params.doModulation
        info.didModulate = false;
        return;
    end

    % Shift if signal is mostly below or above the cochlear band
    if bandMax < cochleaMinHz || bandMin > cochleaMaxHz
        fm = targetCenter - bandCenter;
        x_analytic = hilbert(x);
        x_mod = real(x_analytic .* exp(1i*2*pi*fm*t));
        info.didModulate = true;
    else
        % Optional gentle shift to target band if center is far outside sensitivity region
        targetBand = params.targetBandHz;
        if bandCenter < targetBand(1) || bandCenter > targetBand(2)
            fm = targetCenter - bandCenter;
            x_analytic = hilbert(x);
            x_mod = real(x_analytic .* exp(1i*2*pi*fm*t));
            info.didModulate = true;
        else
            info.didModulate = false;
        end
    end
end


% ========================================================================
% PEAK SELECTION
% ========================================================================
function [peakFreqs_Hz, info] = select_dominant_peaks(x, fs, params)

    x = x(:).';
    n = length(x);

    X = fft(x);
    P2 = abs(X).^2 / n;
    P1 = P2(1:floor(n/2)+1);
    f = (0:floor(n/2)) * fs / n;

    % very light smoothing only
    span = max(1, round(params.smoothingSpan));
    if span > 1
        kernel = ones(1, span) / span;
        P1s = conv(P1, kernel, 'same');
    else
        P1s = P1;
    end

    P1n = P1s / (max(P1s) + eps);

    minSepBins = max(1, round(params.minPeakSeparationHz / (fs/n)));

    [pks, locs] = findpeaks(P1n, ...
        'MinPeakHeight', params.sigma, ...
        'MinPeakDistance', minSepBins, ...
        'SortStr', 'descend');

    peakFreqs_Hz = f(locs);

    info = struct();
    info.f = f;
    info.PSD = P1;
    info.PSDsmooth = P1s;
    info.PSDnorm = P1n;
    info.peaks = pks;
    info.locs = locs;
end
function peakFreq = fallback_peak(x, fs)
    x = x(:).';
    n = length(x);
    X = fft(x);
    P1 = abs(X(1:floor(n/2)+1)).^2 / n;
    f = (0:floor(n/2)) * fs / n;
    [~, idx] = max(P1(2:end));  % avoid DC if possible
    idx = idx + 1;
    peakFreq = f(idx);
end
% ========================================================================
% FILTER CONSTRUCTION
% ========================================================================
function [filters, omega, fgrid, info] = build_selected_cochlear_filters(fs, nfft, thetaSel, params)

    fgrid = (-nfft/2:nfft/2-1) * (fs/nfft);
    omega = 2*pi*fgrid;

    nModes = numel(thetaSel);
    filters = zeros(nModes, nfft);

    centerFreq_Hz = zeros(1, nModes);
    ERB_Hz = zeros(1, nModes);

    for j = 1:nModes
        fc_Hz = theta_to_freq(thetaSel(j));
        erb_Hz = 24.7 * (4.37 * fc_Hz / 1000 + 1);

        centerFreq_Hz(j) = fc_Hz;
        ERB_Hz(j) = erb_Hz;

        % Analytic Gaussian cochlear filter
        psi = exp(-0.5 * ((fgrid - fc_Hz) / (erb_Hz + eps)).^2);

        % suppress negative frequencies to make it analytic-like
        psi(fgrid < 0) = 0;

        % remove DC
        [~, idx0] = min(abs(fgrid));
        psi(idx0) = 0;

        if params.normalizeFilters
            psi = psi / (norm(psi) + eps);
        end

        filters(j,:) = psi;
    end

    info = struct();
    info.centerFreq_Hz = centerFreq_Hz;
    info.ERB_Hz = ERB_Hz;
end
% ========================================================================
% MODE EXTRACTION
% ========================================================================
function modes = extract_modes_fft(x, filters)
    x = x(:).';
    nfft = size(filters, 2);

    X = fftshift(fft(x, nfft));
    nModes = size(filters, 1);
    modes = zeros(nModes, nfft);

    for j = 1:nModes
        mode_fft = conj(filters(j,:)) .* X;
        mode_t = ifft(ifftshift(mode_fft));
        modes(j,:) = mode_t;
    end
end

% ========================================================================
% DEMODULATION
% ========================================================================


function modes_demod = demodulate_modes(modes_mod, fs, fm)

    if abs(fm) < eps
        modes_demod = modes_mod;
        return;
    end

    [K, N] = size(modes_mod);
    t = (0:N-1)/fs;

    modes_demod = zeros(K, N);

    for k = 1:K
        % preserve complex mode directly
        modes_demod(k,:) = modes_mod(k,:) .* exp(-1i*2*pi*fm*t);
    end
end

% ========================================================================
% OPTIONAL MODE CLUSTERING / MERGING
% ========================================================================
function [modesOut, peakFreqsOut, thetaOut, info] = cluster_close_modes(modes, peakFreqs, thetaDeg, mergeThrHz)

    info = struct();
    info.didCluster = false;

    if numel(peakFreqs) <= 1
        modesOut = modes;
        peakFreqsOut = peakFreqs;
        thetaOut = thetaDeg;
        return;
    end

    [peakFreqsSorted, idx] = sort(peakFreqs, 'descend');
    modesSorted = modes(idx,:);
    thetaSorted = thetaDeg(idx);

    keep = true(size(peakFreqsSorted));
    mergedModes = {};
    mergedFreqs = [];
    mergedTheta = [];

    i = 1;
    while i <= numel(peakFreqsSorted)
        currentMode = modesSorted(i,:);
        currentFreq = peakFreqsSorted(i);
        currentTheta = thetaSorted(i);

        j = i + 1;
        count = 1;

        while j <= numel(peakFreqsSorted) && abs(peakFreqsSorted(j) - currentFreq) < mergeThrHz
            currentMode = currentMode + modesSorted(j,:);
            currentFreq = currentFreq + peakFreqsSorted(j);
            currentTheta = currentTheta + thetaSorted(j);
            keep(j) = false;
            count = count + 1;
            j = j + 1;
        end

        mergedModes{end+1,1} = currentMode; %#ok<AGROW>
        mergedFreqs(end+1,1) = currentFreq / count; %#ok<AGROW>
        mergedTheta(end+1,1) = currentTheta / count; %#ok<AGROW>
        i = j;
    end

    modesOut = cell2mat(mergedModes);
    peakFreqsOut = mergedFreqs(:).';
    thetaOut = mergedTheta(:).';

    info.didCluster = true;
    info.nModesBefore = size(modes,1);
    info.nModesAfter = size(modesOut,1);
end

% ========================================================================
% DISPLAY NORMALIZATION
% ========================================================================
function modesDisp = normalize_modes_for_display(modes, thetaDeg)
    thetaSafe = thetaDeg;
    thetaSafe(thetaSafe == 0) = 1e-3;

    [K, ~] = size(modes);
    modesDisp = zeros(size(modes));

    for k = 1:K
        mag = abs(modes(k,:));
        normFactor = mean(mag) + eps;
        magNorm = (mag / normFactor) * (1/sqrt(thetaSafe(k)));
        modesDisp(k,:) = magNorm .* exp(1i*angle(modes(k,:)));
    end
end


% ========================================================================
% HILBERT IA / IP / IF
% ========================================================================
function [IA, IP, IF] = analytic_mode_features(modes, fs, ifClipHz)

    [K, N] = size(modes);
    IA = zeros(K, N);
    IP = zeros(K, N);
    IF = zeros(K, N);

    for k = 1:K
        z = modes(k,:);

        % if almost purely real, fallback to analytic signal
        if max(abs(imag(z))) < 1e-10
            z = hilbert(real(z));
        end

        IA(k,:) = abs(z);
        IP(k,:) = unwrap(angle(z));

        dphi = [diff(IP(k,:)) 0];
        ifk = (fs/(2*pi)) * dphi;

        if nargin >= 3 && ~isempty(ifClipHz)
            ifk(ifk < ifClipHz(1)) = ifClipHz(1);
            ifk(ifk > ifClipHz(2)) = ifClipHz(2);
        end

        IF(k,:) = ifk;
    end
end


% ========================================================================
% BUILD CT-TFR
% ========================================================================
function [TFR, fAxis, tAxis] = build_ct_tfr(IF, IA, fs, nFreqBins)

    [K, T] = size(IF);
    fAxis = linspace(0, fs/2, nFreqBins);
    tAxis = (0:T-1)/fs;

    TFR = zeros(nFreqBins, T);

    for k = 1:K
        for t = 1:T
            [~, idx] = min(abs(fAxis - IF(k,t)));
            TFR(idx, t) = TFR(idx, t) + abs(IA(k,t));
        end
    end
end

% ========================================================================
% SIGNAL RECONSTRUCTION 
% ========================================================================

function xhat = reconstruct_from_modes(x_ref, filters, modes)
% Approximate reconstruction using filterbank compensation
%
% x_ref   : original analysis-domain signal used before mode extraction
% filters : selected filters in FFT domain (one per row)
% modes   : extracted modes (one per row)
%
% Output:
%   xhat  : reconstructed signal

    x_ref = x_ref(:).';
    n = length(x_ref);
    nfft = size(filters,2);

    Xref = fftshift(fft(x_ref, nfft));

    % Sum of filter responses
    Hsum = sum(abs(filters).^2, 1);

    % Avoid division by very small values
    Hsum(Hsum < 1e-8) = 1e-8;

    % Reconstruct spectrum from filtered modes
    Ysum = zeros(1, nfft);
    for k = 1:size(modes,1)
        Mk = fftshift(fft(modes(k,:), nfft));
        Ysum = Ysum + Mk;
    end

    % Compensation
    Xhat = Ysum ./ Hsum;

    xhat_full = real(ifft(ifftshift(Xhat)));
    xhat = xhat_full(1:n);
end


function y = demodulate_signal(x, fs, fm)
    x = x(:).';
    if abs(fm) < eps
        y = x;
        return;
    end

    N = length(x);
    t = (0:N-1)/fs;
    y = x .* exp(-1i*2*pi*fm*t);
    y = real(y);
end


function xhat = reconstruct_modes_lsq(modes, xref)
% Stable reconstruction from modes using least-squares weights
%
% modes : K x N matrix
% xref  : 1 x N or N x 1 original signal
%
% xhat  : reconstructed signal

    xref = xref(:);
    M = real(modes).';   % N x K

    % Solve least-squares weights
    a = M \ xref;

    % Reconstruct
    xhat = (M * a).';
end
