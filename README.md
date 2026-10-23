# Cochlear Transform: MATLAB implementation

Based on H. Alfalahi, A. Khandoker and L. Hadjileontiadis,
"Cochlear transform: A novel bio-inspired digital stethoscope for biomedical
signals," Measurement 292 (2026), 123573, DOI 10.1016/j.measurement.2026.123573.
Source: the PDF supplied in the parent folder. This is an independently
written implementation of the published equations, not a recovery of the
authors' original MATLAB code or a claim to reproduce all published figures.
The existing files in the parent folder are preserved.

## Run

In MATLAB, change to this folder and run:

```matlab
test_cochlear_transform
out = demo_cochlear_transform;
```

For your own single-channel recording:

```matlab
addpath('cochlear_transform')
opts = struct('Order',4,'Sigma',0.05,'SmoothingHz',1, ...
              'MinPeakSeparationHz',2,'Normalization','none');
out = cochlear_transform(x,fs,opts);
xhat = ict(out);                  % does not require out.x
plot(out.time,x,out.time,xhat)
imagesc(out.tfrTime,out.tfrFrequencyHz,out.TFR); axis xy
```

Requires base MATLAB with implicit expansion (R2016b or later). No Signal
Processing, Wavelet, Statistics, or Text Analytics Toolbox is required.
Tests were prepared for the installed R2022b; see validation.log for the
actual execution results. The functions do not call any parent-folder code.

The complete numerical suite passed on MATLAB R2022b Update 4. Measured
NRE was 3.57e-13 for the in-band two-tone case, 2.13e-11 for the low-rate
modulation case, and 4.10e-11 for the translated high-band case. MATLAB
Code Analyzer reported no issues in the supplied .m files. These are
controlled numerical checks, not replications of the paper's experiments.

## Equation-to-code map

| Paper | Function / implementation |
|---|---|
| Eq. (1), cochlear position-frequency map | `cochlear_map`, plus its exact algebraic inverse |
| Eqs. (6), (12), Appendix 6.3 | `ct_wavelet`: `j*w*Gamma(N)/(lambda+j*(w-wc))^N` |
| Eq. (14), ERB | `24.7*(4.37*fc/1000+1)`; `lambda=2*pi*ERB` |
| Eqs. (19)-(20), (23) | FFT analytic signal, complex spectral translation and demodulation |
| Spectral selection, Section 2.3 | Normalized smoothed periodogram, local maxima, sigma and minimum spacing |
| Eq. (21), correlation coefficients | `ifft(fft(analyticSignal).*conj(H))` |
| Eqs. (24)-(27) | `ct_features`: modulus, unwrapped atan2 phase, finite difference IF |
| Eq. (29) | `ct_tfr`: amplitude accumulated at the nearest IF bin |
| Appendix 6.1 | `ict`: discrete canonical-dual frame synthesis with a coverage threshold |
| Eqs. (62)-(74) | `Family='scaledMother'`; scale `a=fref/fc`, response `sqrt(a)*Phi_ref(a*w)` |
| Eqs. (70)-(71) | Trapezoid angular weights times `abs(fc'(theta))/fref`; optional quadrature synthesis |

The continuous Fourier convention is `F(w)=integral f(t)*exp(-j*w*t) dt`.
H samples this continuous response on FFT frequencies. FFT multiplication
and IFFT implement the corresponding circular discrete filtering operator;
an extra `1/fs` is not applied to this frequency-response representation.
Coefficients use Eq. (21)/Appendix (64) inner-product normalization. The
extra normalization appearing in Eq. (2) is not applied a second time.

## Explicit choices where the paper is incomplete or inconsistent

1. **Wavelet formula.** The denominator means the entire complex quantity
   raised to N. The time-domain generating function is
   `g(t)=t^(N-1)*exp((-lambda+j*wc)*t)*u(t)`.
   Its first derivative is
   `((N-1)*t^(N-2)+(-lambda+j*wc)*t^(N-1))*exp((-lambda+j*wc)*t)`
   for positive t. This is the internally consistent derivative construction
   from Eq. (12) and Appendix 6.3; no additional cosine is multiplied into
   that complex derivative. N must be at least 2. One vanishing moment is
   guaranteed; order N is not a count of vanishing moments.
2. **Analytic projection.** The causal complex gammatone is a quadrature
   approximation, with a nonzero negative-frequency tail. The implementation
   explicitly sets H(f<=0)=0, giving exactly analytic discrete coefficients
   and exactly zero DC. This changes the strictly causal impulse response;
   it does not have compact spectral support at positive frequencies.
   The even-length Nyquist bin is omitted, so exact recovery of arbitrary
   Nyquist-only signals is not promised.
3. **Two wavelet families.** Default `localERB` recomputes ERB at each center
   frequency, as in Section 2.2. ERB has an additive constant, so these are
   not exact dilations of one mother wavelet. `scaledMother` implements the
   appendix dilation family separately; it does not also recompute local
   ERB. This distinction matters for the continuous inverse derivation.
4. **Adaptive selection.** Sigma is a fraction of normalized *power*, not
   a fraction of FFT magnitude. Smoothing is an average of the PSD with a
   width specified in Hz. Peaks are selected strongest first under a minimum
   separation and then ordered by decreasing frequency. There is no forced
   fallback peak for a zero or unrepresented signal. `Theta` bypasses
   selection and requests a dense/manual grid in the physical range 0:990.
   Out-of-range centers are never silently clamped to a cochlear endpoint.
5. **Modulation.** An occupied band is estimated relative to its maximum
   power. If it extends outside the physical cochlear map, analytic Fourier
   components are split at the occupied band's midpoint. Each nonempty
   branch is shifted if its minimum frequency lies below the map or above
   it. Shifted branch centers go to the midpoint of `TargetBandHz`.
   Wide branches that cannot fit the target band raise an error.
   This is an explicit ideal spectral split, since filter order and cutoff
   conventions for the paper's low/high split are unspecified. Weak branches
   below the occupancy threshold can be omitted; their error remains visible.
   A branch straddling the upper map limit can retain unrepresented content;
   use a suitable external band split for broadband signals spanning that
   boundary. Disable modulation for direct in-band analysis.
6. **Sampling.** Before modulation, Fourier interpolation raises sampling
   by an integer factor with a Nyquist guard, and inverse Fourier resampling
   returns results to the input grid. Merely multiplying 200 Hz data by a
   1700 Hz carrier would alias and is not used. Fourier interpolation assumes
   a periodic padded interval; reflection padding mitigates endpoint jumps.
   Arbitrary non-bin-centered carrier shifts can still cause edge errors.
7. **Clustering.** The paper gives no reproducible distance, linkage, or
   harmonic rule. Clustering is off by default. `ClusterSeparationHz>0`
   joins connected adjacent centers. `ClusterLabels` explicitly groups the
   selected channels (in the output's descending original-frequency order).
   For example, `[1 2 2]` keeps channel 1 and sums channels 2 and 3.
   Unmerged coefficients remain available for synthesis. Automatic harmonic
   source separation, including the two sawtooth sources in Fig. 5, cannot
   be inferred uniquely from the paper and is not claimed.
8. **Reconstruction.** The default is a thresholded canonical-dual inverse:
   `Zhat=sum(w.*H.*FFT(C))./sum(w.*abs(H).^2)` on covered bins. It is a
   discrete frame inverse, not a literal evaluation of the infinite-scale
   continuous integral. No reference-fitted least-squares weights or input
   residual are inserted. A stored scalar mean is restored separately.
   The optional `ict(out,'quadrature')` uses the finite-grid Jacobian-weighted
   adjoint divided by the numerically integrated positive-frequency
   admissibility constant. It requires unnormalized `scaledMother` responses
   and is approximate on a finite grid. A sparse adaptive bank is not
   automatically tight, orthogonal, or complete. Coverage is reported.
9. **Amplitude and IF.** Default raw responses use k=1, since k is not
   specified by the paper. `peak` and `energy` normalization are explicit
   alternatives; they change coefficient amplitudes and low-frequency gain.
   Use a consistent normalization for comparisons. IF is NaN at negligible
   amplitude and is not clipped; out-of-range trajectories are excluded from
   TFR accumulation. A central phase difference is used internally and a
   one-sided difference at each endpoint, avoiding a spurious final zero.
   The derivative factor shifts the maximum response slightly above the
   gammatone carrier fc; peak normalization does not imply unit gain at fc.
10. **Metrics.** NRE, RMSE, SNR and mode correlation are included. OI is
    the absolute off-diagonal mode cross-energy normalized by total input
    energy, avoiding pointwise division by zero in Eq. (30). It is a stated
    stable diagnostic convention, not a promise of the table's exact OI.

## Options

| Option | Default | Meaning |
|---|---:|---|
| Order | 4 | Gammatone order, integer 2 through 16 |
| Sigma | 0.2 | Minimum normalized power for a local peak |
| SmoothingHz | 1 | PSD averaging width; 0 disables smoothing |
| MinPeakSeparationHz | 1 | Minimum selected center spacing in analysis domain |
| MaxModes | 128 | Maximum channels per branch |
| PaddingSeconds | 0.1 | Reflection on each side; 0 gives a periodic interval |
| Modulate | true | Enable branch translation when outside cochlear range |
| TargetBandHz | [380 3400] | Target range for shifted branch |
| OccupancyThreshold | 1e-6 | Relative power used for band estimation |
| Theta | [] | Manual angle grid, in degrees; [] enables adaptive selection |
| Family | localERB | localERB or scaledMother |
| ReferenceHz | f(0) | Reference frequency of scaled mother |
| Normalization | none | none, peak, or energy |
| ClusterSeparationHz | 0 | Automatic proximity grouping; 0 disables |
| ClusterLabels | [] | Explicit group label per output channel |
| PhaseFloor | 1e-8 | Per-mode relative amplitude floor for IF validity |
| InverseTolerance | 1e-12 | Relative frame-response cutoff for inversion |
| ComputeTFR | true | Allocate dense frequency-by-time TFR |
| FrequencyBins | 512 | Number of TFR bins over original 0:fs/2 |

## Outputs and practical limits

`channels` contains unmerged complex channels on the original time grid.
`centerHz`, `theta`, and `branch` correspond to those rows. `modes` contains
the optionally grouped rows; `groups` maps them back to channels.
`IA`, `IP`, and `IF` correspond to grouped modes. Frequencies in `centerHz`
and IF are demodulated to the original domain. Each `banks{b}` stores its
analysis-domain filters, full padded coefficients, carrier, and peak details.
Filters use unshifted FFT order, indicated by `frequencyHz`. `xhat` and
`residual` retain the original mean. `inverseDiagnostics` reports frame
responses and covered frequency bins. A small reconstruction error does not
establish accurate source separation: overlapping filters can reconstruct
well while individual channels still contain several oscillations.

Memory is O(K*Nanalysis + FrequencyBins*Ninput); for long recordings, turn
off ComputeTFR, reduce the angle grid / MaxModes, or process windows.
This is an offline implementation, not a streaming causal stethoscope.

The paper's first synthetic signal can be entered literally as:

```matlab
fs=200; t=0:1/fs:1;
x=6*t+cos(8*pi*t)+0.5*cos(40*pi*t);
r=cochlear_transform(x,fs,struct('Sigma',0.05));
```

The linear trend is not a zero-frequency constant. It remains in the
unrepresented residual to the extent not recovered by selected wavelets;
mean restoration alone does not recover a trend. The supplied demonstration
uses two periodic tones so frequency and amplitude conventions can be
checked without that ambiguity. Original experimental recordings,
unreported clustering settings and original source code would be needed
to establish numerical replication of the paper's figures and tables.
