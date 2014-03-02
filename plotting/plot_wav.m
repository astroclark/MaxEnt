% plot_wav.m
%
% Simple plotting script for maxent output
%
% USAGE: summary = plot_wav(maxent_file,fs)
%
%   maxent_file: maxent output file
%   fs: sample frequency
%
% James Clark, james.clark@ligo.org
% 2014-03-01

function summary = plot_wav(maxent_file,opts)

% --------------------------------------
%% Load data
% Loads mat file and extracts quantities.  Note that the reconstructed
% waveform is held in maxEntH (See line 142 of mdcer.m and line 83 of 
% maxentdrvr4.m)

data = load(maxent_file);

fs=opts.fs;
% overlap=opts.overlapFrac;

nsegs = length(data.maxEntH{1});
nsamps = length(data.maxEntH{1}{1});
deltaT = 1/fs;
time = 0:deltaT:nsegs*nsamps*deltaT - deltaT;


% --- Stitch together segments --- 
% FIXME: This needs serious modification for non-zero overlaps
hplus=zeros(1,nsegs*nsamps);
hcross=zeros(1,nsegs*nsamps);
k=1;
for i = 1:nsegs
    hplus(k:nsamps+k-1)  = data.maxEntH{1}{i};
    hcross(k:nsamps+k-1) = data.maxEntH{2}{i};
    k = k+nsamps;
end

% make them time series objects.  This will probably be useful later, there
% are a bunch of handy methods defined on these.
hplus = timeseries(hplus);
hcross = timeseries(hcross);

%% Spectral analysis
% Compute PSD, ...

% can use different estimators here
hspec = spectrum.welch;

% computation
Hplus = psd(hspec, hplus.Data, 'Fs', fs);
Hcross = psd(hspec, hcross.Data, 'Fs', fs);


%% Plotting
% Plots the time series, PSD, ...

%% Populate summary object
% This is just a structure for now, but we could define a useful class for
% this.

summary.hplus = hplus;
summary.hplus = hplus;
summary.hplus_psd = Hplus;
summary.hcross_psd = Hcross;
summary.sngl_snr = snrDet;
summary.net_snr = data.snrComb;

%% Subroutines

    % Get the instantaneous amplitude / frequency with a hilbert transform
    function [amp,freq] = hilbert_content(data, deltaT)
        % Hilbert transform
        h=hilbert(data);
        unrolled_phase = unwrap(angle(h));
        freq = squeeze(diff(unrolled_phase)/deltaT/(2*pi));
        amp = squeeze(abs(h));
    end

end


