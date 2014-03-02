% plot_wav.m
%
% Simple plotting script for maxent output
%
% USAGE: summary = plot_wav(maxent_file,fs)
%
%   maxent_file: maxent output file
%   fft_config:  structure with FFT configuration.  Only requires:
%   fft_config.fs for now.
%
% James Clark, james.clark@ligo.org
% 2014-03-01

function summary = plot_wav(maxent_file,fft_config)

% --------------------------------------
%% Load data
% Loads mat file and extracts quantities.  Note that the reconstructed
% waveform is held in maxEntH (See line 142 of mdcer.m and line 83 of 
% maxentdrvr4.m)

% argh, why no strsplit??
idx=strfind(maxent_file,'/');
outdir=maxent_file(idx(end)+1:end);
outdir=strrep(outdir,'.mat','');
mkdir('./',outdir); % this going to create a new directory in the PWD.
                    % Should prob make the parent an argument for easy
                    % condor deployment

data = load(maxent_file);

fs=fft_config.fs;
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
hplus.timeseries = timeseries(hplus);
hcross.timeseries = timeseries(hcross);

% --------------------------------------
%% Spectral analysis
% Compute PSD, ...

% can use different estimators here
hspec = spectrum.welch;

% computation
hplus.psd  = psd(hspec, hplus.timeseries.Data, 'Fs', fs);
hcross.psd = psd(hspec, hcross.timeseries.Data, 'Fs', fs);

% Hilbert
[hplus.instamp, hplus.instfreq] = ...
    hilbert_content(hplus.timeseries.Data, deltaT);

[hcross.instamp, hcross.instfreq] = ...
    hilbert_content(hcross.timeseries.Data, deltaT);

% --------------------------------------
%% Plotting
% Plots the time series, PSD, ...

% --- Time series
f=figure('Position', [10, 50, 900, 200], 'visible','off');
clf;
subplot(121)
plot(hplus.timeseries.Time,squeeze(hplus.timeseries.Data))

subplot(122)
plot(hcross.timeseries.Time,squeeze(hcross.timeseries.Data))

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], ...
    'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])

saveas(f, [outdir,'/timeseries'], 'm')
exportfig([outdir,'/timeseries.png'],'width',8,'height',4,'color','rgb')
close(f)

% TODO: labels, consistent axis limits etc

% --- PSD
f=figure('Position', [10, 50, 900, 400], 'visible','off');
clf;
subplot(121)
plot(hplus.psd.Frequencies, hplus.psd.Data)

subplot(122)
plot(hcross.psd.Frequencies, hcross.psd.Data)

saveas(f, [outdir,'/psd'], 'm')
exportfig([outdir,'/psd.png'],'width',8,'height',4,'color','rgb')
close(f)

% --- Instaneous amp/freq
f=figure('Position', [10, 50, 900, 600], 'visible','off');
clf;
subplot(221)
plot(hplus.instamp.Time, hplus.instamp.Data)

subplot(222)
plot(hplus.instfreq.Time, hplus.instfreq.Data)

subplot(223)
plot(hcross.instamp.Time, hcross.instamp.Data)

subplot(224)
plot(hcross.instfreq.Time, hcross.instfreq.Data)

saveas(f, [outdir,'/hilbert'], 'm')
exportfig([outdir,'/hilbert.png'],'width',8,'height',4,'color','rgb')
close(f)

% --------------------------------------
%% Populate summary object
% This is just a structure for now, but we could define a useful class for
% this.

summary.sngl_snr = data.snrDet;
summary.net_snr = data.snrComb;
summary.hplus = hplus;
summary.hcross = hcross;

%% Subroutines

    % Get the instantaneous amplitude / frequency with a hilbert transform
    function [amp,freq] = hilbert_content(data, deltaT)
        % Hilbert transform
        h=hilbert(data);
        unrolled_phase = unwrap(angle(h));
        freq = timeseries(squeeze(diff(unrolled_phase)/deltaT/(2*pi)));
        amp = timeseries(squeeze(abs(h)));
    end

end


