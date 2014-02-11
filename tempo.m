%%%%%%%%%%%%%%
% MUMT307 Final Project
% Ben Deverett
% McGill University
% April 17, 2012

% fd = frequency domain
% td = time domain
%%%%%%%%%%%%%%

%%%%%%%%%%%% Parameters %%%%%%%%%%%%
filename = '128.wav'
band_separation = 'yes'
halfhan_win_length = 0.2 %seconds
comb_n_impulses = 3
tempo_min = 60 %bpm
tempo_max = 160 %bpm
tempo_resolution = 1 %bpm

%%%%%%%%%%%% Read file %%%%%%%%%%%%
[signal, fs, nbits] = wavread(filename);
signal = signal(:,1); %mono

%%%%%%%%%%%% Frequency Band Separation %%%%%%%%%%%%
fd_signal = fft(signal);
if strcmp(band_separation,'yes')
    bands_f = [0, 200; 200, 400; 400, 800; 800, 1600; 1600, 3200; 3200, fs/2];
elseif strcmp(band_separation,'no')
    bands_f = [0,fs/2];
end
fd_sig_bands = [];
for i=1:numel(bands_f(:,1))
    lower_lim = bands_f(i,1) * length(signal) / fs + 1;
    upper_lim = bands_f(i,2) * length(signal) / fs;
    fd_sig_band = zeros(length(fd_signal),1);
    fd_sig_band(lower_lim : upper_lim) = fd_signal(lower_lim : upper_lim);
    fd_sig_bands{i} = fd_sig_band;
end

%%%%%%%%%%%% Envelope Extraction %%%%%%%%%%%%
han_window = zeros(length(signal),1);
han_window_len = halfhan_win_length * 2; %seconds
han_window_len = han_window_len * fs; %samples
for i = 1:han_window_len/2
    han_window(i) = 0.5 - 0.5*cos(2*pi*i / han_window_len);
end
han_window_fft = fft(han_window);
td_sig = [];
for i = 1:numel(fd_sig_bands)
    td_sig_band = abs(ifft(fd_sig_bands{i})); %rectify (full wave)
    fd_sig_band = fft(td_sig_band);
    convolved = fd_sig_band.*han_window_fft(1:length(fd_sig_band));
    td_sig{i} = ifft(convolved);
end

%%%%%%%%%%%% Signal Differentiation and Rectification %%%%%%%%%%%%
for i = 1:numel(td_sig)
   sig = td_sig{i};
   diff_sig = zeros(numel(sig),1);
   %Differentiate:
   for j = 2:numel(sig)
      derivative = sig(j)-sig(j-1);
      diff_sig(j) = derivative;
   end
   %Half wave rectify: 
   for j = 1:numel(diff_sig)
     if diff_sig(j) < 0
        diff_sig(j) = 0;
     end
   end
   diff_sigs{i} = diff_sig;
end
% To Frequency Domain:
fd_sig = [];
for i=1:numel(diff_sigs)
   fd_sig{i} = fft(diff_sigs{i}); 
end

%%%%%%%%%%%% Resonators / Comb Filters %%%%%%%%%%%%
max_energy = 0;
max_tempo = 0;

for tempo = tempo_min:tempo_resolution:tempo_max
   energy = 0;
   samples_per_impulse = floor(60 / tempo * fs);
   for i = 1:numel(fd_sig)
       band = fd_sig{i};
       comb = zeros(length(band),1);
       for c = 1:comb_n_impulses %make impulse train
          comb(c*samples_per_impulse) = 1;
       end
       fd_comb = fft(comb);
       convolved = abs(fd_comb.*band);
       convolved_energy = sum(convolved.^2);
       energy = energy + convolved_energy;
   end
   if energy > max_energy
      max_energy = energy;
      max_tempo = tempo;
   end
end

max_tempo