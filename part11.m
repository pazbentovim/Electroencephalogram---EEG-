%all the data
close all; clear; clc;
data = readmatrix('sharon_2.csv'); 
col = 5:11;
mes = data(:, col); 
elect = {'FC5', 'O1', 'O2', 'P8', 'FC6', 'F4', 'F8'};
Fs = 125; %[Hz]
Ts = 1/Fs; %[sec]

% FIR Bandpass (1-40 Hz) 
b = designfilt('bandpassfir', 'FilterOrder', 100, ...
               'CutoffFrequency1', 1, 'CutoffFrequency2', 40, ...
               'SampleRate', Fs);

for i=1:size(mes,2)
    current_raw = mes(:,i) - mean(mes(:,i)); 
    current_filt = filtfilt(b, current_raw);
    
    L = length(current_raw);
    Tmax = (0:L-1)*Ts;
        f = Fs*(0:floor(L/2))/L; 
    
    fft_raw = abs(fft(current_raw))/L;
    fft_filt = abs(fft(current_filt))/L;
    
    data_to_plot_raw = 20*log10(fft_raw(1:length(f)));
    data_to_plot_filt = 20*log10(fft_filt(1:length(f)));
    
    %plot
    figure('Name', ['EEG Analysis: ' elect{i}], 'Color', 'w');
    
    %time
    subplot(2,2,1);
    plot(Tmax, current_raw, 'b');
    title(['Time: ', elect{i}]);xlabel('time [sec]'); ylabel('Amp [muV]');
    grid on; xlim([40 85]);
    
    subplot(2,2,2);
    plot(Tmax, current_filt, 'b');
    title('Filtered Time (1-40 Hz)');xlabel('time [sec]'); ylabel('Amp [muV]');
    grid on; xlim([40 85]);
    
    %frequency
    subplot(2,2,3);
    plot(f, data_to_plot_raw, 'r'); 
    title('Raw FFT Spectrum');
    xlabel('f [Hz]'); ylabel('H [dB]'); xlim([0 60]); grid on;
    
    subplot(2,2,4);
    plot(f, data_to_plot_filt, 'r');
    title('Filtered FFT Spectrum');
    xlabel('f [Hz]'); ylabel('H [dB]'); xlim([0 60]); grid on;
end
%%
%calculationg the waves
bands = {
    'Delta', [0.5, 4];
    'Theta', [4, 8];
    'Alpha', [8, 13];
    'Beta',  [13, 30];
};

start_sample = 40 * Fs + 1; %without calibration
mes_cut = mes(start_sample:end, :); 

L_cut = size(mes_cut, 1);
t = (start_sample-1)*Ts + (0:L_cut-1)' * Ts;

for i = 1:size(mes_cut, 2)
    sig = mes_cut(:, i) - mean(mes_cut(:, i)); 
    
    figure('Name', ['Brain Waves Study: ' elect{i}], 'Color', 'w');
    
    for b = 1:size(bands, 1)
        bandName = bands{b, 1};
        freqRange = bands{b, 2};
        
        fir_coeff = fir1(100, freqRange/(Fs/2), 'bandpass');
        wave = filtfilt(fir_coeff, 1, sig);
        
        subplot(4, 1, b);
        plot(t, wave, 'LineWidth', 0.5);
        title([bandName, ' Wave (', elect{i}, ')']);
        ylabel('A [muV]');
        grid on;
        
        xlim([40 70]); 
    end
    xlabel('Time [sec]');
end
%%
%plot waves
%FC5 delta wave
fc5_raw_full = data(:, 5);
fc5_raw_full = fc5_raw_full(~isnan(fc5_raw_full)); 

fc5_raw = fc5_raw_full(start_sample:end) - mean(fc5_raw_full(start_sample:end)); 
L = length(fc5_raw);
f_axis = Fs*(0:floor(L/2))/L;

%delta filter [0.5 4]
b_delta = fir1(200, [0.5, 4]/(Fs/2), 'bandpass');
fc5_delta = filtfilt(b_delta, 1, fc5_raw);

fft_raw = abs(fft(fc5_raw))/L;
fft_filt = abs(fft(fc5_delta))/L;

figure('Name', 'Electrode FC5: Delta Wave Analysis', 'Color', 'w');
sgtitle('FC5 - Delta wave')
%time
subplot(2,2,1);
plot(t, fc5_raw, 'Color', 'b');
title('Time Domain');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]); 

%time filtered
subplot(2,2,2);
plot(t, fc5_delta, 'b', 'LineWidth', 1);
title('Time Filtered');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]);

%FFT
subplot(2,2,3);
plot(f_axis, 20*log10(fft_raw(1:length(f_axis))), 'Color', 'r');
title('Frequency');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 20]); grid on;

%FFT filtered
subplot(2,2,4);
plot(f_axis, 20*log10(fft_filt(1:length(f_axis))), 'r', 'LineWidth', 1.2);
title('Frequency filtered');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 20]); grid on;

%STD
samples_per_seg = 15 * Fs;
seg1 = fc5_delta(1 : samples_per_seg);
seg2 = fc5_delta(samples_per_seg + 1 : 2 * samples_per_seg);
seg3 = fc5_delta(2 * samples_per_seg + 1 : 3 * samples_per_seg);

std1 = std(seg1);
std2 = std(seg2);
std3 = std(seg3);

fprintf('Delta STD (Open 1): %.4f\n', std1);
fprintf('Delta STD (Closed): %.4f\n', std2);
fprintf('Delta STD (Open 2): %.4f\n', std3);

all_f_min_delta = zeros(1,3);
all_f_max_delta = zeros(1,3);
all_f_avg_delta = zeros(1,3); 

segs_delta = {seg1, seg2, seg3};
names_delta = {'Open 1', 'Closed', 'Open 2'};

fprintf('\n--- Delta Range Analysis (Electrode FC5) ---\n');
for i = 1:3
    current_seg = segs_delta{i};
    L_seg = length(current_seg);
    f_axis_seg = Fs*(0:floor(L_seg/2))/L_seg;
    
    fft_res = abs(fft(current_seg))/L_seg;
    psd = fft_res(1:length(f_axis_seg)).^2; 
    
    [max_val, ~] = max(psd);
    
    threshold = 0.5 * max_val;
    
    active_indices = find(psd >= threshold);
    
    current_f_avg = 0; 
    
    if ~isempty(active_indices)
        all_f_min_delta(i) = f_axis_seg(active_indices(1));
        all_f_max_delta(i) = f_axis_seg(active_indices(end));
        
        current_f_avg = mean(f_axis_seg(active_indices));
        
        all_f_avg_delta(i) = current_f_avg; 
    end
    
    fprintf('Active Delta Range (%s): %.2f Hz - %.2f Hz | Avg Freq: %.2f Hz\n', names_delta{i}, all_f_min_delta(i), all_f_max_delta(i), current_f_avg);
end
    
avg_f_min_delta = mean(all_f_min_delta);
avg_f_max_delta= mean(all_f_max_delta);
overall_avg_f_delta = mean(all_f_avg_delta); 

fprintf('Average Active Delta Range: %.2f Hz - %.2f Hz\n', avg_f_min_delta, avg_f_max_delta);
fprintf('Overall Average Delta Frequency (Avg of Avg Freqs): %.2f Hz\n', overall_avg_f_delta); 

%%
%F8 theta wave
f8_raw_full = data(:, 11);
f8_raw_full = f8_raw_full(~isnan(f8_raw_full)); 
f8_raw = f8_raw_full(start_sample:end) - mean(f8_raw_full(start_sample:end));
L = length(f8_raw);
t = (start_sample-1)*Ts + (0:L-1)' * Ts;
f_axis = Fs*(0:floor(L/2))/L;

% bpf 4-8 Hz
b_theta = designfilt('bandpassfir', 'FilterOrder', 200, ...
               'CutoffFrequency1', 4, 'CutoffFrequency2', 8, ...
               'SampleRate', Fs);

f8_theta = filtfilt(b_theta, f8_raw);
fft_raw = abs(fft(f8_raw))/L;
fft_filt = abs(fft(f8_theta))/L;
figure('Name', 'Electrode F8: Theta Wave Analysis', 'Color', 'w');

sgtitle('F8 - Theta Wave', 'FontSize', 14);

subplot(2,2,1);
plot(t, f8_raw, 'Color', 'b');
title('Time Domain');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]); 

subplot(2,2,2);
plot(t, f8_theta, 'b');
title('Time filtered');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]);

subplot(2,2,3);
plot(f_axis, 20*log10(fft_raw(1:length(f_axis))), 'Color','r');
title('Frequency');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 30]); grid on;

subplot(2,2,4);
plot(f_axis, 20*log10(fft_filt(1:length(f_axis))), 'r', 'LineWidth', 1.2);
title('Frequency filtered');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 30]); grid on;

%STD
samples_per_seg = 15 * Fs;
seg1 = f8_theta(1 : samples_per_seg);
seg2 = f8_theta(samples_per_seg + 1 : 2 * samples_per_seg);
seg3 = f8_theta(2 * samples_per_seg + 1 : 3 * samples_per_seg);

std1 = std(seg1);
std2 = std(seg2);
std3 = std(seg3);

fprintf('Theta STD (Open 1): %.4f\n', std1);
fprintf('Theta STD (Closed): %.4f\n', std2);
fprintf('Theta STD (Open 2): %.4f\n', std3);

all_f_min_theta = zeros(1,3);
all_f_max_theta = zeros(1,3);
all_f_avg_theta = zeros(1,3); 

segs_theta = {seg1, seg2, seg3};
names_theta = {'Open 1', 'Closed', 'Open 2'};

fprintf('\n--- Theta Range Analysis (Electrode F8) ---\n');
for i = 1:3
    current_seg = segs_theta{i};
    L_seg = length(current_seg);
    f_axis_seg = Fs*(0:floor(L_seg/2))/L_seg;
    
    fft_res = abs(fft(current_seg))/L_seg;
    psd = fft_res(1:length(f_axis_seg)).^2; 
    
    [max_val, ~] = max(psd);
    
    threshold = 0.5 * max_val;
    
    active_indices = find(psd >= threshold);
    current_f_avg = 0; 
    
    if ~isempty(active_indices)
        all_f_min_theta(i) = f_axis_seg(active_indices(1));
        all_f_max_theta(i) = f_axis_seg(active_indices(end));
        
        current_f_avg = mean(f_axis_seg(active_indices));
        
        all_f_avg_theta(i) = current_f_avg; 
    end
    
    fprintf('Active Theta Range (%s): %.2f Hz - %.2f Hz | Avg Freq: %.2f Hz\n', names_theta{i}, all_f_min_theta(i), all_f_max_theta(i), current_f_avg);
end

avg_f_min_theta = mean(all_f_min_theta);
avg_f_max_theta = mean(all_f_max_theta);
overall_avg_f_theta = mean(all_f_avg_theta); 

fprintf('Average Active Theta Range (F8): %.2f Hz - %.2f Hz\n', avg_f_min_theta, avg_f_max_theta);
fprintf('Overall Average Theta Frequency (Avg of Avg Freqs): %.2f Hz\n', overall_avg_f_theta); 

%%
%O1 alpha wave
o1_raw_full = data(:, 6);
o1_raw_full = o1_raw_full(~isnan(o1_raw_full)); 
o1_raw = o1_raw_full(start_sample:end) - mean(o1_raw_full(start_sample:end));

L = length(o1_raw);
t = (start_sample-1)*Ts + (0:L-1)' * Ts;
f_axis = Fs*(0:floor(L/2))/L;

b_alpha = designfilt('bandpassfir', 'FilterOrder', 200, ...
               'CutoffFrequency1', 8, 'CutoffFrequency2', 13, ...
               'SampleRate', Fs);

o1_alpha = filtfilt(b_alpha, o1_raw);

fft_raw = abs(fft(o1_raw))/L;
fft_filt = abs(fft(o1_alpha))/L;

figure('Name', 'Electrode O1: Alpha Wave Analysis', 'Color', 'w');

sgtitle('O1 - Alpha Wave', 'FontSize', 14);

subplot(2,2,1);
plot(t, o1_raw, 'Color', 'b');
title('Time Domain');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]); 

subplot(2,2,2);
plot(t, o1_alpha, 'b');
title('Time filtered');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]);

subplot(2,2,3);
plot(f_axis, 20*log10(fft_raw(1:length(f_axis))), 'Color', 'r');
title('Frequency');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 40]); grid on;

subplot(2,2,4);
plot(f_axis, 20*log10(fft_filt(1:length(f_axis))), 'r');
title('Alpha Wave Spectrum');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 40]); grid on;

%STD
samples_per_seg = 15 * Fs;
seg1 = o1_alpha(1 : samples_per_seg);
seg2 = o1_alpha(samples_per_seg + 1 : 2 * samples_per_seg);
seg3 = o1_alpha(2 * samples_per_seg + 1 : 3 * samples_per_seg);

std1 = std(seg1);
std2 = std(seg2);
std3 = std(seg3);

fprintf('Alpha STD (Open 1): %.4f\n', std1);
fprintf('Alpha STD (Closed): %.4f\n', std2);
fprintf('Alpha STD (Open 2): %.4f\n', std3);

all_f_min_alpha = zeros(1,3);
all_f_max_alpha = zeros(1,3);
all_f_avg_alpha = zeros(1,3); 

segs_alpha = {seg1, seg2, seg3};
names_alpha = {'Open 1', 'Closed', 'Open 2'};

fprintf('\n--- Alpha Range Analysis (Electrode O1) ---\n');
for i = 1:3
    current_seg = segs_alpha{i};
    L_seg = length(current_seg);
    Fs = 128; 
    f_axis_seg = Fs*(0:floor(L_seg/2))/L_seg;
    
    fft_res = abs(fft(current_seg))/L_seg;
    psd = fft_res(1:length(f_axis_seg)).^2; 
    
    [max_val, ~] = max(psd);
    
    threshold = 0.5 * max_val;
    
    active_indices = find(psd >= threshold);
    
    current_f_avg = 0; 
    
    if ~isempty(active_indices)
        all_f_min_alpha(i) = f_axis_seg(active_indices(1));
        all_f_max_alpha(i) = f_axis_seg(active_indices(end));
        
        current_f_avg = mean(f_axis_seg(active_indices));
        
        all_f_avg_alpha(i) = current_f_avg; 
    end
    
    fprintf('Active Alpha Range (%s): %.2f Hz - %.2f Hz | Avg Freq: %.2f Hz\n', names_alpha{i}, all_f_min_alpha(i), all_f_max_alpha(i), current_f_avg);
end

avg_f_min_alpha = mean(all_f_min_alpha);
avg_f_max_alpha = mean(all_f_max_alpha);
overall_avg_f_alpha = mean(all_f_avg_alpha); 

fprintf('Average Active Alpha Range: %.2f Hz - %.2f Hz\n', avg_f_min_alpha, avg_f_max_alpha);
fprintf('Overall Average Alpha Frequency (Avg of Avg Freqs): %.2f Hz\n', overall_avg_f_alpha); 
%%
% F4 beta wave
f4_raw_full = data(:, 10); 
f4_raw_full = f4_raw_full(~isnan(f4_raw_full)); 
f4_raw = f4_raw_full(start_sample:end) - mean(f4_raw_full(start_sample:end));

L = length(f4_raw);
t = (start_sample-1)*Ts + (0:L-1)' * Ts;
f_axis = Fs*(0:floor(L/2))/L;

b_beta = designfilt('bandpassfir', 'FilterOrder', 200, ...
               'CutoffFrequency1', 13, 'CutoffFrequency2', 30, ...
               'SampleRate', Fs);

f4_beta = filtfilt(b_beta, f4_raw);
fft_raw = abs(fft(f4_raw))/L;
fft_filt = abs(fft(f4_beta))/L;

figure('Name', 'Electrode F4: Beta Wave Analysis', 'Color', 'w');

sgtitle(' F4 - Beta Wave ', 'FontSize', 14);

subplot(2,2,1);
plot(t, f4_raw, 'Color', 'b');
title('Time Domain');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]); 

subplot(2,2,2);
plot(t, f4_beta, 'b'); 
title('Time filtered');
xlabel('Time [sec]'); ylabel('Amplitude [muV]'); grid on;
xlim([40 85]);

subplot(2,2,3);
plot(f_axis, 20*log10(fft_raw(1:length(f_axis))), 'Color', 'r');
title('Frequency');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 50]); grid on;

subplot(2,2,4);
plot(f_axis, 20*log10(fft_filt(1:length(f_axis))), 'r');
title('Frequency filtered');
xlabel('Freq [Hz]'); ylabel('Magnitude [dB]'); xlim([0.3 50]); grid on;

%STD
samples_per_seg = 15 * Fs;
seg1 = f4_beta(1 : samples_per_seg);
seg2 = f4_beta(samples_per_seg + 1 : 2 * samples_per_seg);
seg3 = f4_beta(2 * samples_per_seg + 1 : 3 * samples_per_seg);

std1 = std(seg1);
std2 = std(seg2);
std3 = std(seg3);

fprintf('Beta STD (Open 1): %.4f\n', std1);
fprintf('Beta STD (Closed): %.4f\n', std2);
fprintf('Beta STD (Open 2): %.4f\n', std3);

all_f_min_beta = zeros(1,3); 
all_f_max_beta = zeros(1,3);
all_f_avg_beta = zeros(1,3); 

segs = {seg1, seg2, seg3};
names = {'Open 1', 'Closed', 'Open 2'};

for i = 1:3
    current_seg = segs{i};
    L_seg = length(current_seg);
    f_axis_seg = Fs*(0:floor(L_seg/2))/L_seg;
    
    fft_res = abs(fft(current_seg))/L_seg;
    psd = fft_res(1:length(f_axis_seg)).^2; 
    
    [max_val, ~] = max(psd);
    
    threshold = 0.5 * max_val;
    
    active_indices = find(psd >= threshold);
    
    current_f_avg = 0; 
    
    if ~isempty(active_indices)
        all_f_min_beta(i) = f_axis_seg(active_indices(1));
        all_f_max_beta(i) = f_axis_seg(active_indices(end));
        
        current_f_avg = mean(f_axis_seg(active_indices));
        
        all_f_avg_beta(i) = current_f_avg; 
    end
    
    fprintf('Active Beta Range (%s): %.2f Hz - %.2f Hz | Avg Freq: %.2f Hz\n', names{i}, all_f_min_beta(i), all_f_max_beta(i), current_f_avg);
end

avg_f_min = mean(all_f_min_beta);
avg_f_max = mean(all_f_max_beta);
overall_avg_f_beta = mean(all_f_avg_beta); 

fprintf('Average Active Beta Range: %.2f Hz - %.2f Hz\n', avg_f_min, avg_f_max);
fprintf('Overall Average Frequency (Avg of Avg Freqs): %.2f Hz\n', overall_avg_f_beta);
