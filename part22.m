% Salt Bridge Artifact Detection 
clear all; close all; clc;

data = load('partIII_Group25.txt');
Fs = 125; 
chanels = {'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'};
num_electrodes = length(chanels);
L = length(data);
t = (0:L-1)/Fs;

ED_norm = zeros(num_electrodes-1, 1);
titles = cell(num_electrodes-1, 1);

for i = 1:(num_electrodes - 1)
    pair_name = [chanels{i}, ' VS ', chanels{i+1}];    
    titles{i} = pair_name;
    
  
    Pij = data(:,i) - data(:,i+1);
    
    % var(Pij) = (1/T) * sum((Pij - mean(Pij)).^2)
    raw_ED = var(Pij);
    
    mean_var = mean([var(data(:,i)), var(data(:,i+1))]);
    ED_norm(i) = raw_ED / mean_var;
 
    Pij_centered = Pij - mean(Pij);      
    Pij_abs = abs(Pij_centered);   
    
    figure();
    plot(t, data(:,i), 'b'); hold on;
    plot(t, data(:,i+1), 'r');
    plot(t, Pij_abs, 'k', 'LineWidth', 1.2);     
    
    title(pair_name);
    legend(chanels{i}, chanels{i+1}, '|Difference|');
    xlabel('Time [sec]'); ylabel('Amplitude [muV]');
    grid on; xlim([10 15]); 
end

DetectionTable = table(ED_norm, 'RowNames', titles, 'VariableNames', {'Normalized_ED'});
disp(DetectionTable);

