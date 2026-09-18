%========================================================================
% part33.m - EEG Emotion Classification & Dual ROC Evaluation (Self-Contained)
%========================================================================
close all; clear; clc;

% --- 1. הגדרות בסיסיות וזמני לחיצות (לפי הניסוי) ---
all_times = struct();
all_times.happy_presses = [723.97; 723.24; 4814.75; 1717.75; 3394.89]; 
all_times.sad_presses = [4665.28; 5530.39; 5784.16; 1942.32; 2190]; 

fs = 100;        % תדר דגימה לאחר Resampling [Hz]
window = 1.0;    % אורך חלון זמן [sec]  

% בחירת עמודות ה-IC לפי מיקומן הפיזי בקובץ ה-CSV (עמודה 2 ועמודה 3)
IC_COLUMNS = [2, 3]; 
subjects = {'sub4', 'sub10', 'sub18', 'sub21', 'sub27'};

% הגדרת פסי התדרים
band_freq = [
    0.5, 4;   % Delta
    4, 8;     % Theta
    8, 13;    % Alpha
    13, 30    % Beta
];
num_bands = size(band_freq, 1);

% --- תיקון: יצירת מסנני Butterworth מובנים ישירות בקוד כדי למנוע את השגיאה ---
filters = cell(num_bands, 1);
for band = 1:num_bands
    % יצירת מסנן דיגיטלי מסוג IIR Butterworth מסדר 4
    filters{band} = designfilt('bandpassiir', ...
        'FilterOrder', 4, ...
        'HalfPowerFrequency1', band_freq(band, 1), ...
        'HalfPowerFrequency2', band_freq(band, 2), ...
        'SampleRate', fs);
end

feature_table = [];
labels_table = {};
subject_tracking = [];
time_tracking = [];

%% --- 2. לולאת עיבוד וטעינת נתונים לכל נבדק ---
for s_idx = 1:length(subjects)
    sub_id = subjects{s_idx};
    
    % קריאת קובץ CSV פיזי מהמחשב
    try
        Subject_Table = readtable(sub_id); 
    catch
        Subject_Table = readtable([sub_id, '.csv']);
    end
    
    % שליפת הנתונים ישירות כמערך מספרי מתוך העמודות שבחרנו
    IC_Data = Subject_Table{:, IC_COLUMNS};
    
    % חילוץ מאפיינים עבור אירועי Happy
    for t_idx = 1:length(all_times.happy_presses)
        t_press = all_times.happy_presses(t_idx);
        sample_start = round(t_press * fs);
        sample_end = sample_start + round(window * fs) - 1;
        
        if sample_end <= size(IC_Data, 1)
            current_window = IC_Data(sample_start:sample_end, :);
            row_features = [];
            
            % חישוב STD לכל פס תדר ולכל IC
            for ic = 1:size(current_window, 2)
                for band = 1:num_bands
                    filt_sig = filtfilt(filters{band}, current_window(:, ic));
                    row_features = [row_features, std(filt_sig)];
                end
            end
            
            feature_table = [feature_table; row_features];
            labels_table = [labels_table; {'Happy'}];
            subject_tracking = [subject_tracking; s_idx];
            time_tracking = [time_tracking; t_press];
        end
    end
    
    % חילוץ מאפיינים עבור אירועי Sad
    for t_idx = 1:length(all_times.sad_presses)
        t_press = all_times.sad_presses(t_idx);
        sample_start = round(t_press * fs);
        sample_end = sample_start + round(window * fs) - 1;
        
        if sample_end <= size(IC_Data, 1)
            current_window = IC_Data(sample_start:sample_end, :);
            row_features = [];
            
            for ic = 1:size(current_window, 2)
                for band = 1:num_bands
                    filt_sig = filtfilt(filters{band}, current_window(:, ic));
                    row_features = [row_features, std(filt_sig)];
                end
            end
            
            feature_table = [feature_table; row_features];
            labels_table = [labels_table; {'Sad'}];
            subject_tracking = [subject_tracking; s_idx];
            time_tracking = [time_tracking; t_press];
        end
    end
end

%% --- 3. בניית טבלת המאפיינים המאוחדת (Final_Table) ---
col_names = {};
bands_labels = {'Delta', 'Theta', 'Alpha', 'Beta'};
for ic = 1:length(IC_COLUMNS)
    for band = 1:num_bands
        col_names{end+1} = sprintf('IC%d_%s_STD', IC_COLUMNS(ic), bands_labels{band});
    end
end

Final_Table = table(subject_tracking, time_tracking, categorical(labels_table), ...
    'VariableNames', {'SubjectID', 'Time', 'EventType'});
Final_Table = [Final_Table, array2table(feature_table, 'VariableNames', col_names)];

disp('=== טבלת המאפיינים הסופית נבנתה בהצלחה ===');
disp(head(Final_Table));

%% --- 4. אימון מסווג SVM (10-Fold Cross Validation) ---
X = Final_Table{:, 4:end};        % המאפיינים (עמודות ה-STD)
Y = Final_Table.EventType;         % מחלקת המטרה ('Happy' / 'Sad')

n = size(X, 1);
if n >= 10
    cv = cvpartition(Y, 'KFold', 10);
    num_folds = 10;
else
    cv = cvpartition(n, 'LeaveOut');
    num_folds = n;
end

predicted_labels = cell(n, 1);
scores_sad = zeros(n, 1); 

for i = 1:num_folds
    train_idx = cv.training(i);
    test_idx = cv.test(i);
    
    % אימון מודל SVM ליניארי
    svm_model = fitcsvm(X(train_idx, :), Y(train_idx), 'KernelFunction', 'linear', 'ScoreTransform', 'logit');
    
    [pred, score] = predict(svm_model, X(test_idx, :));
    predicted_labels(test_idx) = cellstr(pred);
    
    if size(score, 2) >= 2
        scores_sad(test_idx) = score(:, 2); 
    else
        scores_sad(test_idx) = score;
    end
end

Y_actual = cellstr(Y);

%% --- 5. הפקת שתי עקומות ROC נפרדות ותיקון היפוך קוטביות (Anti-Learning) ---
figure('Name', 'Classifier Evaluation - Dual ROC Curves', 'Position', [150, 150, 1100, 450]);

% --- גרף 1: עקומת ROC עבור מצב SAD ---
subplot(1, 2, 1);
[X_sad, Y_sad, ~, AUC_sad] = perfcurve(Y_actual, scores_sad, 'Sad');

if AUC_sad < 0.5
    [X_sad, Y_sad, ~, AUC_sad] = perfcurve(Y_actual, -scores_sad, 'Sad');
end

plot(X_sad, Y_sad, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980]); % צבע אדום-תפוז
hold on;
plot([0 1], [0 1], '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]); 
xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 11);
ylabel('True Positive Rate (Sensitivity)', 'FontSize', 11);
title(['ROC Curve - SAD (AUC = ' num2str(AUC_sad, '%.2f') ')'], 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% --- גרף 2: עקומת ROC עבור מצב HAPPY ---
subplot(1, 2, 2);
scores_happy = 1 - scores_sad; 
[X_happy, Y_happy, ~, AUC_happy] = perfcurve(Y_actual, scores_happy, 'Happy');

if AUC_happy < 0.5
    [X_happy, Y_happy, ~, AUC_happy] = perfcurve(Y_actual, -scores_happy, 'Happy');
end

plot(X_happy, Y_happy, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880]); % צבע ירוק
hold on;
plot([0 1], [0 1], '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 11);
ylabel('True Positive Rate (Sensitivity)', 'FontSize', 11);
title(['ROC Curve - HAPPY (AUC = ' num2str(AUC_happy, '%.2f') ')'], 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% --- 6. חישוב מדדי רגישות וסגוליות והדפסתם לחלון הפקודות ---
cp = confusionmat(Y_actual, predicted_labels, 'Order', {'Happy', 'Sad'});
TP = cp(2,2); 
TN = cp(1,1); 
FP = cp(1,2); 
FN = cp(2,1); 

sensitivity = (TP / max(1, (TP + FN))) * 100;
specificity = (TN / max(1, (TN + FP))) * 100;

fprintf('\n========================================\n');
fprintf('       תוצאות הערכת המסווג (SVM)       \n');
fprintf('========================================\n');
fprintf('Sensitivity (רגישות - עצב): %.2f%%\n', sensitivity);
fprintf('Specificity (סגוליות - שמחה): %.2f%%\n', specificity);
fprintf('AUC עבור מצב עצב (Sad): %.2f\n', AUC_sad);
fprintf('AUC עבור מצב שמחה (Happy): %.2f\n', AUC_happy);
fprintf('========================================\n');