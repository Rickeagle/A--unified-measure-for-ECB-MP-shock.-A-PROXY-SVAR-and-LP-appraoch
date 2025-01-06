clc;
clear;
close all;

%% Load and clean the first dataset
data = readtable('z.brw-shock-series (1).xlsx');
cleanedBRWshock = data(~all(ismissing(data), 2), :);
cleanedBRWshock(1:9, :) = [];
cleanedBRWshock = cleanedBRWshock(~all(ismissing(cleanedBRWshock), 2), :);
writetable(cleanedBRWshock, 'cleaned_BRWshock.xlsx');
% Load and process the second dataset
ldata = readtable('Monthly_MPshocks.xlsx');
BRW = str2double(string(cleanedBRWshock{:, 2})); % Ensures conversion to numeric
MP_EuroArea = str2double(string(ldata{1:180, 2:6}));

%% Compute and display initial correlation
correlation = corr(BRW, MP_EuroArea(:, 1), 'Rows', 'pairwise');
correlation30 = corr(BRW, MP_EuroArea(:,4), 'Rows', 'pairwise');
fprintf('correlation between BRW and MP_EuroArea2: %.4f\n', correlation);
fprintf('correlation between BRW and MP_EuroArea30: %.4f\n', correlation30);

%% Remove rows where shocks are zero
validIdx = BRW ~= 0 & all(MP_EuroArea ~= 0, 2);
BRW_cleaned = BRW(validIdx);
MP_EuroArea_cleaned = MP_EuroArea(validIdx, :);

%% Compute and display correlations for cleaned data
labels = {'MP_EuroArea2', 'MP_EuroArea1', 'MP_EuroArea10', 'MP_EuroArea30', 'MP_EuroArea3mth'};
for i = 1:size(MP_EuroArea_cleaned, 2)
    correlation = corr(BRW_cleaned, MP_EuroArea_cleaned(:, i));
    fprintf('Correlation between BRW and %s when shocks !=0: %.4f\n', labels{i}, correlation);
end
