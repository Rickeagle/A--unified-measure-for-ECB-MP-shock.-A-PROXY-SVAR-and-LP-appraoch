clear
clc
close all

%% 1) Load Data
ldata=xlsread('MacroData.xlsx');
% Loan Rate, Infl, Unemp rate, IP, House Price rate, Labor cost index
Y=[ldata(:,[1 2 3 4 5 6] )];
var_names={'Loan Rate','INFL','UR','IP','House Price','Labor cost'};
d=size(Y,2); 
%% 2) Parameters
n=size(Y,2); 
H=48; % horizon of the IRF
c=1;
MaxBoot=200;
cl=0.9;
ind=[];
opt=1;
signif_level=0.32;

% Lag selection
maxLag = 12; 
criteria = zeros(maxLag, 3); % Columns for AIC, BIC, HQIC
for lag = 1:maxLag
    % Estimate VAR model with lag 'lag'
    [~, ~, ~, logL] = estimate_VAR(Y, lag, c);
    % Compute the criteria
    T = size(Y, 1) - lag; % Effective sample size after lags
    k = n^2 * lag + n; 
    criteria(lag, 1) = -2 * logL / T + 2 * k / T;       % AIC
    criteria(lag, 2) = -2 * logL / T + log(T) * k / T;  % BIC
    criteria(lag, 3) = -2 * logL / T + 2 * log(log(T)) * k / T; % HQIC
end

% Select the optimal lag based on minimum BIC (or choose AIC/HQIC)
[~, p_bic] = min(criteria(:, 2)); 
disp(['Optimal lag (BIC): ', num2str(p_bic)]);
[~, p_hqic] = min(criteria(:, 3)); 
disp(['Optimal lag (HQIC): ', num2str(p_hqic)]);

[~, p_aic] = min(criteria(:, 1));
disp(['Optimal lag (AIC): ', num2str(p_aic)]);
p=p_aic;
%using AIC optimal p since we are working in small samples
%addpath('C:\Users\Utente\Desktop\A-UNITO\Second Year\Monetary\Gambetti\MATLAB\USED FILE')
[wirf,res] = WoldBoot(Y,p,H,c);

%% Proxy-SVAR SVAR-IV approach
Z=ldata(:,7);
Z=Z(p+1:end);
Z=Z(~isnan(Z));


% Regression on epsilon onto the instrument
b=(inv(Z'*Z)*Z'*res)';
b=b/b(4);

for j=1:H
    mpirf(:,j)=wirf(:,:,j)*b;
end

%[irfs,irfsboot] = PROXYSVARRecoverabilityCase(Y,Z,p,H,1,ind,MaxBoot,4,1,0,Y);

%% Plot
figure(1);
plot(mpirf', 'LineWidth', 1.5); % Plot all impulse responses with thicker lines
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('Impulse Response Functions (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);
saveas(figure(1),'Proxy-SVAR.png')


%% ROBUSTNESS Proxy-SVAR SVAR-IV approach using MP10 instead of MP2
% Rembvoe NaN
Z10=ldata(:,9); %FOR MP10
Z10=Z10(p+1:end);
Z10=Z10(~isnan(Z10));

% Regression on epsilon onto the instrument
b10=(inv(Z10'*Z10)*Z10'*res)';
b10=b10/b10(4);

for j=1:H
    mpirf10(:,j)=wirf(:,:,j)*b10;
end
r=1;
ind=[2 3 4 5 6 7];
% Plot for MP10
figure(2);
plot(mpirf10', 'LineWidth', 1.5);
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('Impulse Response Functions with MP10 (IRFs)', 'FontSize', 14);

legend(var_names, 'Location', 'best', 'FontSize', 10);

saveas(figure(2),'Proxy-SVAR_MP10.png')
%Not good at all
%% ROBUSTNESS Proxy-SVAR SVAR-IV approach using MP1 instead of MP2
% Rembvoe NaN
Z1=ldata(:,8); %FOR MP1
Z1=Z1(p+1:end);
Z1=Z1(~isnan(Z1));

% Regression on epsilon onto the instrument
b1=(inv(Z1'*Z1)*Z1'*res)';
b1=b1/b1(4);

for j=1:H
    mpirf1(:,j)=wirf(:,:,j)*b1;
end
r=1;
ind=[2 3 4 5 6 7];
% Plot for MP1
figure(3);
plot(mpirf1', 'LineWidth', 1.5);
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('Impulse Response Functions with MP1 (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);

saveas(figure(3),'Proxy-SVAR_MP1.png')

%% ROBUSTNESS Proxy-SVAR SVAR-IV approach using MP30 instead of MP2
Z30=ldata(:,10); %FOR MP30
Z30=Z30(p+1:end);
Z30=Z30(~isnan(Z30));

% Regression on epsilon onto the instrument
b30=(inv(Z30'*Z30)*Z30'*res)';
b30=b30/b30(4);

for j=1:H
    mpirf30(:,j)=wirf(:,:,j)*b30;
end
r=1;
ind=[2 3 4 5 6 7];
% Plot for MP30
figure(4);
plot(mpirf30', 'LineWidth', 1.5);
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('Impulse Response Functions with MP30 (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);

saveas(figure(4),'Proxy-SVAR_MP30.png')
%% ROBUSTNESS Proxy-SVAR SVAR-IV approach using MP3mth instead of MP2
% Rembvoe NaN
Z3=ldata(:,11); %FOR MP3mth
Z3=Z3(p+1:end);
Z3=Z3(~isnan(Z3));

% Regression on epsilon onto the instrument
b3=(inv(Z3'*Z3)*Z3'*res)';
b3=b3/b3(4);

for j=1:H
    mpirf3(:,j)=wirf(:,:,j)*b3;
end
r=1;
ind=[2 3 4 5 6 7];
% Plot for MP3mth
figure(5);
plot(mpirf3', 'LineWidth', 1.5);
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('Impulse Response Functions with MP3mth (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);

saveas(figure(5),'Proxy-SVAR_MP3mth.png')


%% Plagborg Moller code 

plot_var=3;
plot_xticks      = 0:12:H; 
plot_band_xlabel = 'months after shock';
plot_band_ylabel = 'percent';

numdraws_supt        = 1e5; % Number of normal draws used to compute plug-in sup-t crit. val.

numdraws_boot        = 1e4; % Set to 0 if bootstrap inference undesired

verbose = true;             % Print progress of bootstrap procedure

rng(20170114);              % Seed for random number generator

% Bands

band_list            = {'Pwise',...
                        'supt',...
                        'Sidak',...
                        'Bonferroni'}; %, 'thetaproj', 'muproj'}; 
                            % Cell array of bands to be plotted, 
                            % can be any combination of: 'Pwise', 'supt', 'Sidak', 'Bonferroni', 'thetaproj', 'muproj'

legend_bands         = {'Pointwise',...
                        'Sup-t: plug-in',...
                        'Sidak',...
                        'Bonferroni'}; %, '\theta-projection', '\mu-projection'}; 
                            % Legend for bands

linestyle_supt_bands = {'-', '--'};

%% Reduced-form VAR estimation

redf                 = iv_RedForm(Y, Z, p); 
                                      % Struct redf contains reduced-form VAR
                                      % and IV objects
%% 
for plot_var = 1:d  % Loop through all variables
    
    % IV estimation of IRFs and delta method variance
    [Theta, Sigmahat, p] = iv_estim(redf, H);

    % Confidence bands
    sel = select_IRF(d, 1, H, plot_var, 1); % Selection vector for IRF of interest

    bands_plugin = SimInference.bands_plugin(Theta(sel)',...
                   Sigmahat(sel,sel), p, band_list, numdraws_supt,...
                   1-signif_level); 

    % Plot plug-in bands
    fig1 = figure('Unit', 'normalize', 'Position', [0.2 0.2 0.6 0.6], ...
                  'Name', ['Bands: IV, plug-in, Var ' num2str(plot_var)]);

    plot_compare({Theta(sel)'}, bands_plugin, plot_band_xlabel, plot_band_ylabel, plot_xticks, legend_bands);

    % Save figure 1 with updated name
    saveas(fig1, ['Var' num2str(plot_var) '_Bands_IV_Plug_in.png']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL PROJECTION

% ) Estimation of LP

[Ytemp, Xtemp] = LagY(Y, p); 
L= ldata(:, 7); 
zz = L(p+1:end); 

for i = 1:n
    for j = 1:H
        Y_sub = Ytemp(j:end, i);
        X = [ones(size(Y_sub, 1), 1), zz(1:end-j+1), Xtemp(1:end-j+1, 2:end)];
        b = (X' * X) \ (X' * Y_sub); % Regression coefficients
        irf(j, i) = b(2); 
    end
end
irf = irf / irf(1, 4); % Normalize IRF

% Plot LP
plot (irf)

%% Plot SVAR-IV and LP IRFs Side by Side
figure;

% SVAR-IV Plot
subplot(1, 2, 1); % Create a subplot (1 row, 2 columns, position 1)
plot(mpirf', 'LineWidth', 1.5); % Plot all impulse responses with thicker lines
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('SVAR-IV: Impulse Response Functions (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);

% LP Plot
subplot(1, 2, 2); 
plot(irf, 'LineWidth', 1.5); 
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('LP: Impulse Response Functions (IRFs)', 'FontSize', 14);
legend(var_names, 'Location', 'best', 'FontSize', 10);

% Adjust layout if needed
set(gcf, 'Position', [100, 100, 1200, 500]); 
saveas(gcf, 'ProxySVAR-LP.png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WHAT ABOUT USING THE GARCH(1,1) AS AN UNCERTAINTY MEASURE?
U=ldata(:,14); %FOR Uncertainty
U=U(p+1:end);
U=U(~isnan(U));

% Regression on epsilon onto the instrument
bU=(inv(U'*U)*U'*res)';
bU=bU/bU(4);

for j=1:H
    mpirfU(:,j)=wirf(:,:,j)*bU;
end

% Plot for Uncertainty as an instrument
figure(10);
plot(mpirfU', 'LineWidth', 1.5);
grid on;
xlabel('Horizon (months)', 'FontSize', 12);
ylabel('Impulse Response', 'FontSize', 12);
title('IRF using UNCERTAINTY', 'FontSize', 14);

legend(var_names, 'Location', 'best', 'FontSize', 10);

saveas(figure(10),'Proxy-SVAR_Uncertainty.png')