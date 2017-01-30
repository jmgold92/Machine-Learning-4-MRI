clear; clc;

% Import data from csv
cest_tumor = dlmread('cest_tumor.csv');
offset_tumor = dlmread('offset_tumor.csv');
cest_muscle = dlmread('cest_muscle.csv');
offset_muscle = dlmread('offset_muscle.csv');

% Preallocate
cest_fit_tumor = zeros(size(cest_tumor,1),201);
cest_fit_muscle = zeros(size(cest_tumor,1),201);
cest_rsq_tumor = zeros(size(cest_tumor,1),1);
cest_rsq_muscle = zeros(size(cest_tumor,1),1);

% Prescribe number of pools for lorentzian fit
method.Npools=3;

% Prescribe initial guesses and lower and upper bounds for amplitudes,
% widths, and offsets of lorentzian peaks from the 3 pools
method.x0=[0, 0.2, 0.3; 0.4, 1, 1; 0, 0.1, 0.3; 0, 2, 6; 0.5, 4, 6; 0, 2, 6; -6, -4, -2; -1.5, 0, 1.5; 1.5, 3, 4];

% Choose some space over which to display lorentzian fit
xpred=linspace(-8,8,201);

for j=1:size(cest_tumor,1)
    % TOI
    zspec=cest_tumor(j,:);
    ppm=offset_tumor(j,:);
    if isnan(zspec(end)) % Check size of z spectrum (some are 49, some are 101)
        zspec=zspec(1:49);
        ppm=ppm(1:49);
    else
        zspec=zspec(2:end); % If z spectrum has 101 points, the first point is from achieving steady state and is removed
        ppm=ppm(2:end);
    end
    method.range=[1,length(zspec)];
    
    fit=cf_Lorentzian(zspec,ppm,method,zspec(1)); % Lorentzian fitting
    % fit.pars contains amplitude, width, and offset information which can
    % be optionally stored
    fit.pars(7:9)=fit.pars(7:9)-fit.pars(8);
    cest_rsq_tumor(j)=fit.rsq;
    cest_fit_tumor(j,:)=lorentzian(fit.pars,xpred)';
    
    % Muscle
    zspec=cest_muscle(j,:);
    ppm=offset_muscle(j,:);
    if isnan(zspec(end)) % Check size of z spectrum (some are 49, some are 101)
        zspec=zspec(1:49);
        ppm=ppm(1:49);
    else
        zspec=zspec(2:end); % If z spectrum has 101 points, the first point is from achieving steady state and is removed
        ppm=ppm(2:end);
    end
    method.range=[1,length(zspec)];
    
    fit=cf_Lorentzian(zspec,ppm,method,zspec(1)); % Lorentzian fitting
    % fit.pars contains amplitude, width, and offset information which can
    % be optionally stored
    fit.pars(7:9)=fit.pars(7:9)-fit.pars(8);
    cest_rsq_muscle(j)=fit.rsq;
    cest_fit_muscle(j,:)=lorentzian(fit.pars,xpred)';
end

%% Classification

% Optionally screen data for bad fits (r-squared coefficient)
i_tumor = cest_rsq_tumor>0.9;
i_muscle = cest_rsq_muscle>0.9;

% Load class information
cest_classes = readtable('cest_classes.csv');

% Prepare data for classification
cest_table_tumor = array2table(cest_fit_tumor(i_tumor,:));
cest_table_muscle = array2table(cest_fit_muscle(i_muscle,:));
cest_table_tumor.response = cest_classes.Cell_Lines(i_tumor); % Response here was chosen to be the cancer cell line
cest_table_muscle.response = cest_classes.Cell_Lines(i_muscle);

% Pass data to a classifier or use the Classification Learner App
[metrics_tumor_PCA,CM_tumor_PCA]=trainLDA_PCA(cest_table_tumor(:,1:51),cest_table_tumor.response,3,size(cest_table_tumor,1),'Hs766T');
[metrics_tumor_noPCA,CM_tumor_noPCA]=trainLDA(cest_table_tumor(:,1:51),cest_table_tumor.response,size(cest_table_tumor,1),'Hs766T');
    