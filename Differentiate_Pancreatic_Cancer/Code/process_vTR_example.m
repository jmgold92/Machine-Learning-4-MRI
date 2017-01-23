clear; clc;

% Import data from csv
vtr_tumor = dlmread('vtr_tumor.csv');
tr_tumor = dlmread('tr_tumor.csv');
vtr_muscle = dlmread('vtr_muscle.csv');
tr_muscle = dlmread('tr_muscle.csv');

% Preallocate
T1_tumor = zeros(size(vtr_tumor,1),1);
T1_rsq_tumor = T1_tumor;
T1_muscle = zeros(size(vtr_muscle,1),1);
T1_rsq_muscle = T1_muscle;

% Fitting Method
options = optimoptions('lsqcurvefit');  
% Max. number of iterations
options.MaxIter=100E3; 
% Max. number of func. evaluations
options.MaxFunEvals=200;                
% Tolerance for NLSQ method
options.TolX=1e-4;             options.TolFun=1e-4;
% Turn iterative display off
options.Display='off';
options.Algorithm='trust-region-reflective';

for j=1:length(T1_tumor)
    s = vtr_tumor(j,:)./vtr_tumor(j,end);
    tr = tr_tumor(j,:);
    if tr(1)==1.2101
        s = [mean(s(1:7)),s(8:end)];
        tr = [tr(1),tr(8:end)];
    end
    pars=lsqcurvefit(@vTR,[2.5,1.0],tr,s,[0,0],[10,10],options);
    T1_tumor(j)=pars(2);
    T1_rsq_tumor(j)=rsq(s,vTR(pars,tr));
    
    s = vtr_muscle(j,:)./vtr_muscle(j,end);
    tr = tr_muscle(j,:);
    if tr(1)==1.2101
        s = [mean(s(1:7)),s(8:end)];
        tr = [tr(1),tr(8:end)];
    end
    pars=lsqcurvefit(@vTR,[2.5,1.0],tr,s,[0,0],[10,10],options);
    T1_muscle(j)=pars(2);
    T1_rsq_muscle(j)=rsq(s,vTR(pars,tr));
end

%% Classification

% Optionally screen data for bad fits (r-squared coefficient)
i_tumor = T1_rsq_tumor>0.9;
i_muscle = T1_rsq_muscle>0.9;

% Load class information
T1_classes = readtable('vtr_classes.csv');

% Prepare data for classification
T1_table_tumor = array2table(T1_tumor(i_tumor,:));
T1_table_muscle = array2table(T1_muscle(i_muscle,:));
T1_table_tumor.response = T1_classes.Cell_Lines(i_tumor); % Response here was chosen to be the cancer cell line
T1_table_muscle.response = T1_classes.Cell_Lines(i_muscle);

% Pass data to a classifier or use the Classification Learner App
[metrics_tumor,CM_tumor]=trainLDA(T1_table_tumor(:,1),T1_table_tumor.response,size(T1_table_tumor,1),'Hs766T');
[metrics_muscle,CM_muscle]=trainLDA(T1_table_muscle(:,1),T1_table_muscle.response,size(T1_table_muscle,1),'Hs766T');