clear; clc;

% Import data from csv
dce_tumor = dlmread('dce_tumor.csv');
dce_muscle = dlmread('dce_muscle.csv');

% Import class information
dce_classes = readtable('dce_classes.csv');

% Transform data
dce_tumor = dce_tumor./repmat(mean(dce_tumor(:,1:4),2),1,150)-1; % Scale to baseline
dce_muscle = dce_muscle./repmat(mean(dce_muscle(:,1:4),2),1,150)-1;

% Screen data (here based on percent enhancement)
i_tumor = max(dce_tumor,[],2)-abs(median(dce_tumor(:,1:4),2))>=0.35;
i_tumor = logical(i_tumor-isnan(dce_tumor(:,150)));
i_tumor(14) = 0; % By visual inspection, this data is erroneous
i_muscle = max(dce_muscle,[],2)-abs(median(dce_muscle(:,1:4),2))>=0.32;
i_muscle = logical(i_muscle-isnan(dce_muscle(:,150)));
i_muscle(14) = 0; % By visual inspection, this data is erroneous

%% Find cumulative integrals

int_tumor = dce_tumor(i_tumor,:);
int_muscle = dce_muscle(i_muscle,:);

for j=1:size(int_tumor,1)
    int_tumor(j,:)=cumtrapz(int_tumor(j,:));
end

for j=1:size(int_muscle,1)
    int_muscle(j,:)=cumtrapz(int_muscle(j,:));
end

%% Classification

% Prepare data for classification
dce_table_tumor = array2table(int_tumor);
dce_table_muscle = array2table(int_muscle);
dce_table_tumor.response = dce_classes.Cell_Lines(i_tumor); % Response here was chosen to be the cancer cell line
dce_table_muscle.response = dce_classes.Cell_Lines(i_muscle);

% Pass data to a classifier or use the Classification Learner App
[metrics_tumor_PCA,CM_tumor_PCA]=trainLDA_PCA(dce_table_tumor(:,1:150),dce_table_tumor.response,5,size(dce_table_tumor,1),'Hs766T');
[metrics_tumor_noPCA,CM_tumor_noPCA]=trainLDA(dce_table_tumor(:,1:150),dce_table_tumor.response,size(dce_table_tumor,1),'Hs766T');