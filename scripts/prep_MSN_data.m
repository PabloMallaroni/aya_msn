% Ayahuasca MSA Data Preparation

%This script will:
% Extract all drug and control freesurfer stats
% Extract all MRIQC +freesurfer QC data.
% Generate subject,study,age and sex vectors. 
% Generate metric name vectors. 

%this is contigent on you running fslaverage to native script.


% Setup
close all;
clearvars -except paths;
clc;

% General paths
paths.home = 'D:\Aya_struct';
paths.data1 = fullfile(paths.home, 'aya_sourcedata');
paths.data2 = fullfile(paths.home, 'controls_sourcedata');
paths.out = fullfile(paths.home, 'data_clean');
paths.qc_checks = fullfile(paths.home,'MRIQC','QA_checks');


% Flags
flags.prep_QC = 1; % Merge quality checks and retain selected controls.
flags.extract_vals = 1; % Get morphometric features.
flags.save = 1;

% Configurations
configs.input = {paths.data1, paths.data2};
configs.num_parcels = 308;
configs.num_perps = 48;

% Data Labels
configs.data_labels = {'number of vertices', 'total surface area (mm^2)', ...
    'total gray matter volume (mm^3)', 'average cortical thickness +- standard deviation (mm)', ...
    'integrated rectified mean curvature', 'integrated rectified Gaussian curvature', ...
    'folding index', 'intrinsic curvature index', 'structure name'};

% Create identifiers
aya_info = xlsread(fullfile(paths.home, 'demo_match.xlsx'), 'aya_data');
control_info = xlsread(fullfile(paths.home, 'demo_match.xlsx'), 'control_data');
subject = [aya_info(:,1); control_info(:,1)];
group = [aya_info(:,2); control_info(:,2)];
sex = [aya_info(:,3); control_info(:,3)];
age = [aya_info(:,4); control_info(:,4)];

% Extract Data
if flags.extract_vals == 1
    % Initialize variables
    num_subjects = length(subject);
    data_labels = {'SA', 'GM', 'CT', 'MC', 'GC', 'FI', 'CI'};
    num_measures = length(data_labels);
    data = NaN(num_subjects, configs.num_parcels, num_measures);

    for folder = 1:length(configs.input)
        subjectList = dir(configs.input{folder});
        subjectList(1:3) = [];

        for sub = 1:length(subjectList)
            data_dir = fullfile(configs.input{folder}, subjectList(sub).name, 'stats');
            if exist(fullfile(data_dir, 'lh.500.aparc.log'), 'file') && ...
                    exist(fullfile(data_dir, 'rh.500.aparc.log'), 'file') && ...
                    exist(fullfile(data_dir, 'aseg.stats'), 'file')
                % Extract total intracranial number and euler number.
                seg_stats = readtable(fullfile(data_dir, 'aseg.stats'), 'FileType', 'text', 'ReadVariableNames', false);
                data(sub, :, 1) = seg_stats{contains(seg_stats{:, 2}, 'Measure EstimatedTotalIntraCranialVol,'), 3:end};
                data(sub, :, 2) = seg_stats{contains(seg_stats{:, 2}, 'Measure SurfaceHoles, SurfaceHoles,'), 3:end};

                % Extract all 7 measures of interest from log files
                for i = 1:num_measures
                    lh_stats = readtable(fullfile(data_dir, ['lh.500.aparc.log']), 'FileType', 'text', 'ReadVariableNames', false);
                    rh_stats = readtable(fullfile(data_dir, ['rh.500.aparc.log']), 'FileType', 'text', 'ReadVariableNames', false);
                    lh_stats(1, :) = [];
                    rh_stats(1, :) = [];
                    data(sub, :, i + 2) = [table2array(lh_stats(:, i + 1))'; table2array(rh_stats(:, i + 1))'];
                end
            else
                disp(['Missing data for: ', subjectList(sub).name])
            end
        end
    end
    disp('Data extraction completed.');
end

% Clean and merge MRIQC outputs with Euler and TCIF
if flags.prep_QC == 1
    data_QC = [];
    [~, ~, raw] = tsvread(fullfile(paths.qc_checks, 'ayahuasca', 'group_T1w.tsv'));
    data_QC_labels = raw(1, 2:end);
    data_QC = str2double(raw(2:end, 2:end));
    
    [~, ~, raw] = tsvread(fullfile(paths.qc_checks, 'controls', 'group_T1w.tsv'));
    raw = raw(2:end, 2:end);
    raw = str2double(raw(control_info(:, 1), :));
    data_QC = [data_QC; raw];
    
    data_QC = [data_QC, data(:, :, 2), data(:, :, 1)];
    data_QC_labels{end + 1} = 'Euler_raw';
    data_QC_labels{end + 1} = 'TIV';
    disp('QC data preparation completed.');
end

% Export Data as .dat Files
if flags.save == 1
    if ~exist(paths.out, 'dir')
        mkdir(paths.out);
    end
    
    data_labels = [configs.data_labels, 'Euler_raw', 'TIV'];
    dat_vars = {subject, group, sex, age, data(:, :, 1), data(:, :, 2), data(:, :, 3), ...
        data(:, :, 4), data(:, :, 5), data(:, :, 6), data(:, :, 7), data(:, :, 8), data_QC};
    dat_filenames = {'subject.dat', 'group.dat', 'sex.dat', 'age.dat', 'PARC500_SA.dat', 'PARC500_GM.dat', ...
        'PARC500_CT.dat', 'PARC500_MC.dat', 'PARC500_GC.dat', 'PARC500_FI.dat', 'PARC500_CI.dat', ...
        'QC.dat'};
    
    for i = 1:length(dat_vars)
        dlmwrite(fullfile(paths.out, dat_filenames{i}), dat_vars{i});
    end
    
    save(fullfile(paths.out, 'QC_labels.mat'), 'data_QC_labels');
    disp('Data saved successfully.');
end
