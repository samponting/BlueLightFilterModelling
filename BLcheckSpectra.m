% check spectra

clear;clc;close all
dataKey = readtable('dataKey.csv');

fieldLeng = 0;
modes = ["Dark","Lab","Field"];

for folders = 1:length(dataKey.FolderName)
    figure();
    tiledlayout(1,3);
    for m = 1:length(modes)
        nexttile();
        files = dir([pwd,'/spectra_calibrated/',cell2mat(dataKey{folders,"FolderName"}),'/spectra/*.csv']);
        numFiles = length(files(dataKey{folders,sprintf("First%s",modes(m))}:dataKey{folders,sprintf("Last%s",modes(m))}));
    
        fieldFiles = files(dataKey{folders,sprintf("First%s",modes(m))}:dataKey{folders,sprintf("Last%s",modes(m))});

        for spectra = 1:length(fieldFiles)
            fprintf('processing spectrum %s of %s\n',num2str(spectra),num2str(length(fieldFiles)))
            spectrum = readmatrix(fieldFiles(spectra).name);
            plot(spectrum(:,1),spectrum(:,2),'LineWidth',0.5);hold on
        end
        title(modes(m))
        ylim([0 1]);
    end
end