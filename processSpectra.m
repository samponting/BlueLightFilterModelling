addpath(genpath([pwd,'/spectra_calibrated']));
numSpectra = 113;
folder = cell(113,1);


folder = dir([pwd,'/spectra_calibrated']);folder = folder(4:end);


for spectra = 1:length(folder)
    tic
    fprintf('processing Spectrum %d\n',spectra)
    files = dir([pwd,'/spectra_calibrated/',folder(spectra).name,'/spectra/*.csv']);
    means = zeros(length(files),1);
    for file = 1:length(files)
        s = readmatrix(files(file).name);
        means(file) = mean(s(:,2));
    end
    [~,maxIdx] = max(means);
    spectrum = readmatrix(files(maxIdx).name);
    save(sprintf('SP%03d.mat',spectra),"spectrum");
    toc
end

%%
files = dir([pwd,'/*.mat']);
for s = 1:length(files)
    load(files(s).name)

    plot(spectrum(:,1),spectrum(:,2));hold on

end

