function C = getFilterCircadianResponse(transmission)

addpath(genpath([pwd '/spectra_calibrated']))
DRcurve = readmatrix("DRcurve.csv");
d65 = readmatrix('illuminantd65.csv');d65 = d65(101:401,:);
dataKey = readtable('dataKey.csv');
fieldLeng = 0;
for folders = 1:length(dataKey.FolderName)
    files = dir([pwd,'/spectra_calibrated/',cell2mat(dataKey{folders,"FolderName"}),'/spectra/*.csv']);
    numFiles = length(files(dataKey{folders,"FirstField"}:dataKey{folders,"LastField"}));

    fieldFiles(fieldLeng+1:fieldLeng+numFiles) = files(dataKey{folders,"FirstField"}:dataKey{folders,"LastField"});
    fieldLeng = numFiles+fieldLeng;
end
load('T_CIE_Y10.mat');
T_CIE_Y10 = [T_CIE_Y10(:,11:311)];

addpath('~/Documents/Github/AsanoObserverModel/')
age = 32;
CF = generateConeFundamentals(age,10)';

%calculate d65 lux
lumd65 = d65(:,2)'*T_CIE_Y10'*683.002;
%calculate d65 melanopic irradiance
mIrrd65 = d65(:,2)'*CF(5,:)';
%calculate d65 melanopic effective luminous ratio (from CIES026)
mELRd65 = 1.3262;
%calculate d65 melanopic equivalent daylighy illuminance
mEDId65 = 1000*mIrrd65/mELRd65;


lum = zeros(fieldLeng,2);
mIrr = zeros(fieldLeng,2);
mEDI = zeros(fieldLeng,2);
% column 1 = unfiltered, column 2 = filtered

spectra = load('spectrum.mat');spectrum = spectra.spectrum;
filteredSpectrum = spectrum .* transmission;
%calculate test lux
lum(:,1) = spectrum'*T_CIE_Y10'*683.002;
lum(:,2) = filteredSpectrum'*T_CIE_Y10'*683.002;
%calculate test melanopic illuminance
mIrr(:,1) = spectrum'*CF(5,:)';
mIrr(:,2) = filteredSpectrum'*CF(5,:)';
%calculate test mEDI
mEDI(:,1) = log(mIrr(:,1)/mELRd65);
mEDI(:,2) = log(mIrr(:,2)/mELRd65);


f = fit(DRcurve(:,1),DRcurve(:,2),'d+a/(1+exp(-b*(x-c)))','StartPoint',[97, 2, 1, 1]);

response_noFilter = feval(f,mEDI(:,1));
response_filter = feval(f,mEDI(:,2));
delta = response_noFilter - response_filter;
mEDI(isinf(mEDI)) = -10;
[~,sIdx] = sort(mEDI(:,1));

order = 3;
framelen = 21;
sgf = sgolayfilt(delta(sIdx),order,framelen);

C = max(sgf);


end