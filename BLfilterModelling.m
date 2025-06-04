%% Single Blue Light Blocking Filter Analysis
% Circadian effects of blue light blocking filters - technical note.
% Methodology:
% 1. Load in dose-response curve
% 2. Calculate melanopic lux of SPD sample set
% 3. Filter SPD sample set by filters
% 4. Calculate new melanopic lux

clear;clc;
addpath(genpath([pwd '/spectra_calibrated']))
%% Load in dose-response curve
showFigs = 1;
colors = colorcet('C9');

DRcurve = readmatrix("DRcurve.csv");
mEDIlin = 10.^DRcurve(:,1);
%1 - melanopic lux (log)
%2 - normalised response
%3 - neg 95% CI
%4 - pos 95% CI

if showFigs
    fig=figure('Position',[0,0,1920,1080]);
    tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    nexttile([2,1]);
    fsize = 20;
    p = semilogx(mEDIlin,DRcurve(:,2),'LineWidth',3,'Color',colors(1,:));hold on
    fill([mEDIlin(:,1)' flip(mEDIlin(:,1))' mEDIlin(1,1)], [(DRcurve(:,2)-DRcurve(:,3))', flip(DRcurve(:,2)+DRcurve(:,4))', DRcurve(1,2)-DRcurve(1,3)],colors(1,:),'FaceAlpha',0.3)
    box on;grid on;grid minor;
    xlim([10.^-3 10.^4])
    ax = gca;ax.FontSize = fsize;ax.FontName = 'Arial';
    title('Global Dose Response Curve')
    xlabel('Melanopic EDI (log)');
    ylabel('Normalised Response (%)')
end

%% Model filter transmission

%load CIE spectra
% CF = GetCIES026;CF(isnan(CF)) = 0;
addpath('~/Documents/Github/AsanoObserverModel/')
age = 32;
CF = generateConeFundamentals(age,10)';

wl = 1:301;
cutoff = 480-380;
asymtote = 0.95;
slope = 0.5;

T =readtable('opo12648-sup-0003-datas2.csv');
f1 = table2array(T(51,23:323))'/100;

% filt = readmatrix('filter_transmittances.xlsx');
% f1(:,1) = interp1(filt(:,1),filt(:,2),400:700);
transmission = f1;
% transmission = asymtote./(1+exp(-slope*(wl-cutoff)));
if showFigs
    nexttile([2,1]);
    plot(400:700,CF(5,:),'LineWidth',2,'Color',[0,0.5,0.5]);hold on
    plot(400:700,transmission,'LineWidth',2,'Color',[0.8,0.8,0.2])
    ax = gca;ax.FontSize = fsize;ax.FontName = 'Arial';
    title('Example Blue Light Filter')
    xlabel('Wavelength (nm)');
    ylabel('Transmittance (%)')
    xlim([380 700]);
    ylim([0 1.2])
end
%% Calculate mEDI values

d65 = readmatrix('illuminantd65.csv');d65 = d65(101:401,:);%d65 = d65;%/80726;
dataKey = readtable('dataKey.csv');
 
fieldLeng = 0;
for folders = 1:length(dataKey.FolderName)
    files = dir([pwd,'/spectra_calibrated/',cell2mat(dataKey{folders,"FolderName"}),'/spectra/*.csv']);
    numFiles = length(files(dataKey{folders,"FirstField"}:dataKey{folders,"LastField"}));

    fieldFiles(fieldLeng+1:fieldLeng+numFiles) = files(dataKey{folders,"FirstField"}:dataKey{folders,"LastField"});
    fieldLeng = numFiles+fieldLeng;
end

%load Vlambda
load('T_CIE_Y10.mat');
T_CIE_Y10 = [T_CIE_Y10(:,11:311)];

%calculate d65 lux
lumd65 = d65(:,2)'*T_CIE_Y10'*683.002;
%calculate d65 melanopic irradiance
mIrrd65 = d65(:,2)'*CF(5,:)';
%calculate d65 melanopic effective luminous ratio (from CIES026)
mELRd65 = 1.3262;
%calculate d65 melanopic equivalent daylighy illuminance
mEDId65 = 1000*mIrrd65/mELRd65;

lum = zeros(length(files),2);
mIrr = zeros(length(files),2);
mEDI = zeros(length(files),2);
% column 1 = unfiltered, column 2 = filtered

for spectra = 1:length(fieldFiles)
    fprintf('processing spectrum %s of %s\n',num2str(spectra),num2str(length(fieldFiles)))
    spectrum = readmatrix(fieldFiles(spectra).name);spectrum = spectrum(21:321,:);
    filteredSpectrum = spectrum;filteredSpectrum(:,2) = filteredSpectrum(:,2) .* transmission;
    %calculate test lux
    lum(spectra,1) = spectrum(:,2)'*T_CIE_Y10'*683.002;
    lum(spectra,2) = filteredSpectrum(:,2)'*T_CIE_Y10'*683.002;

    %calculate test melanopic illuminance
    mIrr(spectra,1) = spectrum(:,2)'*CF(5,:)';
    mIrr(spectra,2) = filteredSpectrum(:,2)'*CF(5,:)';

    %calculate test mEDI
    mEDI(spectra,1) = log(mIrr(spectra,1)/mELRd65);
    mEDI(spectra,2) = log(mIrr(spectra,2)/mELRd65);
   
end

% save('spectrumMetrics.mat','mEDI','mIrr','lum');

%% Calculate responses

f = fit(DRcurve(:,1),DRcurve(:,2),'d+a/(1+exp(-b*(x-c)))','StartPoint',[97, 2, 1, 1]);

load('spectrumMetrics.mat');
response_noFilter = feval(f,mEDI(:,1));
response_filter = feval(f,mEDI(:,2));

% if showFigs
%     scatter(mEDI(:,1),response_noFilter);hold on
%     scatter(mEDI(:,2),response_Filter);hold on
% end


%% Calculate response delta

nexttile([2,1]);

maxCutoff = 29;

mEDIlin = 10.^mEDI(:,1);

delta = response_noFilter - response_filter;
% delta = (1-response_filter./response_noFilter).*100;


scatter(mEDIlin(:,1),delta,'+','LineWidth',0.3,'MarkerEdgeColor',colors(100,:));hold on
% scatter(mEDIlin(delta>maxCutoff,1),delta(delta>maxCutoff),'+','LineWidth',0.3,'MarkerEdgeColor',[0.8 0.2 0.2])
% scatter(mEDIlin((mEDIlin(:,1)<80 & mEDIlin(:,1) > 40 & delta<26.5),1),delta((mEDIlin(:,1)<80 & mEDIlin(:,1) > 40 & delta<26.5),1),'+','LineWidth',0.3,'MarkerEdgeColor',[0.2 0.2 0.8])


set(gca,'xscale','log')
ax = gca;
xlim([10^-3 10^4])
ylim([0 15])
box on;grid on;grid minor;
title('Response Delta')
ylabel('Percentage Change')
xlabel('Melanopic EDI (log)')
ax.FontSize = fsize;
ax.FontName = 'Arial';

%% Histogram of mEDI

nexttile()

histogram(delta,'FaceColor',[0.3 0.7 0.7],'BinWidth',1)
ax = gca;
% title('Maximally affected spectra')
ax.FontSize = fsize;
ax.FontName = 'Arial';
xlabel('Delta');
xlim([0 55])
ylabel('Frequency')
ylim([0 10000])
set(ax,'yscale','log')

nexttile()

histogram(mEDIlin,'FaceColor',[0.7 0.3 0.3],'BinWidth',100)
ax = gca;
% title('Maximally affected spectra')
ax.FontSize = fsize;
ax.FontName = 'Arial';
xlabel('mEDI');
xlim([0 7000])
ylabel('Frequency')
ylim([0 100000])

set(ax,'yscale','log')



% %% Segment data - maximally affected
% nexttile();
% targetIdx = delta>maxCutoff;
% 
% targetFiles = fieldFiles(targetIdx);
% 
% for t = 1:length(targetFiles)
%     spd = readmatrix(targetFiles(t).name);
%     plot(spd(:,1),spd(:,2),'Color',[0.8,0.2,0.2],'LineWidth',0.01);hold on
% end
% ax = gca;
% title('Maximally affected spectra')
% xlim([380 780]);
% ax.FontSize = fsize;
% ax.FontName = 'Arial';
% xlabel('Wavelength (nm)');
% ylabel('Irradiance (W/m^2)')
% 
% %% Segment data - minimally affected
% nexttile();
% targetXIdx = (mEDIlin(:,1)<80 & mEDIlin(:,1) > 40);
% targetYIdx = delta<26.5;
% 
% targetIdx = targetXIdx & targetYIdx;
% 
% targetFiles = fieldFiles(targetIdx);
% 
% for t = 1:length(targetFiles)
%     spd = readmatrix(targetFiles(t).name);
%     plot(spd(:,1),spd(:,2),'Color',[0.2,0.2,0.8],'LineWidth',0.1);hold on
% end
% ax = gca;
% title('Minimally affected spectra')
% xlim([380 780]);
% ax.FontSize = fsize;
% ax.FontName = 'Arial';
% xlabel('Wavelength (nm)');
% ylabel('Irradiance (W/m^2)')
% 
% % saveas(gcf,'doseResponseDelta.png');
% exportgraphics(fig,sprintf('doseResponseDeltaEssilor_%s.png',num2str(age)),'Resolution',600)
% 
% 