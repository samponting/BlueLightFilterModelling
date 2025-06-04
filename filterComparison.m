%UoG filter analysis
clear;close all;clc;
filters = readmatrix('filter_transmittances.xlsx');
f(:,1) = interp1(filters(:,1),filters(:,2),400:700);
f(:,2) = interp1(filters(:,1),filters(:,3),400:700);

fig = figure();
fig.Units = 'centimeters';
fig.Position = [0,0,12.5,6];
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 12.5 6];
fig.PaperType = 'a4';
t = tiledlayout(1,2);
t.Padding = 'compact';
% nexttile();

colors = colorcet('C9');

% plot(400:700,f(:,1),'LineWidth',2,'Color',colors(100,:));hold on
% plot(400:700,f(:,2),'LineWidth',2,'Color',colors(1,:))
% ax = gca;ax.FontSize = 16;ax.FontName = 'Arial';
% title('Filter Transmittances')
% xlabel('Wavelength (nm)');
% ylabel('Transmittance (%)')
% xlim([380 700]);
% ylim([0 1])
% legend('Essilor','Reticare','Location','southeast')

% nexttile();

showFigs = 1;

DRcurve = readmatrix("DRcurve.csv");
% mEDIlin = 10.^DRcurve(:,1);
%1 - melanopic lux (log)
%2 - normalised response
%3 - neg 95% CI
%4 - pos 95% CI

addpath('~/Documents/Github/AsanoObserverModel/')


for age = [20,50]
    CF = generateConeFundamentals(age,10)';% 380:780 (1nm step)
    
    
    d65 = readmatrix('illuminantd65.csv');d65 = d65(101:401,:);%d65 = d65;%/80726;
    
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
    
%     load('24h_radiances.mat');radiances1 = radiances.*120;%380:2:1080
%     load('24h_radiances_Eva.mat');radiances2 = radiances.*120;
%     radiances = [radiances1 radiances2];
    load('Granada_daylight_2600_161.mat');
%     wls = 300:5:1100;
    radiances = SplineSpd([300 5 ((1100-300)/5)+1],final,[400 1 (700-400)+1]);
    
    lum = zeros(size(radiances,2),2,size(f,2));
    mIrr = zeros(size(radiances,2),2,size(f,2));
    mEDI = zeros(size(radiances,2),2,size(f,2));
    
    
    for filter = 1:size(f,2)
        transmission = f(:,filter);
        for spectra = 1:size(radiances,2)
            fprintf('processing spectrum %s of %s\n',num2str(spectra),num2str(size(radiances,2)))
            spectrum = radiances(11:161,spectra);spectrum = interp1(400:2:700,spectrum,400:700);
            filteredSpectrum = spectrum .* transmission';
            %calculate test lux
            lum(spectra,1,filter) = spectrum*T_CIE_Y10'*683.002;
            lum(spectra,2,filter) = filteredSpectrum*T_CIE_Y10'*683.002;
        
            %calculate test melanopic illuminance
            mIrr(spectra,1,filter) = spectrum*CF(5,:)';
            mIrr(spectra,2,filter) = filteredSpectrum*CF(5,:)';
        
            %calculate test mEDI
            mEDI(spectra,1,filter) = log(mIrr(spectra,1,filter)/mELRd65);
            mEDI(spectra,2,filter) = log(mIrr(spectra,2,filter)/mELRd65);
           
        end
    end
    
    
    DRfit = fit(DRcurve(:,1),DRcurve(:,2),'d+a/(1+exp(-b*(x-c)))','StartPoint',[97, 2, 1, 1]);
    
    
    response_noFilter(:,1) = feval(DRfit,mEDI(:,1,1));response_noFilter(:,2) = feval(DRfit,mEDI(:,1,2));
    response_filter(:,1) = feval(DRfit,mEDI(:,2,1));response_filter(:,2) = feval(DRfit,mEDI(:,2,2));
    
    
    
    mEDIlin(:,1) = 10.^mEDI(:,1,1);
    mEDIlin(:,2) = 10.^mEDI(:,1,2);
    
    
    delta(:,1) = response_noFilter(:,1) - response_filter(:,1);
    delta(:,2) = response_noFilter(:,2) - response_filter(:,2);
    
    % delta = (1-response_filter./response_noFilter).*100;
    
    nexttile();
    scatter(mEDIlin(:,1),delta(:,1),2,'^','LineWidth',0.1,'MarkerEdgeColor',colors(1,:));hold on
    scatter(mEDIlin(:,2),delta(:,2),2,'^','LineWidth',0.1,'MarkerEdgeColor',colors(100,:));hold on
    
    % scatter(mEDIlin(delta>maxCutoff,1),delta(delta>maxCutoff),'+','LineWidth',0.3,'MarkerEdgeColor',[0.8 0.2 0.2])
    % scatter(mEDIlin((mEDIlin(:,1)<80 & mEDIlin(:,1) > 40 & delta<26.5),1),delta((mEDIlin(:,1)<80 & mEDIlin(:,1) > 40 & delta<26.5),1),'+','LineWidth',0.3,'MarkerEdgeColor',[0.2 0.2 0.8])
    
    
    set(gca,'xscale','log')
    xticks([.1 10 1000])
    ax = gca;
    xlim([10^-3 10^3.5])
    ylim([0 15])
    box on;grid on;grid minor;
    title('Response Delta')
    ylabel('Percentage Change')
    xlabel('Melanopic EDI')
    ax.FontSize = 8;
    ax.FontName = 'Arial';
    annotation('textbox','String',sprintf('Age: %s',num2str(age)),'Position',[ax.Position(1) ax.Position(2) ax.Position(3) ax.Position(4)],'HorizontalAlignment','right')
% annotation('textbox','String',sprintf('Age: %s',num2str(60)),'Position',[0.6 0.31,0.1,0.1])
    
    if age ==20
        legend('Essilor','Reticare','Location','best')
    end

end



% exportgraphics(fig,['~/Desktop/ALLillumination.png'],'Resolution',600)






