% plot melanopic / luminous transmittance vs. peak circadian reduction
%% Init
clc;close all;clear
DRcurve = readmatrix("DRcurve.csv");
f = fit(DRcurve(:,1),DRcurve(:,2),'d+a/(1+exp(-b*(x-c)))');

%% Dose-Response Curve

T =readtable('opo12648-sup-0003-datas2.csv');
L = zeros(size(T,1)-1,1);
M = zeros(size(T,1)-1,1);
W = zeros(size(T,1)-1,1);
G = zeros(size(T,1)-1,1);
C = zeros(size(T,1)-1,1);

for filterNumber = 2:121
    fprintf('filter number = %d\n',filterNumber)
    load(sprintf('case1_%s_SpectrumMetrics.mat',num2str(filterNumber)));
    transmission = table2array(T(filterNumber,23:323))'/100;transmission(isnan(transmission)) = 0;


    mEDIlin = 10.^mEDI(:,1);
    response_noFilter = feval(f,mEDI(:,1));
    response_filter = feval(f,mEDI(:,2));
    delta = response_noFilter - response_filter;
    mEDI(isinf(mEDI)) = -10;
    [mEDIsort,sIdx] = sort(mEDI(:,1));

    order = 3;
    framelen = 21;
    
    sgf = sgolayfilt(delta(sIdx),order,framelen);
    C(filterNumber-1) = max(sgf);

    [L(filterNumber-1),M(filterNumber-1),W(filterNumber-1),G(filterNumber-1)] = getFilterParameters(transmission);

%     filterParams.L = L;
%     filterParams.M = M;
%     filterParams.W = W;
%     filterParams.G = G;
%     filterParams.C = circPeak;
%     save([pwd,'/filtParams/',num2str(filterNumber),'_filterParams.mat'],'filterParams');

end
save([pwd,'/allFilterParams.mat'],'C','G','W','M','L')
