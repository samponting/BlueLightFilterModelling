
T =readmatrix('opo12648-sup-0003-datas2.csv');
for i = 2:121
transmittance(:,i-1) = T(i,23:323)'/100;transmittance(isnan(transmittance)) = 0;
end
% transmittance = transmittance(:,1:end-1);

[coeff, score, latent, tsquared, explained,PCAmean] = pca(transmittance');


% plot(400:700,coeff(:,1:3));

%%

% Xcentered = score*coeff' + repmat(PCAmean,120,1);
% 
% plot(Xcentered(1:4,:)');hold on
% plot(transmittance(:,1:4),'b');
% 
% 
% tmp = score(1,1:50) * coeff(:,1:50)' + PCAmean;

%%
% initialise 'score' to... 0.5s?
numVars = 5;
x0 = 0*ones(1,numVars);
ub = max(score(:,1:numVars));
lb = min(score(:,1:numVars));
wl = 400:700;

DRcurve = readmatrix("DRcurve.csv");
d65 = readmatrix('illuminantd65.csv');d65 = d65(101:401,:);
dataKey = readtable('dataKey.csv');
load('T_CIE_Y10.mat');
vLam = readmatrix('linCIE2008v10e_1.csv');
T_xyz = readmatrix('ciexyz64_1.csv');
ref = readmatrix('99Reflectances.xlsx');
spectra = load('spectrum.mat');

optim = createOptimProblem('fmincon', ...
    'objective',@(x)objfunc_PCAmaxCircValues(x,coeff(:,1:numVars),PCAmean,DRcurve,d65,dataKey,T_CIE_Y10,vLam,T_xyz,ref,spectra),...
    'x0', x0,'ub', ub,'lb', lb,...
    'options',optimoptions(@fmincon,'Algorithm','interior-point','Display','iter'));

gs = GlobalSearch('Display', 'on');
[xs,fval] = run(gs,optim);
