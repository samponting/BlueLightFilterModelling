
T =readtable('opo12648-sup-0003-datas2.csv');

for i = 2:121
transmittance(:,i-1) = table2array(T(i,23:323))'/100;transmittance(isnan(transmittance)) = 0;
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

optim = createOptimProblem('fmincon', ...
    'objective',@(x)objfunc_PCAmaxCircValues(x,coeff(:,1:numVars),PCAmean),...
    'x0', x0,'ub', ub,'lb', lb,...
    'options',optimoptions(@fmincon,'Algorithm','interior-point','Display','iter'));

gs = GlobalSearch('Display', 'on');
[xs,fval] = run(gs,optim);
