% circadian metric

% clear;clc;
load('allFilterParams.mat');
G = 100-G;
fig = figure();
fig.Units = "centimeters";
fig.Position = [0 0 8 8];
scatter(C,L,30,'x','LineWidth',1.5);hold on
plot(0:100,100:-1:0,'k--')
xlabel('Maximum Circadian Reduction')
ylabel('Gamut Size')
axis square
% set(gca,'yscale','log','xscale','log')
ax = gca;
ax.FontSize = 8;
%line = y = -x + 100, so 1(y) 1(x) - 100 = 0
%closest point on the line:

% a = 1
% b = 1
% c = -100
% xs = ((C-G)+100)./2;
% ys = ((-C+G)+100)./2;

gamutDistance = (C + G - 100)./sqrt(2);
[~,maxIdx] = max(gamutDistance);
scatter(C(gamutDistance>0),G(gamutDistance>0),30,'x','LineWidth',1.5);hold on
scatter(C(maxIdx),G(maxIdx),30,'x','LineWidth',1.5);hold on

%only 52 of 120 filters perform better than unity.

% exportgraphics(gca,'~/Desktop/gamutCorrelation.png')
