% Filter parameters
%%
fig = figure();
t = tiledlayout(2,2);
t.Padding = 'tight';
t.TileSpacing = 'compact';
nexttile();
fig.Units = 'centimeters';
fig.Position = [0,0,12.5,12.5];
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 12.5 12.5];
fig.PaperType = 'a4';
ax = gca;
colors = colorcet('C9');
f = readmatrix('filter_transmittances.xlsx');
plot(f(:,1),f(:,2),'LineWidth',2,'Color',colors(1,:));hold on
plot(f(:,1),f(:,3),'LineWidth',2,'Color',colors(100,:));
legend('Essilor','Reticare','Location','southeast')
xlim([380 1080])
xlabel('Wavelength (nm)')
ylim([0 1])
ylabel('Transmittance')
box on
set(ax,'FontSize',10,'FontName','Arial')
ax.PlotBoxAspectRatio = [1 1 0.5];

%%

nexttile();
d65 = readmatrix('Illuminantd65.csv');
vLam = readmatrix('linCIE2008v10e_1.csv');
vLam = vLam(1:2:end,:);
d65 = d65(1:2:end,:);
d65 = d65(46:end,:);
f2 = f(6:226,:);
[CF,wls] = GetCIES026;
CF = CF(:,11:2:end);

nFilLum = vLam(:,2)'*d65(:,2);
filSpec = d65(:,2) .* f2(:,2:3);
filLum = vLam(:,2)'*filSpec;
lumRatio = 100*filLum/nFilLum;

nFilMel = CF(5,:) * d65(1:196,2);
filMel = CF(5,:)* filSpec(1:196,:);
melRatio = 100*filMel/nFilMel;

Filters = {'Essilor';'Reticare'};
% T = table(lumRatio',melRatio','RowNames',Filters);
% 
ax = gca;
% ax.YTick = [];
% ax.XTick = [];
box off
axis off
ax.Alphamap=0;
% 
% ax.BoxStyle = "full";
% ax.LineWidth = 1;
% ax.PlotBoxAspectRatio = [1 1 1];
% ax = gca;
% 
% fs = 15;
% annotation(gcf,'Textbox','String',sprintf('%s%%',num2str(round(lumRatio(1)))),'Position',[ax.Position(1) ax.Position(2)+0.5*ax.Position(4) 0.5*ax.Position(3) 0.5*ax.Position(4)],'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontName','Arial','FontSize',fs,'FontWeight','bold','LineWidth',2)
% annotation(gcf,'Textbox','String',sprintf('%s%%',num2str(round(lumRatio(2)))),'Position',[ax.Position(1) ax.Position(2) 0.5*ax.Position(3) 0.5*ax.Position(4)],'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontName','Arial','FontSize',fs,'FontWeight','bold','LineWidth',2)
% annotation(gcf,'Textbox','String',sprintf('%s%%',num2str(round(melRatio(1)))),'Position',[ax.Position(1)+0.5*ax.Position(3) ax.Position(2)+0.5*ax.Position(4) 0.5*ax.Position(3) 0.5*ax.Position(4)],'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontName','Arial','FontSize',fs,'FontWeight','bold','LineWidth',2)
% annotation(gcf,'Textbox','String',sprintf('%s%%',num2str(round(melRatio(2)))),'Position',[ax.Position(1)+0.5*ax.Position(3) ax.Position(2) 0.5*ax.Position(3) 0.5*ax.Position(4)],'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontName','Arial','FontSize',fs,'FontWeight','bold','LineWidth',2)
% 
% annotation(gcf,'Textbox','String','Essilor','Position',[ax.Position(1)+0.015 ax.Position(2)+0.52*ax.Position(4) 0.25*ax.Position(3) 0.25*ax.Position(4)], ...
%     'FontName','Arial','FontSize',8,'LineWidth',2,'Rotation',90,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','left','VerticalAlignment','top')
% annotation(gcf,'Textbox','String','Reticare','Position',[ax.Position(1)+0.015 ax.Position(2) 0.25*ax.Position(3) 0.25*ax.Position(4)], ...
%     'FontName','Arial','FontSize',8,'LineWidth',2,'Rotation',90,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','left','VerticalAlignment','top')
% 
% annotation(gcf,'Textbox','String',sprintf('Luminous\nTransmittance'),'Position',[ax.Position(1) ax.Position(4) 0.5*ax.Position(3) 0.65*ax.Position(4)], ...
%     'FontName','Arial','FontSize',8,'LineWidth',2,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','top')
% annotation(gcf,'Textbox','String',sprintf('Melanopic\nTransmittance'),'Position',[ax.Position(1)+0.5*ax.Position(3) ax.Position(4) 0.5*ax.Position(3) 0.65*ax.Position(4)], ...
%     'FontName','Arial','FontSize',8,'LineWidth',2,'LineStyle','none','FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','top')

%%

nexttile();
T_xyz = readmatrix('ciexyz64_1.csv');
T_xyz = 683*[T_xyz(find(T_xyz(:, 1) == 390):2:find(T_xyz(:, 1) == 780), 2:4)];
T_xyz = T_xyz';

XYZ = T_xyz*[filSpec(1:196,:) d65(1:196,2)];
uvY = XYZTouvY(XYZ);


theWls = [780 390:2:780];
for ii = 1:length(theWls)
    spd = zeros(size(390:2:780));
    spd(theWls(ii) == 390:2:780) = 1;
    
    % Get ref
    XYZref = (T_xyz*spd');
    uvYref = XYZTouvY(XYZref);
    uvref(:, ii) = uvYref(1:2, :);
end
% DrawChromaticity('upvp'); hold on
plot(uvref(1,:),uvref(2,:),':k');hold on
scatter(uvY(1,1),uvY(2,1),30,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(1,:))
scatter(uvY(1,2),uvY(2,2),30,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(100,:))
scatter(uvY(1,3),uvY(2,3),30,'filled','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0.5,0.5,0.5]);
% legend('','Essilor','Reticare','Location','southeast')
xlim([0.1 0.35])
ylim([0.4 0.6])
xlabel(sprintf("u'"))
ylabel(sprintf("v'"))
ax = gca;
box on
set(ax,'FontSize',10,'FontName','Arial')
ax.PlotBoxAspectRatio = [1 1 1];

sqrt(sum(uvY(1:2,3)-uvY(1:2,2).^2))
sqrt(sum(uvY(1:2,3)-uvY(1:2,1).^2))

pdist2(uvY(1:2,3)',uvY(1:2,1)','euclidean')
pdist2(uvY(1:2,3)',uvY(1:2,2)','euclidean')
%%

nexttile();

ref = readmatrix('99Reflectances.xlsx');
ref2 = SplineSpd([380 5 81],ref,(390:2:780)');

refXYZ = T_xyz*ref2;
refuvY = XYZTouvY(refXYZ);
[idx,area] = convhull(refuvY(1,:),refuvY(2,:));

scatter(refuvY(1,:),refuvY(2,:),'o','MarkerEdgeColor',[0.7,0.7,0.7],'LineWidth',1.5);hold on
plot(refuvY(1,idx),refuvY(2,idx),'Color',[0.7,0.7,0.7]);

gamFil(:,:,1) = ref2.*f2(1:196,2);
gamFil(:,:,2) = ref2.*f2(1:196,3);

filXYZ1 = T_xyz*gamFil(:,:,1);
filuvY1 = XYZTouvY(filXYZ1);
[idx1,area1] = convhull(filuvY1(1,:),filuvY1(2,:));

filXYZ2 = T_xyz*gamFil(:,:,2);
filuvY2 = XYZTouvY(filXYZ2);
[idx2,area2] = convhull(filuvY2(1,:),filuvY2(2,:));

scatter(filuvY1(1,:),filuvY1(2,:),'o','MarkerEdgeColor',colors(1,:),'LineWidth',1.5);hold on
plot(filuvY1(1,idx),filuvY1(2,idx),'Color',colors(1,:));
scatter(filuvY2(1,:),filuvY2(2,:),'o','MarkerEdgeColor',colors(100,:),'LineWidth',1.5);hold on
plot(filuvY2(1,idx),filuvY2(2,idx),'Color',colors(100,:));
xlabel(sprintf("u'"))
ylabel(sprintf("v'"))
ax = gca;
box on
set(ax,'FontSize',10,'FontName','Arial')
ax.PlotBoxAspectRatio = [1 1 1];


% exportgraphics(fig,'~/Desktop/filterSpec.png','Resolution',600)



