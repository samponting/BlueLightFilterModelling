function [L,M,W,G] = getFilterParameters(filter)

% Pass in filter, return Luminance percentage, melanopic percentage, white
% shift and gamut reduction.

d65 = readmatrix('Illuminantd65.csv');
vLam = readmatrix('linCIE2008v10e_1.csv');
vLam = vLam(11:311,:);
d65 = d65(101:401,:);

[CF,wls] = GetCIES026;
CF = CF(:,21:321);
%% L
nFilLum = vLam(:,2)'*d65(:,2);
d65Filtered = d65(:,2) .* filter;
filLum = vLam(:,2)'*d65Filtered;
L = 100*filLum/nFilLum;

%% M
nFilMel = CF(5,:) * d65(:,2);
filMel = CF(5,:)* d65Filtered;
M = 100*filMel/nFilMel;

%% W

T_xyz = readmatrix('ciexyz64_1.csv');
T_xyz = 683*[T_xyz(find(T_xyz(:, 1) == 400):find(T_xyz(:, 1) == 700), 2:4)];
T_xyz = T_xyz';

XYZ = T_xyz*[d65Filtered d65(:,2)];
uvY = XYZTouvY(XYZ);

W = pdist2(uvY(1:2,1)',uvY(1:2,2)','euclidean');

%% G


ref = readmatrix('99Reflectances.xlsx');
ref2 = SplineSpd([380 5 81],ref,(400:700)');

refXYZ = T_xyz*ref2;
refuvY = XYZTouvY(refXYZ);
[~,area] = convhull(refuvY(1,:),refuvY(2,:));

gamFil = ref2.*filter;

filXYZ = T_xyz*gamFil;
filuvY = XYZTouvY(filXYZ);
[~,filArea] = convhull(filuvY(1,:),filuvY(2,:));

G = 100-100*filArea/area; %this is equal to amount reduced. lower is better.

end
