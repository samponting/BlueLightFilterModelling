function x = objfunc_maximiseCircValue(xs,wl)

% xs(xs>1) = 1;
% xs(xs<0) = 0;
% wl = 400:700;
xs(3) = (xs(3)*(wl(end)-wl(1))) + wl(1);

transmission = xs(1)./(1+exp(-xs(2).*(wl-xs(3))))';
[L,~,~,G] = getFilterParameters(transmission);
G = 100-G;
C = getFilterCircadianResponse(transmission);

gamutDiff = (C + G - 100)./sqrt(2);

lumDiff = (C + L - 100)./sqrt(2);

x =  gamutDiff + lumDiff;

end