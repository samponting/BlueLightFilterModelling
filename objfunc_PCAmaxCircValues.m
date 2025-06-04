function x = objfunc_PCAmaxCircValues(xs,coeff,PCAmean)
disp(datetime('now','Format','HH:mm:ss.SSS'))
transmission = (xs * coeff' + PCAmean)';
transmission(transmission>1) = 1;
transmission(transmission<0) = 0;
[L,~,~,G] = getFilterParameters(transmission);
G = 100-G;
C = getFilterCircadianResponse(transmission);

gamutDiff = -(C + G - 100)./sqrt(2);
lumDiff = -(C + L - 100)./sqrt(2);


% x =  gamutDiff + lumDiff;
% x = gamutDiff;
x = lumDiff;

plot(0:100,100:-1:0,'k--');hold on
scatter(C,L,30,'x','LineWidth',1.5);drawnow
% disp('end pass')
end