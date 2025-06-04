%% Parametric filter optimisation 

% 1. Define filter using a logistic function
wl = 400:700;
% asymtote = 1;
% slope = .5;
% cutoff = 450;

% transmission = asymtote./(1+exp(-slope.*(wl-cutoff)));
% plot(wl,transmission)

% 2. Optimise
wlsp = 400;

optim = createOptimProblem('fmincon', ...
    'objective',@(x)objfunc_maximiseCircValue(x,wl),...
    'x0', [1, .1, (wlsp-wl(1))/(wl(end)-wl(1))],'ub',[1, 1, .8],'lb',[0, 0, 0.1],...
    'options',optimoptions(@fmincon,'Algorithm','sqp','Display','iter'));

gs = GlobalSearch('Display', 'on');
[xs,fval] = run(gs,optim);
xs(3) = (xs(3)*(wl(end)-wl(1))) + wl(1);



