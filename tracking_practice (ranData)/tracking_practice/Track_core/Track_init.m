function [x, par] = Track_init()

par.kdown = 2;
par.trialnum = 1;
par.stepsize = [.01 .005 .005 .002];
par.tracklimit = [0 0.3];
par.x0 = 0.05;   %the initial deltaf/f

x = par.x0;
par.cor = [];
par.direction = 0;
par.cor_count = 0;
par.rev_count = 0;
end



