function [dx, rev_flag, par] = Track_update(par)

rev_flag = 0;
dx = 0;

if par.cor == 1
    par.cor_count = par.cor_count+1;
    if par.cor_count == par.kdown
        if par.direction == 1
            par.rev_count = par.rev_count+1;
            rev_flag = 1;
        end
        par.direction = -1;
        dx = - par.stepsize(min(par.rev_count+1,length(par.stepsize)));
        par.cor_count = 0;
    end
elseif par.cor == 0
    if par.direction == -1
        par.rev_count = par.rev_count+1;
        rev_flag = 1;
    end
    par.direction = 1;
    dx = par.stepsize(min(par.rev_count+1,length(par.stepsize)));
    par.cor_count = 0;
end
    
par.trialnum = par.trialnum+1;

end