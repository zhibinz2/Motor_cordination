%% save all variables from the workspace in a matfile named with the date
filename=[num2str(seed) '.mat'];
cd behaviraldata/
save(filename);
