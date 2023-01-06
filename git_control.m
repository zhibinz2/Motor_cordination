% cd /home/zhibinz2/Documents/GitHub/Motor_cordination % hnlb
% cd /home/hnl/Documents/GitHub/Motor_cordination % hnlstim2
! git status
! git add .
CommitName=char(datetime('now'));
!git commit -a -m CommitName
% https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
% https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
% My zhibin bash script — repo token: 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
% !git push https://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@github.com/zhibinz2/Motor_cordination.git

%%
% undo all current changes
% !git restore .

% get rid of untracked files
% !git clean -f

%%
% save changes temporarily
% !git stash --include-untracked

% to restore the stash
% !git stash pop

%% Make Git store the username and password and it will never ask for them.
% !git config --global credential.helper store
 
% % Save the username and password for a session (cache it);
% !git config --global credential.helper cache
% % set a timeout for the above setting
% !git config --global credential.helper 'cache --timeout=600'
 
% !git push

%% open github desktop using SSH from home 
% !pkill github
% !github