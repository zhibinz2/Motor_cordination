% cd /home/zhibin/Documents/GitHub/Motor_cordination
! git status
! git add .
CommitName=char(datetime('now'));
!git commit -a -m "CommitName"
% https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
% https://docs.github.com/en/get-started/getting-started-with-git/caching-your-github-credentials-in-git
% My token: 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
% !git push 
%%
% undo all current changes
% !git restore .

%get rid of untracked files
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