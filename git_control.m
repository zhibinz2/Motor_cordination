! git status
! git add .
CommitName=char(datetime('now'));
!git commit -a -m "CommitName"

!git push
%%
% undo all current changes
!git restore .

%get rid of untracked files
!git clean -f

%%
% save changes temporarily
!git stash --include-untracked

% to restore the stash
!git stash pop

%% Update the URL of origin remote using SSH instead of HTTPS so that it will not ask for password when push
!git remote set-url origin git@github.com:zhibinz2/Motor_cordination.git

!git push