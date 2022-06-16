function sobi = getsobi(covavg,ncomponents)
% function sobi = getsobi(covavg,pvar)
%
%does the joint diagonalization, returns some quantities of interest. 
% 
%inputs: cavavg - your average cov by lag matrix
%        pvar - % of variance to include in components
%                 enter as a percentage - 95
%
%Returns a structure with the following fields
%eig - average eigenvectors of sobi components
%corr - cross-covariance (bendat correlation) matrix of sobi components 
%power - power associated with each SOBI component
%coh - coherence between sobi components
%pcaeig - pca eigenvectors
%pcacorr - correlation matrix among pca eigenvectors

%if nargin <2;
%    pvar = 1;
%else
%    pvar = pvar/100;
%end;
nlags = size(covavg,3);
nchans = size(covavg,2);
cor0 = squeeze(covavg(:,:,1));
cor0 = cor0./(sqrt(diag(cor0)*diag(cor0)'));
[v,d] = eig(cor0);
d = (diag(d));   %NOTE I DID NOT FLIP HERE.  HAVE TO BE
                     %CAREFUL.  HAS TO DO WITH THE VERY SMALL
                     %EIGENVALUES AND DEGENERACY
[dsort,sortorder] = sort(d,'descend');
psum = cumsum(abs(dsort))/sum(abs(dsort));
%ncomponents = sum(psum < pvar)+1;
d = d(sortorder);
v = v(:,sortorder);
pcacorr = d;
pcaeig = v;
v(:,1:ncomponents) = real((v(:,1:ncomponents)))./ ...
    (ones(nchans,1)*sqrt(d(1:ncomponents))');  %PREWHITENING

for k = 1:nlags;    %rotate each covariance matrix by the
                           %PCA 
pcacorrav(:,:,k) = v(:,1:ncomponents)'*squeeze(covavg(:,:,k))*v(:,1:ncomponents);
end;
disp('I am joint diagonalizing')
[vSOBI,dSOBI] = jointdiag(pcacorrav,1e-4);
vfinal = (ones(nchans,1)*sqrt(d(1:ncomponents))').*v(:,1:ncomponents)*vSOBI;
sobi.eig = vfinal;
sobi.corr = dSOBI;
sobi.pca.eig = pcaeig;
sobi.pca.corr = pcacorrav;
sobi.pca.val = d;