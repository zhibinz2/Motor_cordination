
load('chanlocs.mat')
chanlocs=chanlocs(1:128); 

X = [chanlocs.X];
Y = [chanlocs.Y];
Z = [chanlocs.Z];



%%

[azimuth,elevation,r] = cart2sph(X,Y,Z);

[TH,PHI,R] = cart2sph(X,Y,Z);

[theta,PHI,R] = cart2sph(X,Y,Z);

openExample('symbolic/PlottingInSphericalCoordinateSystemExample')