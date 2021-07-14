
%SOIdir = uigetdir;
SOIdir = '/Users/idse/Dropbox/flows/flows_shared/codes/ImSAnE_NatureMethods/examples/wingDisc/discProperApicalSOI';
%SOIdir = '/Users/idse/repos/ImSAnE/examples/wingDisc timeseries/discProperApicalSOI';
%%
SOI2 = surfaceAnalysis.SurfaceOfInterest(SOIdir);

%%
surf(double(SOI2.embedding(3).patches{1}.apply{3}))
shading interp;
axis equal;