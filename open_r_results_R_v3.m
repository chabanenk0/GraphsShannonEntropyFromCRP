
inlogfile='dax_04_logsR_v3_wind_500.txt';
initialfile='dax_04.txt';
wind=500;
tstep=10;
% Import the file
newData1 = importdata(inlogfile);

seriesname=strrep(initialfile,'_',' ');
seriesname=strrep(seriesname,'.txt','');

data=newData1.data;
y_init=dlmread(initialfile);
x_init=1:length(y_init);

y_mera=data(:,2);
x_mera=data(:,1)+wind;%wind:tstep:length(y_init);
meraname='topologicalEntr';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,3);
meraname='bertz';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,4);
meraname='radialCentric';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,5);
meraname='graphVertexComplexity';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,6);
meraname='edgeEqualityMIC';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);


y_mera=data(:,7);
meraname='edgeMagnitudeMIC';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,8);
meraname='symmetryIndex';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,9);
meraname='distanceDegreeMIC';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,10);
meraname='distanceDegreeEquality';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);


