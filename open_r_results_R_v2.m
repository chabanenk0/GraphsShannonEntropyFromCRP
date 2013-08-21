
inlogfile='ux_04_crp_logsR_2.txt';
initialfile='ux_04.txt';
wind=500;
tstep=10;
% Import the file
newData1 = importdata(inlogfile);

seriesname=strrep(initialfile,'_',' ');
seriesname=strrep(seriesname,'.txt','');

data=newData1.data;
y_init=dlmread(initialfile);
x_init=1:length(y_init);
%GraphEntropy	WienerIndex	HararyIndex	Meandistvertdev	EccentricGraph	IndexofTotalAdjacency	ZagrebGroupIndices	ZagrebGroupIndices	ZagrebGroupIndices	ZagrebGroupIndices	RandicConnectivityIndex	TheComplexityIndexB	NormalizedEdgeComplexity	Atom-bondConnectivity	Geometric-arithmeticIndices	Geometric-arithmeticIndices	Geometric-arithmeticIndices	Narumi-KatayamaIndex
y_mera=data(:,2);
x_mera=data(:,1)+wind;%wind:tstep:length(y_init);
meraname='GraphEntropy';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,3);
meraname='WienerInder';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,4);
meraname='HararyIndex';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,5);
meraname='Meandistvertdev';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,6);
meraname='EccentricGraph';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,7);
meraname='IndexofTotalAdjacency';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,8);
meraname='ZagrebGroupIndices';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,9);
meraname='ZagrebGroupIndices2';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,10);
meraname='ZagrebGroupIndices3';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,11);
meraname='ZagrebGroupIndices4';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,12);
meraname='RandicConnectivityIndex';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,13);
meraname='ZagrebGroupIndices2';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,14);
meraname='ZagrebGroupIndices3';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,15);
meraname='ZagrebGroupIndices4';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);


y_mera=data(:,16);
meraname='TheComplexityIndexB';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,17);
meraname='NormalizedEdgeComplexity';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,18);
meraname='AtombondConnectivity';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,16);
meraname='GeometricarithmeticIndices1';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,17);
meraname='GeometricarithmeticIndices2';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,18);
meraname='GeometricarithmeticIndices3';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,19);
meraname='NarumiKatayamaIndex';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);


