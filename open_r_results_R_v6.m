
inlogfile='dax_04_logsR_v6_wind_500_globonly.txt';
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


%1 energy
%2 laplacianEnergy
%3 estrada
%4 laplacianEstrada
%5spectralRadius
%6 toString(subgr1),
%7 toString(subgr2),
%8 toString(onedDel1),
%9 toString(onedDel2),
%10 toString(twoeDel),
%11 localClusteringCoeff
%12 globalClusteringCoeff

% Global version (without local clustering) локальный удален вручную
%2energy	3laplacianEnergy	4estrada	5laplacianEstrada	6spectralRadius	7globalClusteringCoeff

y_mera=data(:,2);
x_mera=data(:,1)+wind;%wind:tstep:length(y_init);
meraname='energy';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,3);
meraname='laplacianEnergy';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,4);
meraname='estrada';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,5);
meraname='laplacianEstrada';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,6);
meraname='spectralRadius';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

% y_mera=data(:,6);
% meraname='subgr1';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);
% 
% 
% y_mera=data(:,7);
% meraname='subgr2';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);
% 
% y_mera=data(:,8);
% meraname='onedDel1';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);
% 
% y_mera=data(:,9);
% meraname='onedDel2';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);
% 
% y_mera=data(:,10);
% meraname='twoeDel';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

%11 localClusteringCoeff
%12 globalClusteringCoeff

% y_mera=data(:,11);
% meraname='localClusteringCoeff';
% plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

y_mera=data(:,7);
meraname='globalClusteringCoeff';
plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname);

