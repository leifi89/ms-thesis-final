
%% Start file by clearing all variables and configure the home directory
clear all; clear session; close all; clc
warning off all

cd ~/Desktop/'MS Thesis'/ms-thesis-final/MatlabCode/SRLeifi/

% Greining 1 tekur svo langan tíma að keyra 
keyraGreiningu1 = false; % Hér er hægt að stilla hvort það eigi að keyra hana eða ekki 



%%  Set path prefixes

versionFolderName = '20240805_Greining2/'

myFigPath2  = strcat('../Exported/', versionFolderName, 'ExportedFigures/')
myFigPath  = strcat('../', myFigPath2)
myTexPath  = strcat('../Exported/', versionFolderName , 'ExportedTex/' )
myDataPath = strcat('../Exported/', versionFolderName , 'Data/' )


%% Load data from file 
load(strcat('../Import/', 'FredData'))

checkForNewData = true;
if checkForNewData
   
url = 'https://fred.stlouisfed.org/';
%url = 'https://research.stlouisfed.org/fred2/';
c = fred(url);
rgdpISFetch = fetch(c,'CLVMNACSAB1GQIS');
gdpdefISFetch = fetch(c,'ISLGDPDEFQISMEI');
cpiISFetch = fetch(c,'CPALTT01ISQ657N');
close(c)

fredHasChanged = 0;
if( strcmp(rgdpISFetch.LastUpdated ,  rgdpIS.LastUpdated)==0 ) disp('rgdpIS data has changed!!'); fredHasChanged = fredHasChanged + 1; end
if( strcmp(gdpdefISFetch.LastUpdated ,  gdpdefIS.LastUpdated)==0 ) disp('gdpdefIS data has changed!!'); fredHasChanged= fredHasChanged + 1; end
if(fredHasChanged > 0 ) return; end 

end

% Update saved data if changed-  Uncomment if we want to update source data 
%rgdpIS = rgdpISFetch;
%gdpdefIS = gdpdefISFetch;
%save( strcat(myDataPath, '20240609/', 'FredData'), "gdpdefIS", "rgdpIS" )



%% 
varNames = ["VLF", "VVLF"] 
minx = min(size(rgdpIS.Data,1),size(gdpdefIS.Data,1));
fredCombined = [ rgdpIS.Data(1:minx,1) rgdpIS.Data(1:minx,2) gdpdefIS.Data(1:minx,2)];
t = fredCombined(:, 1);
Xraw = fredCombined(:, 2:end);
XrawTbl = array2table(Xraw,"VariableNames",varNames);

vnames = {'VLF', 'VVLF'}; 
vnames_long = vnames;

Xlog = log(Xraw);
X = diff(Xlog);
XTbl = array2table(X,"VariableNames",varNames);


%% NEY fra hagstofu

xlsDataCPITbl = readtable('../Import/VIS01000_20240701-231606.xlsx','Sheet','ImportDataOriginal');
xlsDataHagvxtTbl = readtable('../Import/THJ01000_20240705-233617.xlsx','Sheet','ImportData');
dthag = xlsDataHagvxtTbl.Var3
xlsDataMeginVextirTbl = readtable('../Import/Meginvextir.xlsx','Sheet','ImportData');
xlsDataHagstVvlfTbl = readtable('../Import/THJ01104_20240706-225820.xlsx','Sheet','ImportData');

figure
dt2 = xlsDataCPITbl.Dags(21:end);
arsbr = xlsDataCPITbl.Arsbreyting(21:end)
mean(arsbr)
hold on
plot(dt2, arsbr)
plot(dthag, xlsDataHagvxtTbl.Var1)
plot(xlsDataMeginVextirTbl.Var1, xlsDataMeginVextirTbl.Var2)
plot(xlsDataHagstVvlfTbl.Var1, xlsDataHagstVvlfTbl.Var3*100)

legend(["NEY", "Hagvöxtur", "Meginvextir SÍ", "VVLF"],'Orientation','Horizontal', 'Location', 'southoutside', 'Box','off', 'Fontsize', get(gca,'Fontsize')-1,'Interpreter', 'Tex')
axis tight
set(gcf, 'Color', 'w');
xtickangle(90)
ytickformat('percentage')
xlim([datetime(1994,1,1) datetime(2024,01,01)])
export_fig(strcat( myFigPath,  'CPI_HAGSTOFA_ANNUALIZED'),'-pdf','-painters');

xlim([datetime(2016,1,1) datetime(2024,01,01)])
export_fig(strcat( myFigPath,  'CPI_HAGSTOFA_ANNUALIZED_COVID'),'-pdf','-painters');


%% Teikna grunngögn 

for iplot=1:2
%close all
figure

plotYTickFormat = NaN;
if iplot==1
    plotX = Xraw;
    plotX(:,1) = plotX(:,1)/1e3;
    plotFigurePath = strcat( myFigPath,  'Grunngogn_FRED');
    plotVarNames = varNames;
    plotYTickFormat = '%g ma';
end
if iplot==2
    plotX = Xlog;
    plotFigurePath = strcat( myFigPath, 'Grunngogn_FRED_log');
    plotVarNames = varNames; %["Real GDP Log", "GDP Deflator Log"];
    %title('Logri af upprunagögnum')
end
if iplot==3
    plotX = X;
    plotFigurePath =  strcat( myFigPath, 'Grunngogn_FRED_logdiff');
    plotVarNames =  varNames; 
    %title('Fyrstu gráðu mismunur og logri af upprunagögnum')
end

xtix = t(1:4:end,1) ;
xlab=datestr(xtix, 'YYYY-QQ');  
pl = plot(t, plotX(:,1), "Color",cmap(1)); 
xlim([t(1) t(end)])
datetick('x','YYYY-QQ', 'keeplimits')
set(gca,'xtick',xtix)                     
set(gca,'xticklabel',xlab)                
set(gca,'xminortick','on')
if not(isnan(plotYTickFormat))  ytickformat(plotYTickFormat); end
xtickangle(90)

yyaxis right
h = plot(t, plotX(:,2), "Color",cmap(2));
ax = gca;
ax.YAxis(2).Color = 'k';
set(gca,'Layer','top');
legend(plotVarNames,'Orientation','Horizontal', 'Location', 'southoutside', 'Box','off', 'Fontsize', get(gca,'Fontsize')-1,'Interpreter', 'Tex')
%title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
set(gcf, 'Color', 'w');

yyaxis left 
ylabel(plotVarNames(1))
yyaxis right
ylabel(plotVarNames(2))

export_fig(plotFigurePath,'-pdf','-painters');


end


%% Plot Differencing
%close all
figure; hold on;

plotX = X;
plotX  = plotX * 100; ytickformat('percentage');
pl = plot(t(1:end-1), plotX(:,1), "Color",cmap(1)); 
plot(t(1:end-1), plotX(:,2), "Color",cmap(2));
xlim([t(1) t(end-1)])
datetick('x','YYYY-QQ', 'keeplimits')
set(gca,'xtick',xtix)                     
set(gca,'xticklabel',xlab)                
set(gca,'xminortick','on')
xtickangle(90)
set(gca,'Layer','top');
legend(["VLF" "VVLF"],'Orientation','Horizontal', 'Location', 'southoutside', 'Box','off', 'Fontsize', get(gca,'Fontsize')-1,'Interpreter', 'Tex')
set(gcf, 'Color', 'w');

export_fig(strcat(myFigPath, '/Grunngogn_FRED_log_diff'),'-pdf','-painters');
hold off

%%

adfkpssTestArray = [];

for ivar= 1:length(varNames)
iVarNam = strcat('',  varNames(ivar) );
iVarNam = convertStringsToChars(iVarNam);

fprintf('\n==== Augmented Dickey–Fuller -  %s (Original) ==== \n', varNames(ivar));
adftestoriginal =  adftest(XrawTbl, 'DataVariable', varNames(ivar), 'lags',0:8);
adftestoriginal = addvars(adftestoriginal,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','ADF Variable Name');
table2latex(adftestoriginal , strcat( myTexPath, 'ADFTestOriginal_', iVarNam, '.tex' ))

fprintf('\n==== Augmented Dickey–Fuller -  %s (Transformed) ==== \n', varNames(ivar));
adftestWithDiff = adftest(XTbl, 'DataVariable', varNames(ivar), 'lags',0:8);
adftestWithDiff = addvars(adftestWithDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','ADF Variable Name');
table2latex(adftestWithDiff , strcat( myTexPath, 'ADFTestWithDifference_', iVarNam, '.tex' ))

fprintf('\n==== Kwiatkowski–Phillips–Schmidt–Shin (KPSS) tests - %s (Original) ==== \n', varNames(ivar));
kpsstestoriginal = kpsstest(XrawTbl, 'DataVariable', varNames(ivar), 'lags',0:8,'trend',true );
kpsstestoriginal = addvars(kpsstestoriginal,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','KPSS Test Variable Name');
table2latex(kpsstestoriginal , strcat( myTexPath, 'KPSSTestOriginal_', iVarNam, '.tex' ))

fprintf('\n==== Kwiatkowski–Phillips–Schmidt–Shin (KPSS) tests - %s (Transformed) ==== \n', varNames(ivar));
kpsstestWithDiff  = kpsstest(XTbl, 'DataVariable', varNames(ivar), 'lags',0:8,'trend',true );
kpsstestWithDiff = addvars(kpsstestWithDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','KPSS Test Variable Name');
table2latex(kpsstestWithDiff , strcat( myTexPath, 'KPSSTestWithDifference_', iVarNam, '.tex' ))

adfkpssTestArray = [adfkpssTestArray [adftestoriginal.h adftestWithDiff.h kpsstestoriginal.h kpsstestWithDiff.h]]
table2latex( array2table(adfkpssTestArray) , strcat( myTexPath, 'ADFandKPSSTestH', '.tex' ))


end

%% Lag order - Velja hvaða lag á að nota
det = 1;

[lagsVec, AIC, BIC, logL, HQIC] = VARlag_LeifiCustom(X, 12,1);
LagOrderArray = [lagsVec, AIC, BIC, HQIC];
LagOrderTbl =  array2table(LagOrderArray );
LagOrderTbl.Properties.VariableNames(1:end) = {'Lag Order', 'AIC','BIC','HQIC'};
table2latex(LagOrderTbl, strcat(myTexPath, 'LagOrder.tex'))

AIC_LagOrder = find(AIC==min(AIC)); % 4
BIC_LagOrder = find(BIC==min(BIC)); % 1
HQIC_LagOrder = find(HQIC==min(HQIC)); %2


%% Ýmis próf eftir að lag hefur verið valid (úr prófi)

[VAR_AIC, VARopt_AIC] = VARmodel(X,AIC_LagOrder,det);
[VAR_BIC, VARopt_BIC] = VARmodel(X,BIC_LagOrder,det);
[VAR_HQIC, VARopt_HQIC] = VARmodel(X,HQIC_LagOrder,det);

%%

lagOrderRow = [AIC_LagOrder, BIC_LagOrder, HQIC_LagOrder];

% Engle
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = archtest(VAR_AIC.resid(: ,1), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = archtest(VAR_BIC.resid(: ,1), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = archtest(VAR_HQIC.resid(:,1), Lags=1:8); 
engleTest1_h = [h_AIC' h_BIC' h_HQIC'];
engleTest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = archtest(VAR_AIC.resid(: ,2), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = archtest(VAR_BIC.resid(: ,2), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = archtest(VAR_HQIC.resid(:,2), Lags=1:8); 
engleTest2_h = [h_AIC' h_BIC' h_HQIC'];
engleTest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

% Ljung-Box Q-Test
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = lbqtest(VAR_AIC.resid(: ,1), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = lbqtest(VAR_BIC.resid(: ,1), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = lbqtest(VAR_HQIC.resid(:,1), Lags=1:8); 
lbqtest1_h = [h_AIC' h_BIC' h_HQIC'];
lbqtest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = lbqtest(VAR_AIC.resid(: ,2), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = lbqtest(VAR_BIC.resid(: ,2), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = lbqtest(VAR_HQIC.resid(:,2), Lags=1:8); 
lbqtest2_h = [h_AIC' h_BIC' h_HQIC'];
lbqtest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

% One-sample Kolmogorov-Smirnov test
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = kstest(VAR_AIC.resid(: ,1)); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = kstest(VAR_BIC.resid(: ,1)); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = kstest(VAR_HQIC.resid(:,1)); 
kstest1_h = [h_AIC' h_BIC' h_HQIC'];
kstest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = kstest(VAR_AIC.resid(: ,2)); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = kstest(VAR_BIC.resid(: ,2)); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = kstest(VAR_HQIC.resid(:,2)); 
kstest2_h = [h_AIC' h_BIC' h_HQIC'];
kstest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

% [h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = jbtest(VAR_AIC.resid(: ,1)); 
% [h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = jbtest(VAR_BIC.resid(: ,1)); 
% [h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = jbtest(VAR_HQIC.resid(:,1)); 
% jbtest1_h = [h_AIC' h_BIC' h_HQIC'];
% jbtest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];
% [h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = jbtest(VAR_AIC.resid(: ,2)); 
% [h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = jbtest(VAR_BIC.resid(: ,2)); 
% [h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = jbtest(VAR_HQIC.resid(:,2)); 
% jbtest2_h = [h_AIC' h_BIC' h_HQIC'];
% jbtest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];
% 
% [h_AIC, pValue_AIC]     = LeifiCustom_swtest(VAR_AIC.resid(: ,1)); 
% [h_BIC, pValue_BIC]     = LeifiCustom_swtest(VAR_BIC.resid(: ,1)); 
% [h_HQIC, pValue_HQIC] = LeifiCustom_swtest(VAR_HQIC.resid(:,1)); 
% swtest1_h = [h_AIC' h_BIC' h_HQIC'];
% swtest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];
% [h_AIC, pValue_AIC]     = LeifiCustom_swtest(VAR_AIC.resid(: ,2)); 
% [h_BIC, pValue_BIC]     = LeifiCustom_swtest(VAR_BIC.resid(: ,2)); 
% [h_HQIC, pValue_HQIC] = LeifiCustom_swtest(VAR_HQIC.resid(:,2)); 
% swtest2_h = [h_AIC' h_BIC' h_HQIC'];
% swtest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];


% Multivariate ljung-box Q-test
%[h,pValue,stat,cValue] = LeifiCustom_mlbqtest(resid1,[3,6,10]); 
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC, dof_AIC]     = LeifiCustom_mlbqtest(VAR_AIC.resid    , 1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC, dof_BIC]     = LeifiCustom_mlbqtest(VAR_BIC.resid    , 1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC, dof_HQIC] = LeifiCustom_mlbqtest(VAR_HQIC.resid , 1:8); 
mlbqtest_h = [h_AIC' h_BIC' h_HQIC'];
mlbqtest_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];
mlbqtest_dof = [dof_AIC' dof_BIC' dof_HQIC'];
mlbqtestDesc =  'Multivariate Portmanteau test with ';

finalLagOrderArray = [lagOrderRow;
    engleTest1_h; engleTest1_pValue ; engleTest2_h; engleTest2_pValue ;
    lbqtest1_h; lbqtest1_pValue;  lbqtest2_h; lbqtest2_pValue; 
    kstest1_h;kstest1_pValue;kstest2_h;kstest2_pValue;
    mlbqtest_h;mlbqtest_pValue; mlbqtest_dof;
];

LagOrderSelectionTbl = array2table(finalLagOrderArray);
LagOrderSelectionTbl.Properties.VariableNames(1:end) = {'AIC','BIC','HQIC'};

LagOrderSelectionTbl.Properties.RowNames = {'Lag Order' ...
    'Engle Test (GDP) h Lags 1' ...
    'Engle Test (GDP) h Lags 2' ...
    'Engle Test (GDP) h Lags 3' ...
    'Engle Test (GDP) h Lags 4' ...
    'Engle Test (GDP) h Lags 5' ...
    'Engle Test (GDP) h Lags 6' ...
    'Engle Test (GDP) h Lags 7' ...
    'Engle Test (GDP) h Lags 8' ...
    'Engle Test (GDP) p-value Lags 1' ...
    'Engle Test (GDP) p-value Lags 2' ...
    'Engle Test (GDP) p-value Lags 3' ...
    'Engle Test (GDP) p-value Lags 4' ...
    'Engle Test (GDP) p-value Lags 5' ...
    'Engle Test (GDP) p-value Lags 6' ...
    'Engle Test (GDP) p-value Lags 7' ...
    'Engle Test (GDP) p-value Lags 8' ...
    'Engle Test (GDPDEF) h Lags 1' ...
    'Engle Test (GDPDEF) h Lags 2' ...
    'Engle Test (GDPDEF) h Lags 3' ...
    'Engle Test (GDPDEF) h Lags 4' ...
    'Engle Test (GDPDEF) h Lags 5' ...
    'Engle Test (GDPDEF) h Lags 6' ...
    'Engle Test (GDPDEF) h Lags 7' ...
    'Engle Test (GDPDEF) h Lags 8' ...
    'Engle Test (GDPDEF) p-value Lags 1' ...
    'Engle Test (GDPDEF) p-value Lags 2' ...
    'Engle Test (GDPDEF) p-value Lags 3' ...
    'Engle Test (GDPDEF) p-value Lags 4' ...
    'Engle Test (GDPDEF) p-value Lags 5' ...
    'Engle Test (GDPDEF) p-value Lags 6' ...
    'Engle Test (GDPDEF) p-value Lags 7' ...
    'Engle Test (GDPDEF) p-value Lags 8' ...
    ...
    'Ljung-Box Q-Test (GDP) h Lags 1' ...
    'Ljung-Box Q-Test (GDP) h Lags 2' ...
    'Ljung-Box Q-Test (GDP) h Lags 3' ...
    'Ljung-Box Q-Test (GDP) h Lags 4' ...
    'Ljung-Box Q-Test (GDP) h Lags 5' ...
    'Ljung-Box Q-Test (GDP) h Lags 6' ...
    'Ljung-Box Q-Test (GDP) h Lags 7' ...
    'Ljung-Box Q-Test (GDP) h Lags 8' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 1' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 2' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 3' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 4' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 5' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 6' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 7' ...
    'Ljung-Box Q-Test (GDP) p-value Lags 8' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 1' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 2' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 3' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 4' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 5' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 6' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 7' ...
    'Ljung-Box Q-Test (GDPDEF) h Lags 8' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 1' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 2' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 3' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 4' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 5' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 6' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 7' ...
    'Ljung-Box Q-Test (GDPDEF) p-value Lags 8' ...
    ...
'Kolmogorov-Smirnov Test (GDP) h' 'Kolmogorov-Smirnov Test (GDP) P-Value' 'Kolmogorov-Smirnov Test (GDPDEF) h' 'Kolmogorov-Smirnov Test (GDPDEF) P-Value' ...
...
...
    'Multivariate Portmanteau test h Lags 1' ...
    'Multivariate Portmanteau test h Lags 2' ...
    'Multivariate Portmanteau test h Lags 3' ...
    'Multivariate Portmanteau test h Lags 4' ...
    'Multivariate Portmanteau test h Lags 5' ...
    'Multivariate Portmanteau test h Lags 6' ...
    'Multivariate Portmanteau test h Lags 7' ...
    'Multivariate Portmanteau test h Lags 8' ...
    'Multivariate Portmanteau test p-value Lags 1' ...
    'Multivariate Portmanteau test p-value Lags 2' ...
    'Multivariate Portmanteau test p-value Lags 3' ...
    'Multivariate Portmanteau test p-value Lags 4' ...
    'Multivariate Portmanteau test p-value Lags 5' ...
    'Multivariate Portmanteau test p-value Lags 6' ...
    'Multivariate Portmanteau test p-value Lags 7' ...
    'Multivariate Portmanteau test p-value Lags 8' ...
    'Multivariate Portmanteau test DF Lags 1' ...
    'Multivariate Portmanteau test DF Lags 2' ...
    'Multivariate Portmanteau test DF Lags 3' ...
    'Multivariate Portmanteau test DF Lags 4' ...
    'Multivariate Portmanteau test DF Lags 5' ...
    'Multivariate Portmanteau test DF Lags 6' ...
    'Multivariate Portmanteau test DF Lags 7' ...
    'Multivariate Portmanteau test DF Lags 8' ...
    }

table2latex(LagOrderSelectionTbl, strcat(myTexPath, 'LagOrderSelection.tex'));

finalLagOrderArrayShortHValues =  [engleTest1_h;  engleTest2_h;
    lbqtest1_h;  lbqtest2_h; 
    kstest1_h;kstest2_h;
    mlbqtest_h;]
LagOrderSelectionHValuesTbl = array2table(finalLagOrderArrayShortHValues);
LagOrderSelectionHValuesTbl.Properties.VariableNames(1:end) = {'AIC','BIC','HQIC'};

LagOrderSelectionHValuesTbl.Properties.RowNames = {  ...
    'Engle próf (GDP)          Lags 1' ...
    'Engle próf (GDP)          Lags 2' ...
    'Engle próf (GDP)          Lags 3' ...
    'Engle próf (GDP)          Lags 4' ...
    'Engle próf (GDP)          Lags 5' ...
    'Engle próf (GDP)          Lags 6' ...
    'Engle próf (GDP)          Lags 7' ...
    'Engle próf (GDP)          Lags 8' ...
    'Engle próf (GDPDEF)       Lags 1' ...
    'Engle próf (GDPDEF)       Lags 2' ...
    'Engle próf (GDPDEF)       Lags 3' ...
    'Engle próf (GDPDEF)       Lags 4' ...
    'Engle próf (GDPDEF)       Lags 5' ...
    'Engle próf (GDPDEF)       Lags 6' ...
    'Engle próf (GDPDEF)       Lags 7' ...
    'Engle próf (GDPDEF)       Lags 8' ...
    ...
    'Ljung-Box Q-próf (GDP)    Lags 1' ...
    'Ljung-Box Q-próf (GDP)    Lags 2' ...
    'Ljung-Box Q-próf (GDP)    Lags 3' ...
    'Ljung-Box Q-próf (GDP)    Lags 4' ...
    'Ljung-Box Q-próf (GDP)    Lags 5' ...
    'Ljung-Box Q-próf (GDP)    Lags 6' ...
    'Ljung-Box Q-próf (GDP)    Lags 7' ...
    'Ljung-Box Q-próf (GDP)    Lags 8' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 1' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 2' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 3' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 4' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 5' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 6' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 7' ...
    'Ljung-Box Q-próf (GDPDEF) Lags 8' ...
    ...
'Kolmogorov-Smirnov Test (GDP)'          ...
'Kolmogorov-Smirnov Test (GDPDEF)'       ...  
...
...
    'Fjölvítt Portmanteau próf Lags 1' ...
    'Fjölvítt Portmanteau próf Lags 2' ...
    'Fjölvítt Portmanteau próf Lags 3' ...
    'Fjölvítt Portmanteau próf Lags 4' ...
    'Fjölvítt Portmanteau próf Lags 5' ...
    'Fjölvítt Portmanteau próf Lags 6' ...
    'Fjölvítt Portmanteau próf Lags 7' ...
    'Fjölvítt Portmanteau próf Lags 8' ...
    }

table2latex(LagOrderSelectionHValuesTbl, strcat(myTexPath, 'LagOrderSelectionHValues.tex'));




%% Plot test for residuals

for iac = 1:3

if (iac==1)  
    tacVar = VAR_AIC;
    tacStr = 'AIC';
elseif(iac==2)
    tacVar = VAR_BIC;
     tacStr = 'BIC';
elseif(iac==3)
    tacVar = VAR_HQIC;
     tacStr = 'HQIC';
end

figure
set(gcf, 'Color', 'w');
subplot(3,1,1)
tempResid = tacVar.resid(:,1);
[acf,lags,bounds] = autocorr(tempResid);
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(vnames_long(1))

subplot(3,1,2)
tempResid = tacVar.resid(:,2);
[acf,lags,bounds] = autocorr(tempResid);
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(vnames_long(2))

subplot(3,1,3)
[acf,lags,bounds] = crosscorr(tacVar.resid(:,1), tacVar.resid(:,2));
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title('Cross correlation')
%set(gca,'Layer','top');

export_fig(strcat(myFigPath, 'Autocorrelation_' , tacStr),'-pdf','-painters');

figure
set(gcf, 'Color', 'w');
qqplot(tacVar.resid)
%ylabel('')
%xlabel('')
title('')
export_fig(strcat(myFigPath, 'QQPlot' , tacStr),'-pdf','-painters');

figure
set(gcf, 'Color', 'w');
subplot(2,1,1)
histogram(tacVar.resid(:,1))
subplot(2,1,2)
histogram(tacVar.resid(:,2))


export_fig(strcat(myFigPath, 'Histogram' , tacStr),'-pdf','-painters');


figure
hold on
set(gcf, 'Color', 'w');
plot(tacVar.resid(:,1), 'x')
plot(tacVar.resid(:,2), 'o')
export_fig(strcat(myFigPath, 'ScatterOfResiduals' , tacStr),'-pdf','-painters');


end


%%
% VAR ESTIMATION Test 1
% =====================

[nobs, nvar] = size(X);
det = 1;
nlags = 4;  
[VAR, VARopt] = VARmodel(X,nlags,det);
VARopt.vnames = vnames_long; %Variable names
[beta, TABLE] = VARprint(VAR,VARopt);

varRFETable = array2table( TABLE(2:end, 2:end));
varRFETable.Properties.VariableNames =  TABLE(1, 2:end);
varRFETable.Properties.RowNames = TABLE(2:end, 1);
    table2latex(varRFETable, strcat(myTexPath, 'varRFETable.tex'));
varSigmaTable = array2table(VAR.sigma);
    table2latex(varSigmaTable, strcat(myTexPath, 'varSigmaTable.tex'));
varBetaTable = array2table( beta(2:end, 2:end));
varBetaTable.Properties.VariableNames =  beta(1, 2:end);
varBetaTable.Properties.RowNames = beta(2:end, 1);
    table2latex(varBetaTable, strcat(myTexPath, 'varBetaTable.tex'))

(datetime(fredCombined(1,1),'ConvertFrom','datenum','Format','yyyy-QQ'))
%initialDateNumberValue = fredCombined(1,1);
%initialDate = str2double(string(datetime(initialDateNumberValue,'ConvertFrom','datenum','Format','yyyy')))*1.00  + (str2double(string(datetime(initialDateNumberValue,'ConvertFrom','datenum','Format','QQ')))-1)/4.000


initialDate = 1995.00+nlags*0.25; 

VARopt.firstdate = initialDate;      % initial date of the sample in format 1999.75 => 1999Q4 (both for annual and quarterly data)
VARopt.frequency = 'q'; 
VARopt.snames = {'Framboð', 'Eftirspurn'}; %{'Supply Shock', 'Demand Shock'}; % Shock names

VARopt.nsteps = 12; %nsteps;
VARopt.ident   = 'sign'; 
VARopt.sr_mod  = 1; % 0 no Bayesian draw / 1 Bayesian
VARopt.sr_hor = 2;  % 
VARopt.ndraws  = 10; 

table2latex(array2table( [VARopt.ndraws VARopt.nsteps VARopt.sr_hor ], "VariableNames", {'Number of Draws' 'Number of Steps' 'Horizon'}), strcat(myTexPath, '_SD_IS_Parameters'))

R = [ 1   1       % gdp         %Quantity
     -1   1 ];    % gdp_def     %Price

VARopt.sr_draw = 100000*1000;
VARopt.sr_rot = 500*2; %Hámarksfjöldi snúninga



if(keyraGreiningu1)
disp('Begin SRout=SR(...)')
temCurrTime = datetime;
disp(temCurrTime)

[SRout , VARDraws ]= SR(VAR,R,VARopt); % 22 klst

disp(datetime)
disp(datetime-temCurrTime)
disp('End SRout=SR(...)')

pFigName = strcat(myFigPath, 'SD_IS_nsteps', num2str(VARopt.nsteps) , '_hor', num2str(VARopt.sr_hor), '_ndraws', num2str(VARopt.ndraws));
VARopt.figname = pFigName;
VARopt.suptitle = 0;
VARopt.datestype = 1;       % 1 smart labels; 2 less smart labels

end

%% HD plot
%VARhdplot_LeifiCustom(SRout.HD, VARopt) 

%DatesPlot(VARopt.firstdate,nsteps,8*7,VARopt.frequency);
%xlim([4 8*14.5])

if(exist('SRout'))

disp('SRout exists');

HD = SRout.HD;

[nsteps, nvars, nshocks] = size(HD.shock);
for ii=1:nvars
    figure; FigSize(27,15); clf;
    %H = AreaPlot(squeeze(HD.shock(:,:,ii))); hold on; 
    H = BarPlot(squeeze(HD.shock(:,:,ii))); hold on; 
%    h = plot(fredCombined(1:end-1, 1), sum(squeeze(HD.shock(:,:,ii)),2),'-k','LineWidth',2);
    h = plot( sum(squeeze(HD.shock(:,:,ii)),2),'-k','LineWidth',2);
%    if ~isempty(VARopt.firstdate); DatesPlot(VARopt.firstdate,nsteps,8,VARopt.frequency); end    
    xlim([1 nsteps]);
	set(gca,'Layer','top');
    title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
    % Save
       if 0==1 %suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            tempSubtitle = convertCharsToStrings( [Alphabet(ii) ') HD of '  vnames{ii}]);
            SupTitle(tempSubtitle)
        end
        opt = LegOption; opt.handle = [H(1,:) h];
        LegSubplot([VARopt.snames {'Gögn'}],opt);
        set(gcf, 'Color', 'w');
        
%dateVecOut = DatesPlot(VARopt.firstdate,nsteps,8*4,VARopt.frequency);
xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
%tstr = (string(t, 'yyyy-QQ'));
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)

xtickangle(90)
xlim([4 116])
%set(gca,'xTick',tick:tick:tick_max,'xTickLabel',labyear(tick:tick:tick_max),'xLim',[0 nobs+1],'Layer','top');
%set(gca,'xTick',4:114)%,'xTickLabel',labyear(tick:tick:tick_max),'xLim',[0 nobs+1],'Layer','top');
savefig(strcat(myFigPath2, '_HD_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD_', num2str(ii)),'-pdf','-painters');                 

xlim([100 114])
savefig(strcat(myFigPath2, '_HD_Zoom_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD_Zoom_', num2str(ii)),'-pdf','-painters');

figure; FigSize(27,15); clf;
H = AreaPlot(squeeze(HD.shock(:,:,ii))); hold on; 
h = plot( sum(squeeze(HD.shock(:,:,ii)),2),'-k','LineWidth',2);
set(gca,'Layer','top');
title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
opt = LegOption; opt.handle = [H(1,:) h];
LegSubplot([VARopt.snames {'Gögn'}],opt);
set(gcf, 'Color', 'w');
%dateVecOut = DatesPlot(VARopt.firstdate,nsteps,8*4,VARopt.frequency);
xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
%tstr = (string(t, 'yyyy-QQ'));
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)

xtickangle(90)
xlim([4 116])
savefig(strcat(myFigPath2, '_HD_Area_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD_Area_', num2str(ii)),'-pdf','-painters');                 

xlim([100 114])
savefig(strcat(myFigPath2, '_HD_Area_Zoom_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD_Area_Zoom_', num2str(ii)),'-pdf','-painters');


temp1 = squeeze(HD.shock(:,:,ii)); 
temp1Sum = sum(temp1, 2);
temp1SumAbs = sum(abs(temp1), 2);
temp1perc = abs(temp1)./temp1SumAbs;
temp1percWithSign = temp1perc.* sign(temp1);
tempSignalPercOfTotalPower = temp1Sum./sum(abs(temp1), 2);
tempXlsHDData = [(1995.0:0.25:1995+0.25*113)' temp1 temp1Sum temp1SumAbs temp1perc temp1percWithSign tempSignalPercOfTotalPower];
tempXlsHDDataTbl = array2table(tempXlsHDData, 'VariableNames', {'Date', 'Shock 1', 'Shock 2', 'Sum(Shocks) [i.e.Signal]' , 'Sum(abs(shocks))', 'Percentage Shock 1' , 'Percantage Shock 2' 'Signed Percentage Shock 1' , 'Signed Percantage Shock 2' 'Net Percentage Effect' });
writetable(tempXlsHDDataTbl, strcat(myDataPath, '_HD_', '.xls') , 'Sheet', vnames{ii} );
writetable(array2table(VARopt.snames), strcat(myDataPath, '_HD_', '.xls') , 'Sheet', 'Shocks' );

calcInd1 = 5; calcInd2 = 101;
tempWriteTableRatios = [
                [sum(abs( temp1(calcInd1:end, :) ))/ sum(sum(abs( temp1(calcInd1:end, :) )))       NaN];
                mean(abs( temp1(calcInd1:end, :)) ./    sum(abs( temp1(calcInd1:end, :)), 2) )      mean( abs(sum(temp1(calcInd1:end, :), 2)) ./ sum(abs(temp1(calcInd1:end, :)), 2)) ;

                [sum(abs( temp1(calcInd2:end, :) ))/ sum(sum(abs( temp1(calcInd2:end, :) )))       NaN];
                mean(abs( temp1(calcInd2:end, :)) ./     sum(abs( temp1(calcInd2:end, :)), 2) )      mean( abs(sum(temp1(calcInd2:end, :), 2)) ./ sum(abs(temp1(calcInd2:end, :)), 2)) ;
]; 
table2latex(( array2table(tempWriteTableRatios ,'RowNames', {'Ratio Method 1 1995-' , 'Ratio Method 2 1995-' , 'Ratio Method 1 2000-' , 'Ratio Method 2 2000-'  }, 'VariableNames', [vnames 'Efficiency']) ) , strcat(myTexPath , 'HDRatios', '_HD_Var_', num2str(ii), '.tex'));


end

end


save(strcat(myDataPath, 'MyWorkspace') )


%%
% VAR ESTIMATION Test 2
% ======================

[nobs, nvar] = size(X);
det = 1;
nlags = 4;  
[VAR, VARopt] = VARmodel(X,nlags,det);
VARopt.vnames = vnames_long; 
[beta, TABLE] = VARprint(VAR,VARopt);


VARopt.firstdate = 1995.00+nlags*0.25;      
VARopt.frequency = 'q'; 
VARopt.snames = {'Framboð', 'Eftirspurn'};   


VARopt.nsteps = 40; %nsteps;
VARopt.ident   = 'sign'; 
VARopt.sr_mod  = 1; % 0 no Bayesian draw / 1 Bayesian
VARopt.sr_hor = 1;  % 
VARopt.ndraws  = 5000*10; 


R = [ 1   1       % gdp         %Quantity
     -1   1 ];    % gdp_def     %Price

VARopt.sr_draw = 100000*1000;
VARopt.sr_rot = 500*2; %Hámarksfjöldi snúninga
disp('Begin SRout2=SR(...)')
temCurrTime = datetime;
disp(temCurrTime)
[SRout2 , VARDraws2 ] = SR(VAR,R,VARopt); 
disp(datetime)
disp(datetime-temCurrTime)
disp('End SRout2=SR(...)')

pFigName = strcat(myFigPath, 'SD_IS_nsteps', num2str(VARopt.nsteps) , '_hor', num2str(VARopt.sr_hor), '_ndraws', num2str(VARopt.ndraws));
VARopt.figname = pFigName;
VARopt.suptitle = 0;
VARopt.datestype = 1;       % 1 smart labels; 2 less smart labels


%%

HD2 = SRout2.HD;
hd2snames = {'Framboð', 'Eftirspurn'};

[nsteps, nvars, nshocks] = size(HD2.shock);
for ii=1:nvars
    figure; FigSize(27,15); clf;
    H = BarPlot(squeeze(HD2.shock(:,:,ii))); hold on; 
    h = plot( sum(squeeze(HD2.shock(:,:,ii)),2),'-k','LineWidth',2);
%   if ~isempty(VARopt.firstdate); DatesPlot(VARopt.firstdate,nsteps,8,VARopt.frequency); end    
    xlim([1 nsteps]);
	set(gca,'Layer','top');
    title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
    % Save
       if 0==1 %suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            tempSubtitle = convertCharsToStrings( [Alphabet(ii) ') HD of '  vnames{ii}]);
            SupTitle(tempSubtitle)
        end
        opt = LegOption; opt.handle = [H(1,:) h];
        LegSubplot([hd2snames {'Gögn'}],opt);
        set(gcf, 'Color', 'w');
        
%dateVecOut = DatesPlot(VARopt.firstdate,nsteps,8*4,VARopt.frequency);
xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
%tstr = (string(t, 'yyyy-QQ'));
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)
xtickangle(90)

%set(gca,'xTick',tick:tick:tick_max,'xTickLabel',labyear(tick:tick:tick_max),'xLim',[0 nobs+1],'Layer','top');
%set(gca,'xTick',4:114)%,'xTickLabel',labyear(tick:tick:tick_max),'xLim',[0 nobs+1],'Layer','top');

xlim([4 116])
export_fig(strcat(myFigPath, '_HD2_', num2str(ii)),'-pdf','-painters');                 
savefig(strcat(myFigPath2, '_HD2_', num2str(ii)))


xlim([100 114])
export_fig(strcat(myFigPath, '_HD2_Zoom_', num2str(ii)),'-pdf','-painters');
savefig(strcat(myFigPath2, '_HD2_Zoom_', num2str(ii)))

clf
% figure
H = AreaPlot(squeeze(HD2.shock(:,:,ii))); hold on; 
h = plot( sum(squeeze(HD2.shock(:,:,ii)),2),'-k','LineWidth',2);
set(gca,'Layer','top');
title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
opt = LegOption; opt.handle = [H(1,:) h];
LegSubplot([hd2snames {'Gögn'}],opt);
set(gcf, 'Color', 'w');
%dateVecOut = DatesPlot(VARopt.firstdate,nsteps,8*4,VARopt.frequency);
xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
%tstr = (string(t, 'yyyy-QQ'));
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)
xtickangle(90)


xlim([4 116])
savefig(strcat(myFigPath2, '_HD2_Area_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD2_Area_', num2str(ii)),'-pdf','-painters');                 

xlim([100 114])
savefig(strcat(myFigPath2, '_HD2_Area_Zoom_', num2str(ii)))
export_fig(strcat(myFigPath, '_HD2_Area_Zoom_', num2str(ii)),'-pdf','-painters');


temp1 = squeeze(HD2.shock(:,:,ii)); 
temp1Sum = sum(temp1, 2);
temp1SumAbs = sum(abs(temp1), 2);
temp1perc = abs(temp1)./temp1SumAbs;
temp1percWithSign = temp1perc.* sign(temp1);
tempSignalPercOfTotalPower = temp1Sum./sum(abs(temp1), 2);
tempXlsHDData = [(1995.0:0.25:1995+0.25*113)' temp1 temp1Sum temp1SumAbs temp1perc temp1percWithSign tempSignalPercOfTotalPower];
tempXlsHDDataTbl = array2table(tempXlsHDData, 'VariableNames', {'Date', 'Shock 1', 'Shock 2', 'Sum(Shocks) [i.e.Signal]' , 'Sum(abs(shocks))', 'Percentage Shock 1' , 'Percantage Shock 2' 'Signed Percentage Shock 1' , 'Signed Percantage Shock 2' 'Net Percentage Effect' });
writetable(tempXlsHDDataTbl, strcat(myDataPath, '_HD2_', '.xls'), 'Sheet', vnames{ii} );
writetable(array2table(hd2snames), strcat(myDataPath, '_HD2_', '.xls') , 'Sheet', 'Shocks' );

calcInd1 = 5; calcInd2 = 101;
tempWriteTableRatios = [
                [sum(abs( temp1(calcInd1:end, :) ))/ sum(sum(abs( temp1(calcInd1:end, :) )))       NaN];
                mean(abs( temp1(calcInd1:end, :)) ./    sum(abs( temp1(calcInd1:end, :)), 2) )      mean( abs(sum(temp1(calcInd1:end, :), 2)) ./ sum(abs(temp1(calcInd1:end, :)), 2)) ;

                [sum(abs( temp1(calcInd2:end, :) ))/ sum(sum(abs( temp1(calcInd2:end, :) )))       NaN];
                mean(abs( temp1(calcInd2:end, :)) ./     sum(abs( temp1(calcInd2:end, :)), 2) )      mean( abs(sum(temp1(calcInd2:end, :), 2)) ./ sum(abs(temp1(calcInd2:end, :)), 2)) ;
]; 
table2latex(( array2table(tempWriteTableRatios ,'RowNames', {'Ratio Method 1 1995-' , 'Ratio Method 2 1995-' , 'Ratio Method 1 2000-' , 'Ratio Method 2 2000-'  }, 'VariableNames', [vnames 'Efficiency']) ) , strcat(myTexPath , 'HDRatios', '_HD2_Var_', num2str(ii), '.tex'));



end


save(strcat(myDataPath, 'MyWorkspace_AfterHD2') )


%%
% ==============================  
% Fry og Pagan
%  ==============================  

% Set þetta í nýja möppu en ætti að vera hægt að sleppa þessu 

versionFolderName = '20240805_Greining4_FP/'
myFigPath2  = strcat('../Exported/', versionFolderName, 'ExportedFigures/')
myFigPath  = strcat('../', myFigPath2)
myTexPath  = strcat('../Exported/', versionFolderName , 'ExportedTex/' )
myDataPath = strcat('../Exported/', versionFolderName , 'Data/' )


% Hreinsa og sæki gögn 


xlsDataTbl = readtable('../Import/KeyInterestRate.xlsx','Sheet','ImportData');

xlsDataTemp = table2array(xlsDataTbl);
xlsDataTemp2 = [xlsDataTemp(1:end-1, :) diff(log(xlsDataTemp(:,2)))];

xlsData = [xlsDataTemp2(:,1) xlsDataTemp2(:,3)];
percMargfaldari = 1; %100;

ind= find(xlsData(:,1) == 1995.0); 

XfpPercentageWithTime = [xlsData(ind:ind+length(X)-1, 1) X*percMargfaldari xlsData(ind:ind+length(X)-1, 2)];
XfpWithTime = [XfpPercentageWithTime(:,1) XfpPercentageWithTime(:,2:end)/percMargfaldari];
Xfp = XfpWithTime(:,2:end);

vnamesfp = { 'VLF', 'VVLF', 'Vaxtastig'}; %{ 'GDP', 'GDP Deflator', 'Interest Rate',};
vnames_long_fp = vnamesfp;

figure
plot(XfpPercentageWithTime(:,1), XfpPercentageWithTime(:,2:end) * 100) %% Margfalda með 100 til að fá prósent eftir diff breytingu
xticks(XfpPercentageWithTime(1:4:114,1))
xlim([1995 2023.5])
xtickangle(90)
set(gca,'Layer','top');
legend(vnamesfp,'Orientation','Horizontal', 'Location', 'southoutside', 'Box','off', 'Fontsize', get(gca,'Fontsize')-1,'Interpreter', 'Tex')
set(gcf, 'Color', 'w');
ytickformat('percentage')
export_fig(strcat(myFigPath, 'FP_Grunngögn'),'-pdf','-painters');

%% Fry Pagan ADF og KPSS test

iVarNam = vnamesfp(3)
iVarNam = convertStringsToChars(iVarNam);
XfpIntTbl = array2table(Xfp(:,3), 'VariableNames', iVarNam);

FP_adftestWithDiff = adftest(XfpIntTbl, 'lags',0:8);
FP_adftestWithDiff = addvars(FP_adftestWithDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','ADF Variable Name');
table2latex(FP_adftestWithDiff ,  convertStringsToChars(convertCharsToStrings(strcat( myTexPath, 'FP_ADFTestWithDifference_', iVarNam, '.tex' ))))

FP_kpsstestWithDiff  = kpsstest(XfpIntTbl, 'lags',0:8,'trend',true );
FP_kpsstestWithDiff = addvars(FP_kpsstestWithDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','KPSS Test Variable Name');
table2latex(FP_kpsstestWithDiff ,  convertStringsToChars(convertCharsToStrings(strcat( myTexPath, 'FP_KPSSTestWithDifference_', iVarNam, '.tex' ))))


% Próf á móti upprunalegu gögnunum, þ.e. án umbreytingar
FP_adftestNoDiff = adftest(array2table(xlsDataTemp(:,2), 'VariableNames', iVarNam), 'lags',0:8);
FP_adftestNoDiff = addvars(FP_adftestNoDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','ADF Variable Name');
table2latex(FP_adftestNoDiff ,  convertStringsToChars(convertCharsToStrings(strcat( myTexPath, 'FP_ADFTestNoDifference_', iVarNam, '.tex' ))))

FP_kpsstestNoDiff  = kpsstest(array2table(xlsDataTemp(:,2), 'VariableNames', iVarNam), 'lags',0:8,'trend',true );
FP_kpsstestNoDiff = addvars(FP_kpsstestNoDiff,  repmat(iVarNam,9,1) , 'Before','h', 'NewVariableNames','KPSS Test Variable Name');
table2latex(FP_kpsstestNoDiff ,  convertStringsToChars(convertCharsToStrings(strcat( myTexPath, 'FP_KPSSTestNoDifference_', iVarNam, '.tex' ))))




%% Fry Pagan Lag order selection
det = 2;

[lagsVec, AIC, BIC, logL, HQIC] = VARlag_LeifiCustom(Xfp, 12,2);
LagOrderArray = [lagsVec, AIC, BIC, HQIC];
LagOrderTbl =  array2table(LagOrderArray );
LagOrderTbl.Properties.VariableNames(1:end) = {'Lag Order', 'AIC','BIC','HQIC'};
table2latex(LagOrderTbl, strcat(myTexPath, 'FP_LagOrder.tex'))

AIC_LagOrder = find(AIC==min(AIC)); % 7         %With diff 4 
BIC_LagOrder = find(BIC==min(BIC)); % 2         %With diff 1
HQIC_LagOrder = find(HQIC==min(HQIC)); %4       %With diff 1

[VAR_AIC, VARopt_AIC] = VARmodel(Xfp,AIC_LagOrder,det);
[VAR_BIC, VARopt_BIC] = VARmodel(Xfp,BIC_LagOrder,det);
[VAR_HQIC, VARopt_HQIC] = VARmodel(Xfp,HQIC_LagOrder,det);

%%
lagOrderRow = [AIC_LagOrder, BIC_LagOrder, HQIC_LagOrder];

% Engle
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = archtest(VAR_AIC.resid(: ,1), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = archtest(VAR_BIC.resid(: ,1), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = archtest(VAR_HQIC.resid(:,1), Lags=1:8); 
engleTest1_h = [h_AIC' h_BIC' h_HQIC'];
engleTest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = archtest(VAR_AIC.resid(: ,2), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = archtest(VAR_BIC.resid(: ,2), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = archtest(VAR_HQIC.resid(:,2), Lags=1:8); 
engleTest2_h = [h_AIC' h_BIC' h_HQIC'];
engleTest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = archtest(VAR_AIC.resid(: ,3), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = archtest(VAR_BIC.resid(: ,3), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = archtest(VAR_HQIC.resid(:,3), Lags=1:8); 
engleTest3_h = [h_AIC' h_BIC' h_HQIC'];
engleTest3_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];


% Ljung-Box Q-Test
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = lbqtest(VAR_AIC.resid(: ,1), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = lbqtest(VAR_BIC.resid(: ,1), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = lbqtest(VAR_HQIC.resid(:,1), Lags=1:8); 
lbqtest1_h = [h_AIC' h_BIC' h_HQIC'];
lbqtest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = lbqtest(VAR_AIC.resid(: ,2), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = lbqtest(VAR_BIC.resid(: ,2), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = lbqtest(VAR_HQIC.resid(:,2), Lags=1:8); 
lbqtest2_h = [h_AIC' h_BIC' h_HQIC'];
lbqtest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = lbqtest(VAR_AIC.resid(: ,3), Lags=1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = lbqtest(VAR_BIC.resid(: ,3), Lags=1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = lbqtest(VAR_HQIC.resid(:,3), Lags=1:8); 
lbqtest3_h = [h_AIC' h_BIC' h_HQIC'];
lbqtest3_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

% One-sample Kolmogorov-Smirnov test
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = kstest(VAR_AIC.resid(: ,1)); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = kstest(VAR_BIC.resid(: ,1)); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = kstest(VAR_HQIC.resid(:,1)); 
kstest1_h = [h_AIC' h_BIC' h_HQIC'];
kstest1_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = kstest(VAR_AIC.resid(: ,2)); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = kstest(VAR_BIC.resid(: ,2)); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = kstest(VAR_HQIC.resid(:,2)); 
kstest2_h = [h_AIC' h_BIC' h_HQIC'];
kstest2_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];

[h_AIC, pValue_AIC,stat_AIC,cValue_AIC]     = kstest(VAR_AIC.resid(: ,3)); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC]     = kstest(VAR_BIC.resid(: ,3)); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC] = kstest(VAR_HQIC.resid(:,3)); 
kstest3_h = [h_AIC' h_BIC' h_HQIC'];
kstest3_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];


% Multivariate ljung-box Q-test
%[h,pValue,stat,cValue] = LeifiCustom_mlbqtest(resid1,[3,6,10]); 
[h_AIC, pValue_AIC,stat_AIC,cValue_AIC, dof_AIC]     = LeifiCustom_mlbqtest(VAR_AIC.resid    , 1:8); 
[h_BIC, pValue_BIC,stat_BIC,cValue_BIC, dof_BIC]     = LeifiCustom_mlbqtest(VAR_BIC.resid    , 1:8); 
[h_HQIC, pValue_HQIC,stat_HQIC,cValue_HQIC, dof_HQIC] = LeifiCustom_mlbqtest(VAR_HQIC.resid   , 1:8); 
mlbqtest_h = [h_AIC' h_BIC' h_HQIC'];
mlbqtest_pValue = [pValue_AIC' pValue_BIC' pValue_HQIC'];
mlbqtest_dof = [dof_AIC' dof_BIC' dof_HQIC'];
mlbqtestDesc =  'Multivariate Portmanteau test with ';

%%

finalLagOrderArray = [lagOrderRow;
    engleTest1_h;  engleTest1_pValue ; engleTest2_h; engleTest2_pValue ; engleTest3_h; engleTest3_pValue ;
    lbqtest1_h; lbqtest1_pValue;  lbqtest2_h; lbqtest2_pValue; lbqtest3_h; lbqtest3_pValue; 
    kstest1_h;kstest1_pValue;kstest2_h;kstest2_pValue; kstest3_h;kstest3_pValue;
    mlbqtest_h;mlbqtest_pValue; mlbqtest_dof;
    %%omnibustest_h;omnibustest_pValue;
];

LagOrderSelectionTbl = array2table(finalLagOrderArray);
LagOrderSelectionTbl.Properties.VariableNames(1:end) = {'AIC','BIC','HQIC'};
LagOrderSelectionTbl.Properties.RowNames = {'Lag Order' ...
    'Engle Test 1 h Lags 1' ...
    'Engle Test 1 h Lags 2' ...
    'Engle Test 1 h Lags 3' ...
    'Engle Test 1 h Lags 4' ...
    'Engle Test 1 h Lags 5' ...
    'Engle Test 1 h Lags 6' ...
    'Engle Test 1 h Lags 7' ...
    'Engle Test 1 h Lags 8' ...
    'Engle Test 1 p-value Lags 1' ...
    'Engle Test 1 p-value Lags 2' ...
    'Engle Test 1 p-value Lags 3' ...
    'Engle Test 1 p-value Lags 4' ...
    'Engle Test 1 p-value Lags 5' ...
    'Engle Test 1 p-value Lags 6' ...
    'Engle Test 1 p-value Lags 7' ...
    'Engle Test 1 p-value Lags 8' ...
    'Engle Test 2 h Lags 1' ...
    'Engle Test 2 h Lags 2' ...
    'Engle Test 2 h Lags 3' ...
    'Engle Test 2 h Lags 4' ...
    'Engle Test 2 h Lags 5' ...
    'Engle Test 2 h Lags 6' ...
    'Engle Test 2 h Lags 7' ...
    'Engle Test 2 h Lags 8' ...
    'Engle Test 2 p-value Lags 1' ...
    'Engle Test 2 p-value Lags 2' ...
    'Engle Test 2 p-value Lags 3' ...
    'Engle Test 2 p-value Lags 4' ...
    'Engle Test 2 p-value Lags 5' ...
    'Engle Test 2 p-value Lags 6' ...
    'Engle Test 2 p-value Lags 7' ...
    'Engle Test 2 p-value Lags 8' ...
    'Engle Test 3 h Lags 1' ...
    'Engle Test 3 h Lags 2' ...
    'Engle Test 3 h Lags 3' ...
    'Engle Test 3 h Lags 4' ...
    'Engle Test 3 h Lags 5' ...
    'Engle Test 3 h Lags 6' ...
    'Engle Test 3 h Lags 7' ...
    'Engle Test 3 h Lags 8' ...
    'Engle Test 3 p-value Lags 1' ...
    'Engle Test 3 p-value Lags 2' ...
    'Engle Test 3 p-value Lags 3' ...
    'Engle Test 3 p-value Lags 4' ...
    'Engle Test 3 p-value Lags 5' ...
    'Engle Test 3 p-value Lags 6' ...
    'Engle Test 3 p-value Lags 7' ...
    'Engle Test 3 p-value Lags 8' ...
    ...
    'Ljung-Box Q-Test 1 h Lags 1' ...
    'Ljung-Box Q-Test 1 h Lags 2' ...
    'Ljung-Box Q-Test 1 h Lags 3' ...
    'Ljung-Box Q-Test 1 h Lags 4' ...
    'Ljung-Box Q-Test 1 h Lags 5' ...
    'Ljung-Box Q-Test 1 h Lags 6' ...
    'Ljung-Box Q-Test 1 h Lags 7' ...
    'Ljung-Box Q-Test 1 h Lags 8' ...
    'Ljung-Box Q-Test 1 p-value Lags 1' ...
    'Ljung-Box Q-Test 1 p-value Lags 2' ...
    'Ljung-Box Q-Test 1 p-value Lags 3' ...
    'Ljung-Box Q-Test 1 p-value Lags 4' ...
    'Ljung-Box Q-Test 1 p-value Lags 5' ...
    'Ljung-Box Q-Test 1 p-value Lags 6' ...
    'Ljung-Box Q-Test 1 p-value Lags 7' ...
    'Ljung-Box Q-Test 1 p-value Lags 8' ...
    'Ljung-Box Q-Test 2 h Lags 1' ...
    'Ljung-Box Q-Test 2 h Lags 2' ...
    'Ljung-Box Q-Test 2 h Lags 3' ...
    'Ljung-Box Q-Test 2 h Lags 4' ...
    'Ljung-Box Q-Test 2 h Lags 5' ...
    'Ljung-Box Q-Test 2 h Lags 6' ...
    'Ljung-Box Q-Test 2 h Lags 7' ...
    'Ljung-Box Q-Test 2 h Lags 8' ...
    'Ljung-Box Q-Test 2 p-value Lags 1' ...
    'Ljung-Box Q-Test 2 p-value Lags 2' ...
    'Ljung-Box Q-Test 2 p-value Lags 3' ...
    'Ljung-Box Q-Test 2 p-value Lags 4' ...
    'Ljung-Box Q-Test 2 p-value Lags 5' ...
    'Ljung-Box Q-Test 2 p-value Lags 6' ...
    'Ljung-Box Q-Test 2 p-value Lags 7' ...
    'Ljung-Box Q-Test 2 p-value Lags 8' ...
    'Ljung-Box Q-Test 3 h Lags 1' ...
    'Ljung-Box Q-Test 3 h Lags 2' ...
    'Ljung-Box Q-Test 3 h Lags 3' ...
    'Ljung-Box Q-Test 3 h Lags 4' ...
    'Ljung-Box Q-Test 3 h Lags 5' ...
    'Ljung-Box Q-Test 3 h Lags 6' ...
    'Ljung-Box Q-Test 3 h Lags 7' ...
    'Ljung-Box Q-Test 3 h Lags 8' ...
    'Ljung-Box Q-Test 3 p-value Lags 1' ...
    'Ljung-Box Q-Test 3 p-value Lags 2' ...
    'Ljung-Box Q-Test 3 p-value Lags 3' ...
    'Ljung-Box Q-Test 3 p-value Lags 4' ...
    'Ljung-Box Q-Test 3 p-value Lags 5' ...
    'Ljung-Box Q-Test 3 p-value Lags 6' ...
    'Ljung-Box Q-Test 3 p-value Lags 7' ...
    'Ljung-Box Q-Test 3 p-value Lags 8' ...
    ...
'Kolmogorov-Smirnov Test 1 h'           ...
'Kolmogorov-Smirnov Test 1 P-Value'     ...
'Kolmogorov-Smirnov Test 2 h'           ...
'Kolmogorov-Smirnov Test 2 P-Value'     ...
'Kolmogorov-Smirnov Test 3 h'           ...
'Kolmogorov-Smirnov Test 3 P-Value'     ...
...
...
    'Multivariate Portmanteau test h Lags 1' ...
    'Multivariate Portmanteau test h Lags 2' ...
    'Multivariate Portmanteau test h Lags 3' ...
    'Multivariate Portmanteau test h Lags 4' ...
    'Multivariate Portmanteau test h Lags 5' ...
    'Multivariate Portmanteau test h Lags 6' ...
    'Multivariate Portmanteau test h Lags 7' ...
    'Multivariate Portmanteau test h Lags 8' ...
    'Multivariate Portmanteau test p-value Lags 1' ...
    'Multivariate Portmanteau test p-value Lags 2' ...
    'Multivariate Portmanteau test p-value Lags 3' ...
    'Multivariate Portmanteau test p-value Lags 4' ...
    'Multivariate Portmanteau test p-value Lags 5' ...
    'Multivariate Portmanteau test p-value Lags 6' ...
    'Multivariate Portmanteau test p-value Lags 7' ...
    'Multivariate Portmanteau test p-value Lags 8' ...
    'Multivariate Portmanteau test DF Lags 1' ...
    'Multivariate Portmanteau test DF Lags 2' ...
    'Multivariate Portmanteau test DF Lags 3' ...
    'Multivariate Portmanteau test DF Lags 4' ...
    'Multivariate Portmanteau test DF Lags 5' ...
    'Multivariate Portmanteau test DF Lags 6' ...
    'Multivariate Portmanteau test DF Lags 7' ...
    'Multivariate Portmanteau test DF Lags 8' ...
    }


table2latex(LagOrderSelectionTbl, strcat(myTexPath, 'FP_LagOrderSelection.tex'));


%% Fry Pagan - Plot test for residuals

for iac = 1:3

if (iac==1)  
    tacVar = VAR_AIC;
    tacStr = 'AIC';
elseif(iac==2)
    tacVar = VAR_BIC;
     tacStr = 'BIC';
elseif(iac==3)
    tacVar = VAR_HQIC;
    tacStr = 'HQIC';
end

temvnames = vnames_long_fp;

figure
set(gcf, 'Color', 'w');set(gca,'Layer','top');
subplot(3,1,1)
tempResid = tacVar.resid(:,1);
[acf,lags,bounds] = autocorr(tempResid);
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(temvnames(1))

subplot(3,1,2)
tempResid = tacVar.resid(:,2);
[acf,lags,bounds] = autocorr(tempResid);
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(temvnames(2))

subplot(3,1,3)
tempResid = tacVar.resid(:,3);
[acf,lags,bounds] = autocorr(tempResid);
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(temvnames(3))
export_fig(strcat(myFigPath, 'FP_Autocorrelation_' , tacStr),'-pdf','-painters');

figure
set(gcf, 'Color', 'w');set(gca,'Layer','top');
subplot(3,1,1)
[acf,lags,bounds] = crosscorr(tacVar.resid(:,1), tacVar.resid(:,2));
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(strcat('Cross correlation (', temvnames(1), ' vs ', temvnames(2), ')'))

subplot(3,1,2)
[acf,lags,bounds] = crosscorr(tacVar.resid(:,1), tacVar.resid(:,3));
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(strcat('Cross correlation (', temvnames(1), ' vs ', temvnames(3), ')'))

subplot(3,1,3)
[acf,lags,bounds] = crosscorr(tacVar.resid(:,2), tacVar.resid(:,3));
stem(lags,acf); xlabel('Lag'); %ylabel('\rho(k)');
hold on;
set(line(lags,bounds(1)*ones(length(acf),1))    ,'color',[1 0 0]);
set(line(lags,bounds(2)*ones(length(acf),1))    ,'color',[1 0 0]);
title(strcat('Cross correlation (', temvnames(2), ' vs ', temvnames(3), ')'))
export_fig(strcat(myFigPath, 'FP_CrossCorrelation_' , tacStr),'-pdf','-painters');


figure
set(gcf, 'Color', 'w');
qqplot(tacVar.resid)
title('')
export_fig(strcat(myFigPath, 'FP_QQPlot' , tacStr),'-pdf','-painters');

figure
set(gcf, 'Color', 'w');
subplot(3,1,1); histogram(tacVar.resid(:,1))
subplot(3,1,2); histogram(tacVar.resid(:,2))
subplot(3,1,3); histogram(tacVar.resid(:,3))
export_fig(strcat(myFigPath, 'FP_Histogram' , tacStr),'-pdf','-painters');

end


%%

[nobs, nvar] = size(Xfp);

det = 1; % no constant / constant / constant and trend
nlags = 4;  
[VAR, VARopt] = VARmodel(Xfp,nlags,det);
VARopt.vnames = vnamesfp;
VARopt.firstdate = 1995.00+nlags*0.25;       % initial date of the sample in format 1999.75 => 1999Q4 (both for annual and quarterly data)
VARopt.frequency = 'q'; 
[beta, TABLE] = VARprint(VAR,VARopt);

varRFETable = array2table( TABLE(2:end, 2:end));
    varRFETable.Properties.VariableNames =  TABLE(1, 2:end);
    varRFETable.Properties.RowNames = TABLE(2:end, 1);
varSigmaTable = array2table(VAR.sigma);
varBetaTable = array2table( beta(2:end, 2:end));
    varBetaTable.Properties.VariableNames =  beta(1, 2:end);
    varBetaTable.Properties.RowNames = beta(2:end, 1);
table2latex(varRFETable, strcat(myTexPath, 'FP_varRFETable.tex'));
    table2latex(varSigmaTable, strcat(myTexPath, 'FP_varSigmaTable.tex'));
    table2latex(varBetaTable, strcat(myTexPath, 'FP_varBetaTable.tex'))


VARopt.snames = {'Eftirspurn', 'Framboð', 'Peningastefna'}; %{'Demand', 'Cost-Push', 'Monetary Shock'}; % Shock names
VARopt.ident   = 'sign'; 
VARopt.sr_mod  = 1; % 0 no Bayesian draw / 1 Bayesian


VARopt.sr_hor = 1;  % 
VARopt.nsteps = 40;
VARopt.ndraws  = 5000*10; 
  

R = [ 1   -1    -1   
      1   1     -1
      1   1      1   ];  

VARopt.sr_draw = 100000*1000;
VARopt.sr_rot = 500*2; %Hámarksfjöldi snúninga

disp('Begin SRout_FP=SR(...)')
temCurrTime = datetime;
disp(temCurrTime)
[SRout_FP, VARDraws_FP] = SR(VAR,R,VARopt);
disp(datetime)
disp(datetime-temCurrTime)
disp('End SRout_FP=SR(...)')


%%
VARopt.figname = [myFigPath 'FryAndPagan_']; 
% Plot Impulse response
VARirplot(SRout_FP.IRmed, VARopt, SRout_FP.IRinf, SRout_FP.IRsup )


%% Plot IR

quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

IR = SRout_FP.IRmed;
INF = SRout_FP.IRinf;
SUP = SRout_FP.IRsup;
[nsteps, nvars, nshocks] = size(IR);

steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);

SwatheOpt = PlotSwatheOption;
SwatheOpt.marker = '*';
SwatheOpt.trans = 1;

figure; FigSize(26,24); clf;
for jj=1:nshocks                
    for ii=1:nvars
        subplot(3,1,ii)

        plot(steps,IR(:,ii,jj),'LineStyle','-','Color','k','LineWidth',2,'Marker',SwatheOpt.marker); hold on
        if exist('INF','var') && exist('SUP','var')
            PlotSwathe(IR(:,ii,jj),[INF(:,ii,jj) SUP(:,ii,jj)],SwatheOpt); hold on;
        end
        plot(x_axis,'--k','LineWidth',0.5); hold on
        xlim([1 nsteps]);
        title([VARopt.vnames{ii} ' (Rykkur: ' VARopt.snames{jj} ' ) '], 'FontWeight','bold','FontSize',10); 
        set(gca, 'Layer', 'top');
    end
        set(gcf, 'Color', 'w');
        export_fig(strcat(myFigPath, 'FP_IR_Shock', num2str(jj),'_Var_All'),'-pdf','-painters');
    clf('reset');
end

close all


%% Fry and Pagan HD plot

HDFP = SRout_FP.HD;

[nsteps, nvars, nshocks] = size(HDFP.shock);
for ii=1:nvars

    figure; FigSize(27,15); clf;
    H = BarPlot(squeeze(HDFP.shock(:,:,ii))); hold on; 
    h = plot( sum(squeeze(HDFP.shock(:,:,ii)),2),'-k','LineWidth',2);
    xlim([1 nsteps]);
	set(gca,'Layer','top');
    title([vnamesfp{ii}], 'FontWeight','bold','FontSize',10); 
    % Save
       if 0==1 %suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            tempSubtitle = convertCharsToStrings( [Alphabet(ii) ') HD of '  vnamesfp{ii}]);
            SupTitle(tempSubtitle)
        end
        opt = LegOption; opt.handle = [H(1,:) h];
        LegSubplot([VARopt.snames {'Gögn'}],opt);
        set(gcf, 'Color', 'w');

xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)
xtickangle(90)


xtickangle(90)
xlim([4 116])
savefig(strcat(myFigPath2, '_FP_HD_', num2str(ii)))
export_fig(strcat(myFigPath, '_FP_HD_', num2str(ii)),'-pdf','-painters');                 

xlim([100 114])
savefig(strcat(myFigPath2, '_FP_HD_Zoom_', num2str(ii)))
export_fig(strcat(myFigPath, '_FP_HD_Zoom_', num2str(ii)),'-pdf','-painters');

%figure
figure; FigSize(27,15); clf;  
H = AreaPlot(squeeze(HDFP.shock(:,:,ii))); hold on; 
h = plot( sum(squeeze(HDFP.shock(:,:,ii)),2),'-k','LineWidth',2);
set(gca,'Layer','top');
title([vnamesfp{ii}], 'FontWeight','bold','FontSize',10); 
opt = LegOption; opt.handle = [H(1,:) h];
LegSubplot([VARopt.snames {'Gögn'}],opt);
set(gcf, 'Color', 'w');

xticks(1:1:114)
t = [datetime(1995,1,31):calmonths(3):datetime(2023,05,31)];
tstr = (string(t, 'yyyy'));
tstr(2:4:end) = '';
tstr(3:4:end) = '';
tstr(4:4:end) = '';
xticklabels(tstr)
xtickangle(90)


xtickangle(90)
xlim([4 116])
savefig(strcat(myFigPath2, '_FP_HD_Area_', num2str(ii)))
export_fig(strcat(myFigPath, '_FP_HD_Area_', num2str(ii)),'-pdf','-painters');                 
% 
xlim([100 114])
savefig(strcat(myFigPath2, '_FP_HD_Area_Zoom_', num2str(ii)))
export_fig(strcat(myFigPath, '_FP_HD_Area_Zoom_', num2str(ii)),'-pdf','-painters');


temp1 = squeeze(HDFP.shock(:,:,ii)); 
temp1Sum = sum(temp1, 2);
temp1SumAbs = sum(abs(temp1), 2);
temp1perc = abs(temp1)./temp1SumAbs;
temp1percWithSign = temp1perc.* sign(temp1);
tempSignalPercOfTotalPower = temp1Sum./sum(abs(temp1), 2);
tempXlsHDData = [(1995.0:0.25:1995+0.25*113)' temp1 temp1Sum temp1SumAbs temp1perc temp1percWithSign tempSignalPercOfTotalPower];
tempXlsHDDataTbl = array2table(tempXlsHDData, 'VariableNames', {'Date', 'Shock 1', 'Shock 2', 'Shock 3', 'Sum(Shocks) [i.e.Signal]' , 'Sum(abs(shocks))', 'Percentage Shock 1' , 'Percantage Shock 2', 'Percantage Shock 3', 'Signed Percentage Shock 1' , 'Signed Percantage Shock 2', 'Signed Percantage Shock 3' ,'Net Percentage Effect' });
writetable(tempXlsHDDataTbl, strcat(myDataPath, '_FP_HD_', '.xls'), 'Sheet', vnamesfp{ii} );
writetable(array2table(VARopt.snames), strcat(myDataPath, '_FP_HD_', '.xls') , 'Sheet', 'Shocks' );

calcInd1 = 5; calcInd2 = 101;
tempWriteTableRatios = [
                [sum(abs( temp1(calcInd1:end, :) ))/ sum(sum(abs( temp1(calcInd1:end, :) )))       NaN];
                mean(abs( temp1(calcInd1:end, :)) ./    sum(abs( temp1(calcInd1:end, :)), 2) )      mean( abs(sum(temp1(calcInd1:end, :), 2)) ./ sum(abs(temp1(calcInd1:end, :)), 2)) ;

                [sum(abs( temp1(calcInd2:end, :) ))/ sum(sum(abs( temp1(calcInd2:end, :) )))       NaN];
                mean(abs( temp1(calcInd2:end, :)) ./     sum(abs( temp1(calcInd2:end, :)), 2) )      mean( abs(sum(temp1(calcInd2:end, :), 2)) ./ sum(abs(temp1(calcInd2:end, :)), 2)) ;
]; 
table2latex(( array2table(tempWriteTableRatios ,'RowNames', {'Ratio Method 1 1995-' , 'Ratio Method 2 1995-' , 'Ratio Method 1 2000-' , 'Ratio Method 2 2000-'  }, 'VariableNames', [vnamesfp 'Efficiency']) ) , strcat(myTexPath , 'HDRatios', '_FP_HD_Var_', num2str(ii), '.tex'));




end

save(strcat(myDataPath, 'MyWorkspace_AfterFP') )


% FIN