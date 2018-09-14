%% FIGURE2
%   Creates a figure panel showing weighted mutual information and error
%   across different signal strengths, noise levels, number of
%   information sources, and number of observations. This is Figure 2 in
%   the manuscript.
%
% Other m-files required: ensModel
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2017; Last revision: 26-Apr-2018

%% Set Parameters

% Set the seed to produce the exact figure shown in the paper. Remove this
% line or change the argument to allow for stochastic fluctuations.
rng(4)

% Set the number of members of the ensemble
nEns = [10,40,160];

% Set the number of trials used to generate the probability distribution
% (must be a multiple of 4)
nObs = [20,60,100];

% Set the signal strength
sigStr = 0:0.1:1;

% Set the noise strength
noisStr = [0,0.2,0.5,1];

% Set the number of Monte Carlo trials for significance testing for each
% member of the ensemble (calculating the p-value of an individual
% information source's MI value).
nMC = 500;

% Set the p-value threshold to consider an MI distribution significantly
% different from null.
KSpThresh = 0.01;



%% Run the Model

% Preallocate space
MI = cell([length(nEns),length(nObs),length(sigStr),length(noisStr)]);
MIp = cell([length(nEns),length(nObs),length(sigStr),length(noisStr)]);
nullMI = cell([length(nEns),length(nObs),length(sigStr),length(noisStr)]);
nullMIp = cell([length(nEns),length(nObs),length(sigStr),length(noisStr)]);
KSp = NaN([length(nEns),length(nObs),length(sigStr),length(noisStr)]);
MIPlot = NaN([3,length(nEns),length(nObs),length(sigStr),length(noisStr)]);

if exist([pwd,'\Figure2Data.mat'],'file') ~= 2
    for iEns = 1:length(nEns)
        for iTrials = 1:length(nObs)
            for isigStr = 1:length(sigStr)
                for inoisStr = 1:length(noisStr)
                    
                    % Assign parameters
                    nModelSets = 1;
                    
                    % Run the model
                    [MIPlot(1,iEns,iTrials,isigStr,inoisStr),...
                        MIPlot(2,iEns,iTrials,isigStr,inoisStr),...
                        MIPlot(3,iEns,iTrials,isigStr,inoisStr),...
                        MI{iEns,iTrials,isigStr,inoisStr},...
                        MIp{iEns,iTrials,isigStr,inoisStr},...
                        nullMI{iEns,iTrials,isigStr,inoisStr},...
                        nullMIp{iEns,iTrials,isigStr,inoisStr},...
                        KSp(iEns,iTrials,isigStr,inoisStr)] = ensModel(nEns(iEns), nObs(iTrials), sigStr(isigStr), noisStr(inoisStr), nMC, nModelSets);
                    
                    
                end
            end
        end
    end
    
    
    save('Figure2Data')
else
    
    disp('A previously generated data file exists in this folder,')
    disp('so we will load that rather than rerun the model to save time.')
    disp('To rerun the model, delete the file Figure2Data.mat')
    load([pwd,'\Figure2Data.mat'])
    
end





%% Make the Figure

% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.5;
rmargin = 0.1;

% Set the top and bottom margins in inches
tmargin = 0.3;
bmargin = 0.35;

% Set the horizontal distance between plots in inches
hspace = 0.5;

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7 6];
% papersize = [3 5];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures
figdim = [length(nObs),length(nEns)];

% calculate the width and height based on dimensions above in inches
width = (papersize(1) - lmargin - rmargin - (figdim(1) - 1)*hspace)/figdim(1);
height = (papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace)/figdim(2);

LeftCoord = (0:(figdim(1) - 1))*(width + hspace) + lmargin;
BottomCoord = (0:(figdim(2) - 1))*(height + vspace) + bmargin;
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

% Set fontsizes
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 7;
LegFS = 7;


% Set the linewidths
ThinLW = 0.6;
MedLW = 0.8;
ThickLW = 1.2;

% Set the dot sizes
DotSz1 = 15;

% Make some plot colors
PColor = {[200/255,0,255/255],[1,0,0],[0,210/255,50/255],[0,0,1]};

% Calculate the offset to make data points visible
delta = 0.01;
OffSet = (-((length(noisStr) - 1)*delta)/2):delta:(((length(noisStr) - 1)*delta)/2);

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

for inEns = 1:length(nEns)
    for inObs = 1:length(nObs)
        
        subplot('Position',[LeftCoord(inObs),BottomCoord(inEns),width,height]);
        hold on
        
        h = zeros([1,length(noisStr)]);
        M = cell([1,length(noisStr)]);
        for inoisStr = 1:length(noisStr)
            
            h(inoisStr) = plot(sigStr + OffSet(inoisStr),squeeze(MIPlot(1,inEns,inObs,:,inoisStr)),'Color',PColor{inoisStr});
            
            for isigStr = 1:length(sigStr)
                
                line([sigStr(isigStr) + OffSet(inoisStr),sigStr(isigStr) + OffSet(inoisStr)],...
                    [MIPlot(1,inEns,inObs,isigStr,inoisStr) - MIPlot(2,inEns,inObs,isigStr,inoisStr),MIPlot(1,inEns,inObs,isigStr,inoisStr) + MIPlot(2,inEns,inObs,isigStr,inoisStr)],...
                    'Color',PColor{inoisStr})
                if KSp(inEns,inObs,isigStr,inoisStr) < KSpThresh
                    scatter(sigStr(isigStr) + OffSet(inoisStr),MIPlot(1,inEns,inObs,isigStr,inoisStr),DotSz1,...
                        'MarkerEdgeColor',PColor{inoisStr},'MarkerFaceColor',PColor{inoisStr})
                else
                    scatter(sigStr(isigStr) + OffSet(inoisStr),MIPlot(1,inEns,inObs,isigStr,inoisStr),DotSz1,...
                        'MarkerEdgeColor',PColor{inoisStr},'MarkerFaceColor','none')
                end
            end
            
            M{inoisStr} = ['a = ',num2str(noisStr(inoisStr))];
        end
%         legend(h,M,'Location','NorthWest')
        
        if isequal([inEns,inObs],[1,1])
           x1 = 0.05;
           x2 = 0.1;
           y1 = linspace(1.05,0.5,6);
        
           scatter(x1,y1(1),DotSz1,'k','filled')
           text(x2,y1(1),'Sig. Dif. from Null','FontSize',LegFS)
           scatter(x1,y1(2),DotSz1,'k')
           text(x2,y1(2),'Not Sig. Dif. from Null','FontSize',LegFS)
           for inoisStr = 1:length(noisStr)
               scatter(x1,y1(2 + inoisStr),DotSz1,PColor{inoisStr},'filled')
                text(x2,y1(2 + inoisStr),M{inoisStr},'FontSize',LegFS)
           end
        end
        
        xlabel('Signal Strength (s)','FontSize',AxLabelFS)
        ylabel('Weighted MI (bits)','FontSize',AxLabelFS)
        title(['Model MI (N_{ens} = ',num2str(nEns(inEns)),', N_{obs} = ',num2str(nObs(inObs)),')'],'FontSize',TitleFS)
        set(gca,'Layer','Top')
        set(gca,'FontSize',UnitFS)
        
        xlim([-5*delta,1 + 5*delta])
        ylim([0,1.1])
        
        
    end
end

% Save the figure
print(F1,'-dpdf','-painters','Figure2')
