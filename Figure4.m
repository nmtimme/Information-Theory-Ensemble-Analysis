%% FIGURE4
%   Creates a plot of p-values from comparisons between pairs of ensembles
%   with different signal values across multiple parameters. This is Figure
%   4 in the manuscript.
%
% Other m-files required: ensModel, ensCompModel
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% May 2017; Last revision: 5-Jul-2018


%% Set Parameters

% Set the seed to produce the exact figure shown in the paper. Remove this
% line to allow for stochastic fluctuations.
rng(1)

% Set the input parameters to use in the figure (nEns, nTrials, sigStr,
% noisStr, nMC) for each comparison.
setParams = cell([9,1]);
setParams{1} = [10,20,0.2,0.5,100;...
    10,20,0.5,0.5,100;...
    10,20,1,0.5,100];
setParams{2} = [40,20,0.2,0.5,100;...
    40,20,0.5,0.5,100;...
    40,20,1,0.5,100];
setParams{3} = [160,20,0.2,0.5,100;...
    160,20,0.5,0.5,100;...
    160,20,1,0.5,100];
setParams{4} = [10,60,0.2,0.5,100;...
    10,60,0.5,0.5,100;...
    10,60,1,0.5,100];
setParams{5} = [40,60,0.2,0.5,100;...
    40,60,0.5,0.5,100;...
    40,60,1,0.5,100];
setParams{6} = [160,60,0.2,0.5,100;...
    160,60,0.5,0.5,100;...
    160,60,1,0.5,100];
setParams{7} = [10,100,0.2,0.5,100;...
    10,100,0.5,0.5,100;...
    10,100,1,0.5,100];
setParams{8} = [40,100,0.2,0.5,100;...
    40,100,0.5,0.5,100;...
    40,100,1,0.5,100];
setParams{9} = [160,100,0.2,0.5,100;...
    160,100,0.5,0.5,100;...
    160,100,1,0.5,100];

% Set the signal strength (proportion of nTrials/4 trials that will be
% flipped to create mutual information between X and Y) for reference
% models that will be compared to the models listed in setParams.
sigStrRef = 0:0.1:1;

% Set the number of Monte Carlo trials for the comparisons between the
% weighted means of the data sets
nMCcomp = 1000;

% Set the number of models to run for each comparison.
nModels = 50;
    



%% Run the Model

% Preallocate space
compp = NaN([9,3,length(sigStrRef),nModels]);

% Loop through all the sets
if exist([pwd,'\Figure4Data.mat'],'file') ~= 2
    
    % Warn the user that this code takes a while to run
    disp('Please note, this script can take ~1 hour to run.')
    
    for iComp = 1:length(setParams)
        for iSet = 1:3
            for isigStr = 1:length(sigStrRef)
                
                % Assign parameters
                nEns = [setParams{iComp}(iSet,1),setParams{iComp}(iSet,1)];
                nTrials = [setParams{iComp}(iSet,2),setParams{iComp}(iSet,2)];
                sigStr = [setParams{iComp}(iSet,3),sigStrRef(isigStr)];
                noisStr = [setParams{iComp}(iSet,4),setParams{iComp}(iSet,4)];
                nMC = [setParams{iComp}(iSet,5),setParams{iComp}(iSet,5)];
                
                % Run the models
                compp(iComp,iSet,isigStr,:) = ensCompModel(nEns, nTrials, sigStr, noisStr, nMC, nMCcomp, nModels);
                
            end
        end
    end

    save('Figure4Data')
    
else
    
    disp('A previously generated data file exists in this folder,')
    disp('so we will load that rather than rerun the model to save time.')
    disp('To rerun the model, delete the file Figure4Data.mat')
    
    load([pwd,'\Figure4Data.mat'])
    
end


%% Make the figure


% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.5;
rmargin = 0.1;

% Set the top and bottom margins in inches
tmargin = 0.45;
bmargin = 0.35;

% Set the horizontal distance between plots in inches
hspace = 0.5;

% Set the vertical distance between plots in inches
vspace = 0.75;

% Set Paper size
papersize = [7,6];
% papersize = [3 5];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures
figdim = [3,3];

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
MedLW = 1;
ThickLW = 1.2;

% Set the dot sizes
DotSz1 = 30;
DotSz2 = 20;

% Set the transparency factor
TransFact = 0.4;

% Set colors
PColor = {[1,0,0],[0,210/255,50/255],[0,0,1]};

% Find the limits of the data
yMax = -inf;
for iComp = 1:length(setParams)
    for iSet = 1:3
        yCoords = zeros([2,length(sigStrRef)]);
        for isigStr = 1:length(sigStrRef)
            
            temp = squeeze(compp(iComp,iSet,isigStr,:));
            temp = sort(temp,'ascend');
            nEl = length(temp);
            yCoords(:,isigStr) = [temp(ceil(0.25*nEl));temp(floor(0.75*nEl))];
            
        end
        yCoords = [yCoords(1,:),fliplr(yCoords(2,:))];
        yMax = max([yMax,max(yCoords)]);
    end
end
yMin = 1/(2*nMCcomp);
ldy = log10(yMax) - log10(yMin);

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

for iComp = 1:length(setParams)
    
    iCol = rem(iComp,3);
    iCol(iCol == 0) = 3;
    iRow = ceil(iComp/3);
    
    subplot('Position',[LeftCoord(iCol),BottomCoord(iRow),width,height]);
    hold on
    
    for iSet = 1:3
        xCoords = zeros([2,length(sigStrRef)]);
        yCoords = zeros([2,length(sigStrRef)]);
        yMed = zeros(size(sigStrRef));
        for isigStr = 1:length(sigStrRef)
            
            temp = squeeze(compp(iComp,iSet,isigStr,:));
            temp = sort(temp,'ascend');
            nEl = length(temp);
            xCoords(:,isigStr) = sigStrRef(isigStr);
            yCoords(:,isigStr) = [temp(ceil(0.25*nEl));temp(floor(0.75*nEl))];
            yMed(isigStr) = median(temp);
            
        end
        xCoords = [xCoords(1,:),fliplr(xCoords(2,:))];
        yCoords = [yCoords(1,:),fliplr(yCoords(2,:))];
%         patch(xCoords,yCoords,PColor{iSet})
        fill(xCoords,yCoords,PColor{iSet},'EdgeColor','none','facealpha',.3)
        
        plot(sigStrRef,yMed,'Color',PColor{iSet},'LineWidth',MedLW)
        
    end
    
    x1 = 0.05;
    x2 = 0.15;
    x3 = 0.18;
    y1 = 10.^(linspace(log10(yMax) + 0.3*ldy,log10(yMax) + 0.1*ldy,3));
    
    for iSet = 1:3
        line([x1,x2],[y1(iSet),y1(iSet)],'Color',PColor{iSet},'LineWidth',MedLW)
        text(x3,y1(iSet),['s = ',num2str(setParams{iComp}(iSet,3))],'FontSize',LegFS)
    end
    
    
    xlim([0,1])
    ylim([yMin,10^(log10(yMax) + 0.4*ldy)])
    
    set(gca,'Layer','Top')
    set(gca,'YScale','log')
    
    xlabel('Reference Signal Strength (s)','FontSize',AxLabelFS)
    ylabel('p-value','FontSize',AxLabelFS)
    title(['Weighted Mean MI Comparison',char(10),'(N_{ens} = ',num2str(setParams{iComp}(1,1)),', N_{obs} = ',num2str(setParams{iComp}(1,2)),', a = ',num2str(setParams{iComp}(1,4)),')'],'FontSize',TitleFS)
    set(gca,'FontSize',UnitFS)
    
end



% Save the figure
print(F1,'-dpdf','-painters','Figure4')
