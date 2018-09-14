%% FIGURE1
%   Creates an explanatory plot that shows the distributions of information
%   values for three model ensembles, as well as comparisons between the
%   ensembles and significance testing for the ensembles. This is Figure 1
%   in the manuscript.
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

% Set the input parameters to use in the figure (nEns, nObs, sigStr,
% noisStr, nMC)
setParams = [40,60,0.4,0.2,100;...
    40,60,0.4,0.5,100;...
    40,60,0.4,1,100];

% Set the number of Monte Carlo trials for the comparisons between the
% weighted means of the data sets
nMCcomp = 1000;
    



%% Run the Model

% Preallocate variables
MI = cell([1,3]);
MIp = cell([1,3]);
nullMI = cell([1,3]);
nullMIp = cell([1,3]);
KSp = NaN([1,3]);
MIPlot = NaN([6,3]);

% Loop through all the sets
for iSet = 1:3
    
    % Assign parameters
    nEns = setParams(iSet,1);
    nObs = setParams(iSet,2);
    sigStr = setParams(iSet,3);
    noisStr = setParams(iSet,4);
    nMC = setParams(iSet,5);
    nModelSets = 1;
    
    % Run the model
    [MIPlot(1,iSet),MIPlot(2,iSet),MIPlot(3,iSet),MI{iSet},MIp{iSet},nullMI{iSet},nullMIp{iSet},KSp(iSet)] = ensModel(nEns, nObs, sigStr, noisStr, nMC, nModelSets);
    
    % Calculate the weighted mean and standard error of the
    % weighted mean for the null data
    MIVals = nullMI{iSet}(:);
    weights = -log10(nullMIp{iSet}(:));
    if isequal(weights,zeros(size(weights)))
        weights = ones(size(weights));
    end
    weights = weights./sum(weights);
    MIPlot(4,iSet) = sum(MIVals.*weights); % Weighted mean
    MIPlot(5,iSet) = std(MIVals)*sqrt(sum(weights.^2)); % Standard error of the weighted mean
    MIPlot(6,iSet) = sqrt(sum(weights.*(MIVals - sum(MIVals.*weights)).^2)); % Weighted standard deviation
    
end

% Run the comparisons between the weighted means of the sets
compEx = NaN([3,3]);
for iSet = 1:2
    for jSet = (iSet + 1):3
        
        % Get the data
        iMI = MI{iSet};
        iWeights = -log10(MIp{iSet});
        if isequal(iWeights,zeros(size(iWeights)))
            iWeights = ones(size(iWeights));
        end
        iWeights = iWeights./sum(iWeights);
        jMI = MI{jSet};
        jWeights = -log10(MIp{jSet});
        if isequal(jWeights,zeros(size(jWeights)))
            jWeights = ones(size(jWeights));
        end
        jWeights = jWeights./sum(jWeights);
        
        % Calculate weighted means
        iwMean = sum(iMI.*iWeights);
        jwMean = sum(jMI.*jWeights);
        
        % Calculate the difference in the weighted means
        realDif = iwMean - jwMean;
        
        % Perform the Monte Carlo on the weighted means
        MCDifs = zeros([nMCcomp,1]);
        dataList = [iMI,iWeights;jMI,jWeights];
        ni = length(iMI);
        nAll = size(dataList,1);
        for iMC = 1:nMCcomp
            dataList = dataList(randperm(nAll),:);
            MCiwMean = sum(dataList(1:ni,1).*dataList(1:ni,2))/sum(dataList(1:ni,2));
            MCjwMean = sum(dataList((ni + 1):end,1).*dataList((ni + 1):end,2))/sum(dataList((ni + 1):end,2));
            MCDifs(iMC) = MCiwMean - MCjwMean;
        end
        if realDif >= 0
            pBS = nnz(MCDifs >= realDif)/nMCcomp;
        else
            pBS = nnz(MCDifs <= realDif)/nMCcomp;
        end
        pBS(pBS == 0) = 1/(2*nMCcomp);
        
        % Record the result
        compEx(iSet,jSet) = pBS;
    end
end


%% Make the figure


% Set parameters for the figure

% Set the side margins in inches
lmargin = 0.45;
rmargin = 0.1;

% Set the top and bottom margins in inches
tmargin = 0.25;
bmargin = 0.4;

% Set the horizontal distance between plots in inches
hspace = 0.5;

% Set the vertical distance between plots in inches
vspace = 0.6;

% Set Paper size
papersize = [7 6];
% papersize = [3 5];

% Set the figure panel dimensions (number of horizontal figures by number
% of vertical figures
figdim = [4,3];

% Set the width of the comparison figure in inches
widthComp = 1.5;

% calculate the width and height based on dimensions above in inches
width = (papersize(1) - lmargin - rmargin - (figdim(1) - 1)*hspace - widthComp)/(figdim(1) - 1);
height = (papersize(2) - tmargin - bmargin - (figdim(2) - 1)*vspace)/figdim(2);

LeftCoord = (0:(figdim(1) - 1))*(width + hspace) + lmargin;
BottomCoord = (0:(figdim(2) - 1))*(height + vspace) + bmargin;
BottomCoord = fliplr(BottomCoord);

% Convert to fraction of page sizes
LeftCoord = LeftCoord/papersize(1);
BottomCoord = BottomCoord/papersize(2);
width = width/papersize(1);
widthComp = widthComp/papersize(1);
height = height/papersize(2);
hspace = hspace/papersize(1);

% Set fontsizes
TitleFS = 8;
AxLabelFS = 7;
UnitFS = 7;
LegFS = 7;
SubFigLabelFS = 14;


% Set the linewidths
ThinLW = 0.6;
MedLW = 0.7;
ThickLW = 1.3;

% Set the dot sizes
DotSz1 = 30;
DotSz2 = 20;

% Set colors
PColor = {'r',[0,210/255,50/255],'b'};

% Figure out the MI, weight, and proportion ranges to unify axis across plots
MILimits = zeros([2,3]);
weightLimits = zeros([2,3]);
propLimits = zeros([2,3]);
for iSet = 1:3
    MItemp = MI{iSet};
    weightstemp = -log10(MIp{iSet});
    if isequal(weightstemp,zeros(size(weightstemp)))
        weightstemp = ones(size(weightstemp));
    end
    weightstemp = weightstemp./sum(weightstemp);
    MILimits(1,iSet) = min(MItemp);
    MILimits(2,iSet) = max(MItemp);
    weightLimits(1,iSet) = min(weightstemp);
    weightLimits(2,iSet) = max(weightstemp);
    
    MItemp = MI{iSet};
    nullMItemp = nullMI{iSet}(:);
    
    xMI = unique(MItemp);
    yMI = histc(MItemp,xMI)./length(MItemp);
    xnullMI = unique(nullMItemp);
    ynullMI = histc(nullMItemp,xnullMI)./length(nullMItemp);
    propLimits(1,iSet) = 0;
    propLimits(2,iSet) = max([max(yMI),max(ynullMI)]);
end

% Make the figure
F1 = figure('PaperPosition',[0 0 papersize],'PaperSize',papersize);
set(gcf, 'PaperUnits', 'inches');

for iSet = 1:3
    
    % Plot the PMFs
    subplot('Position',[LeftCoord(1),BottomCoord(iSet),width,height]);
    hold on
    
    MItemp = MI{iSet};
    nullMItemp = nullMI{iSet}(:);
    
    xMI = unique(MItemp);
    yMI = histc(MItemp,xMI)./length(MItemp);
    xnullMI = unique(nullMItemp);
    ynullMI = histc(nullMItemp,xnullMI)./length(nullMItemp);
    
    plot(xMI,yMI,'Color',PColor{iSet},'LineWidth',MedLW);
    scatter(xMI,yMI,DotSz1,PColor{iSet})
    plot(xnullMI,ynullMI,'k','LineWidth',MedLW);
    scatter(xnullMI,ynullMI,DotSz1,'k','d')
    
    y1 = 1.05*max(propLimits(2,:));
    y2 = 0.95*max(propLimits(2,:));
    x1 = 0.7*max(MILimits(2,:));
    x2 = 0.85*max(MILimits(2,:));
    x3 = 0.87*max(MILimits(2,:));
    
    line([x1,x2],[y1,y1],'Color',PColor{iSet},'LineWidth',MedLW)
    scatter(mean([x1,x2]),y1,DotSz1,PColor{iSet})
    text(x3,y1,'Real','FontSize',LegFS)
    line([x1,x2],[y2,y2],'Color','k','LineWidth',MedLW)
    scatter(mean([x1,x2]),y2,DotSz1,'k','d')
    text(x3,y2,'Null','FontSize',LegFS)
    
    xlim([0,1.1*max(MILimits(2,:))])
    ylim([0,1.1*max(propLimits(2,:))])
    
    title(['Histogram (a = ',num2str(setParams(iSet,4)),')'],'FontSize',TitleFS)
    ylabel('Proportion','FontSize',AxLabelFS)
    xlabel('Mutual Information (bits)','FontSize',AxLabelFS)
    set(gca,'FontSize',UnitFS)
    
    text(-0.1*max(MILimits(2,:)),1.133*max(propLimits(2,:)),['A',num2str(iSet)],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',SubFigLabelFS)
    
    
    
    
    % Plot the CMFs
    subplot('Position',[LeftCoord(2),BottomCoord(iSet),width,height]);
    hold on
    
    MItemp = MI{iSet};
    nullMItemp = nullMI{iSet}(:);
    
    xMI = unique(MItemp);
    yMI = cumsum(histc(MItemp,xMI)./length(MItemp));
    xnullMI = unique(nullMItemp);
    ynullMI = cumsum(histc(nullMItemp,xnullMI)./length(nullMItemp));
    
    plot(xMI,yMI,'Color',PColor{iSet},'LineWidth',MedLW);
    scatter(xMI,yMI,DotSz1,PColor{iSet})
    plot(xnullMI,ynullMI,'k','LineWidth',MedLW);
    scatter(xnullMI,ynullMI,DotSz1,'k','d')
    
    xlim([0,1.1*max(MILimits(2,:))]);
    ylim([0,1]);
    
    y1 = 0.2;
    y2 = 0.1;
    x1 = 0.7*max(MILimits(2,:));
    x2 = 0.85*max(MILimits(2,:));
    x3 = 0.87*max(MILimits(2,:));
    
    line([x1,x2],[y1,y1],'Color',PColor{iSet},'LineWidth',MedLW)
    scatter(mean([x1,x2]),y1,DotSz1,PColor{iSet})
    text(x3,y1,'Real','FontSize',LegFS)
    line([x1,x2],[y2,y2],'Color','k','LineWidth',MedLW)
    scatter(mean([x1,x2]),y2,DotSz1,'k','d')
    text(x3,y2,'Null','FontSize',LegFS)
    
    title('Cumulative Dist.','FontSize',TitleFS)
    ylabel('Cumulative Proportion','FontSize',AxLabelFS)
    xlabel('MI (bits)','FontSize',AxLabelFS)
    set(gca,'FontSize',UnitFS)
    
    text(-0.1*max(MILimits(2,:)),1.03,['B',num2str(iSet)],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',SubFigLabelFS)
    
    
    
    % Plot the mutual information vs. weights
    subplot('Position',[LeftCoord(3),BottomCoord(iSet),width,height]);
    hold on
    
    MItemp = MI{iSet};
    weightstemp = -log10(MIp{iSet});
    if isequal(weightstemp,zeros(size(weightstemp)))
        weightstemp = ones(size(weightstemp));
    end
    weightstemp = weightstemp./sum(weightstemp);
    
    scatter(MItemp,weightstemp,DotSz1,PColor{iSet})
    
    xlim([0,1.1*max(MILimits(2,:))]);
    ylim([0,1.1*max(weightLimits(2,:))]);
    
    title('MI Weights','FontSize',TitleFS)
    ylabel('Normalized Weights','FontSize',AxLabelFS)
    xlabel('MI (bits)','FontSize',AxLabelFS)
    set(gca,'FontSize',UnitFS)
    
    text(-0.1*max(MILimits(2,:)),1.133*max(weightLimits(2,:)),['C',num2str(iSet)],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',SubFigLabelFS)
    
    
    
    % Plot the weighted mutual information with standard error of the
    % weighted mean
    if iSet == 2
        subplot('Position',[LeftCoord(4),BottomCoord(iSet),widthComp,height]);
        hold on
        OffSet = [-0.03,0.03];
        yLimits = zeros([2,3,2]);
        
        for iEx = 1:3
            
            wMean = MIPlot(1,iEx);
            SEwM = MIPlot(2,iEx);
            wSTD = MIPlot(3,iEx);
            line([setParams(iEx,4) + OffSet(1),setParams(iEx,4) + OffSet(1)],[wMean - wSTD,wMean + wSTD],'Color',PColor{iEx},'LineWidth',MedLW)
            line([setParams(iEx,4) + OffSet(1),setParams(iEx,4) + OffSet(1)],[wMean - SEwM,wMean + SEwM],'Color',PColor{iEx},'LineWidth',ThickLW)
            scatter(setParams(iEx,4) + OffSet(1),wMean,DotSz2,PColor{iEx})
            
            yLimits(1,iEx,1) = min([wMean - wSTD,wMean - SEwM]);
            yLimits(2,iEx,1) = max([wMean + wSTD,wMean + SEwM]);
            
            wMean = MIPlot(4,iEx);
            SEwM = MIPlot(5,iEx);
            wSTD = MIPlot(6,iEx);
            line([setParams(iEx,4) + OffSet(2),setParams(iEx,4) + OffSet(2)],[wMean - wSTD,wMean + wSTD],'Color','k','LineWidth',MedLW)
            line([setParams(iEx,4) + OffSet(2),setParams(iEx,4) + OffSet(2)],[wMean - SEwM,wMean + SEwM],'Color','k','LineWidth',ThickLW)
            scatter(setParams(iEx,4) + OffSet(2),wMean,DotSz2,'k','filled')
            
            yLimits(1,iEx,2) = min([wMean - wSTD,wMean - SEwM]);
            yLimits(2,iEx,2) = max([wMean + wSTD,wMean + SEwM]);
            
            if KSp(iEx) >= 0.01
                text(setParams(iEx,4),min(squeeze(yLimits(1,iEx,:))),['p = ',num2str(KSp(iEx),2)],...
                    'HorizontalAlignment','Center','VerticalAlignment','Top','FontSize',LegFS - 1)
            else
                temp = KSp(iEx);
                temp = ceil(log10(temp));
                text(setParams(iEx,4),min(squeeze(yLimits(1,iEx,:))),['p < 10^{',num2str(temp),'}'],...
                    'HorizontalAlignment','Center','VerticalAlignment','Top','FontSize',LegFS - 1)
            end
            
        end
        
        yMax = max(max(squeeze(yLimits(2,:,:))));
        yMin = min(min(squeeze(yLimits(1,:,:))));
        dy = yMax - yMin;
        xlim([0.05,1.15])
        ylim([yMin - 0.15*dy,yMax + 0.3*dy])
        
        for iEx = 1:2
            for jEx = (iEx + 1):3
                if isequal([iEx,jEx],[1,3])
                    bump = 0.19*dy;
                else
                    bump = 0.03*dy;
                end
                line([setParams(iEx,4) + OffSet(1),setParams(jEx,4) + OffSet(1)],[max(max(squeeze(yLimits(2,[iEx,jEx],:)))) + bump,max(max(squeeze(yLimits(2,[iEx,jEx],:)))) + bump],...
                    'LineWidth',ThinLW,'Color','k')
                if compEx(iEx,jEx) >= 0.01
                    text(mean([setParams(iEx,4) + OffSet(1),setParams(jEx,4) + OffSet(1)]),max(max(squeeze(yLimits(2,[iEx,jEx],:)))) + bump + 0.02*dy,['p = ',num2str(compEx(iEx,jEx),2)],...
                        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontSize',LegFS - 1)
                else
                    temp = compEx(iEx,jEx);
                    temp = ceil(log10(temp));
                    text(mean([setParams(iEx,4) + OffSet(1),setParams(jEx,4) + OffSet(1)]),max(max(squeeze(yLimits(2,[iEx,jEx],:)))) + bump + 0.02*dy,['p < 10^{',num2str(temp),'}'],...
                        'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontSize',LegFS - 1)
                end
            end
        end
        
        y1 = yMin + 0.68*dy;
        y2 = yMin + 0.58*dy;
        x1 = 0.63;
        x2 = 0.69;
        x3 = 0.75;
        x4 = 0.8;
        
        scatter(x1,y1,DotSz2,PColor{1})
        scatter(x2,y1,DotSz2,PColor{2})
        scatter(x3,y1,DotSz2,PColor{3})
        text(x4,y1,'Real','FontSize',LegFS)
        scatter(x3,y2,DotSz2,'k','filled')
        text(x4,y2,'Null','FontSize',LegFS)
        
        title('Weighted MI','FontSize',TitleFS)
        ylabel('Weighted MI (bits)','FontSize',AxLabelFS)
        xlabel('Noise (a)','FontSize',AxLabelFS)
        set(gca,'FontSize',UnitFS)
        
        text(-0.06,yMax + 0.3435*dy,'D','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',SubFigLabelFS)
        
    end
    
    
    
end

% Save the figure
print(F1,'-dpdf','-painters','Figure1')
