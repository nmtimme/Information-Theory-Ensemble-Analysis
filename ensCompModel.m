%% ENSCOMMODEL - ensemble information theory comparison model
% Generates mutual information values for two versions of a simple model of
% an ensemble of information sources and compares their weighted mean 
% mutual information. Each information source is modelled as a system with
% two variables, each with two discrete states. High signal strength
% results in strong relationships between the variables and noise
% randomizes interactions.
%
% Syntax: [MI,MIp,nullMI,nullMIp,KSp] = ensCompModel(nEns, nTrials, sigStr, noisStr, nMC, nMCcomp, nModels)
%
% Input:
%   nEns (2 by 1 integer array) - the number of information sources in the
%     ensemble. The first element corresponds to model 1 and the second
%     element corresponds to model 2.
%   nTrials (2 by 1 integer array) - the number of measurement trials in 
%     the model. For simplicity, this must be a multiple of 4. The first 
%     element corresponds to model 1 and the second element corresponds to 
%     model 2.
%   sigStr (2 by 1 real number array) - the strength of the signal between 
%     the two variables in each information source. The number is greater 
%     than or equal to 0 and less than or equal to 1. The first element 
%     corresponds to model 1 and the second element corresponds to model 2.
%   noisStr (2 by 1 real number array) - the strength of the noise in the 
%     interaction between the two variables in each information source. The
%     number is greater than or equal to 0 and less than or equal to 1. The 
%     first element corresponds to model 1 and the second element 
%     corresponds to model 2.
%   nMC (2 by 1 integer array) - the number of Monte Carlo randomized data 
%     trials used to estimate the p-value of each information source value.
%     The first element corresponds to model 1 and the second element 
%     corresponds to model 2.
%   nMCcomp (integer) - the number of Monte Carlo randomized data trials
%     used to estimate the p-value for the comparison between the weighted
%     means of the two models.
%   nModels (integer) - the number of pairs of models (each with nEns 
%     sources) to run.
%
% Outputs:
%   MI (nEns by nModels double array) - mutual information values in bits
%     for each member of the ensembles of each model.
%   MIp (nEns by nModels double array) - the estimated p-value for each
%     corresponding mutual information value in MI computed via Monte Carlo
%     with randomized data.
%   nullMI (nEns by nMC by nModels double array) - the mutual information
%     values found in the nMC randomized trials used compute the p-value.
%   nullMIp (nEns by nMC by nModels double array) - the p-values of the
%     randomized trials.
%   KSp (1 by nModels double array) - the p-value from the KS-test between
%     the distribution of real MI values for the ensemble members and the
%     distribution of MI values for the randomized data.
%
%
% Other m-files required: ensModel
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% February 2017; Last revision: 5-Jul-2018


function [compp] = ensCompModel(nEns, nTrials, sigStr, noisStr, nMC, nMCcomp, nModels)

% Preallocate space
compp = NaN([nModels,1]);

for iModel = 1:nModels
    
    % Run the models
    [iwMean,~,~,iMI,iMIp] = ensModel(nEns(1), nTrials(1), sigStr(1), noisStr(1), nMC(1), 1);
    [jwMean,~,~,jMI,jMIp] = ensModel(nEns(2), nTrials(2), sigStr(2), noisStr(2), nMC(2), 1);
    
    % Organize the weights
    iWeights = -log10(iMIp);
    if isequal(iWeights,zeros(size(iWeights)))
        iWeights = ones(size(iWeights));
    end
    iWeights = iWeights./sum(iWeights);
    jWeights = -log10(jMIp);
    if isequal(jWeights,zeros(size(jWeights)))
        jWeights = ones(size(jWeights));
    end
    jWeights = jWeights./sum(jWeights);
    
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
    compp(iModel) = pBS;
    
end

