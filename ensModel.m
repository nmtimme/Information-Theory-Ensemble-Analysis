%% ENSMODEL - ensemble information theory model
% Generates mutual information values for a simple model of an ensemble of
% information sources. Each information source is modelled as a system with
% two variables, each with two discrete states. High signal strength
% results in strong relationships between the variables and noise
% randomizes interactions.
%
% Syntax: [wmMI,sewmMI,wsdMI,MI,MIp,nullMI,nullMIp,KSp] = ensModel(nEns, nObs, sigStr, noisStr, nMC, nModelSets)
%
% Input:
%   nEns (integer) - the number of information sources in the ensemble.
%   nObs (integer) - the number of measurement observations in the model. 
%     For simplicity, this must be a multiple of 4.
%   sigStr (real number) - the strength of the signal between the two
%     variables in each information source. The number is greater than or
%     equal to 0 and less than or equal to 1. This variable is refered to
%     as s in the manuscript.
%   noisStr (real number) - the strength of the noise in the interaction
%     between the two variables in each information source. The number is
%     greater than or equal to 0 and less than or equal to 1. This variable
%     is refered to as a in the manuscript.
%   nMC (integer) - the number of Monte Carlo randomized data trials used
%     to estimate the p-value of each information source value.
%   nModelSets (integer) - the number of model ensembles (each with nEns 
%     sources) to produce.
%
% Outputs:
%   wmMI (1 by nModelSets double array) - weighted mean of the mutual
%     information values produced by the ensemble in bits.
%   sewmMI (1 by nModelSets double array) - standard error of the weighted
%     mean of the mutual information values produced by the ensemble in
%     bits.
%   wsdMI (1 by nModelSets double array) - weighted standard deviation of
%     mutual information values produced by the ensemble in bits.
%   MI (nEns by nModelSets double array) - mutual information values in 
%     bits for each member of the ensembles of each model.
%   MIp (nEns by nModelSets double array) - the estimated p-value for each
%     corresponding mutual information value in MI computed via Monte Carlo
%     with randomized data.
%   nullMI (nEns by nMC by nModelSets double array) - the mutual 
%     information values found in the nMC randomized trials used compute 
%     the p-value.
%   nullMIp (nEns by nMC by nModelSets double array) - the p-values of the
%     randomized trials.
%   KSp (1 by nModelSets double array) - the p-value from the KS-test 
%     between the distribution of real MI values for the ensemble members 
%     and the distribution of MI values for the randomized data.
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Nick Timme
% Email: nicholas.m.timme@gmail.com
% February 2017; Last revision: 26-Apr-2018


function [wmMI,sewmMI,wsdMI,MI,MIp,nullMI,nullMIp,KSp] = ensModel(nEns, nObs, sigStr, noisStr, nMC, nModelSets)

% Error check the number of trials
if rem(nObs,4) ~= 0
    error('nObs must be a multiple of 4.')
end
% Error check the signal strength
if (sigStr > 1) || (sigStr < 0)
    error('sigStr must be greater than or equal to 0 and less than or equal to 1.')
end
% Error check the noise strength
if (noisStr > 1) || (noisStr < 0)
    error('noisStr must be greater than or equal to 0 and less than or equal to 1.')
end

% Calculate a rounding factor to correct for system resolution
rFact = 10^ceil(log10(eps) + 2);

% Make the states based on the signal strength
xTrialStates = [ones([1,nObs/2]),2*ones([1,nObs/2])];
nHigh = (nObs/4) + round((nObs/4)*sigStr);
nLow = (nObs/4) - round((nObs/4)*sigStr);
yTrialStates = [ones([1,nHigh]),2*ones([1,nLow]),ones([1,nLow]),2*ones([1,nHigh])];

% Calculate the number of random flips due to noise
nNois = round(nObs*noisStr);

% Preallocate space
MI = NaN([nEns,nModelSets]);
nullMI = NaN([nEns,nMC,nModelSets]);
MIp = NaN([nEns,nModelSets]);
nullMIp = NaN([nEns,nMC,nModelSets]);
wmMI = NaN([1,nModelSets]);
sewmMI = NaN([1,nModelSets]);
wsdMI = NaN([1,nModelSets]);

% Generate the data and calculate the mutual information
for iModelSet = 1:nModelSets
    for iEns = 1:nEns
        
        % Make the counts matrix by adding noise to the general
        % states
        randlist = randperm(nObs);
        yTrialRandStates = yTrialStates;
        yTrialRandStates(randlist(1:nNois)) = yTrialRandStates(randlist(randperm(nNois)));
        counts = accumarray({xTrialStates,yTrialRandStates},ones([1,nObs]),[2,2]);
        
        % Convert joint state counts to a probability distribution
        p = counts./nObs;
        
        % Calculate MI
        px = sum(p,2);
        px = px(:,[1,1]);
        py = sum(p,1);
        py = py([1,1],:);
        MITerms = p.*log2(p./(px.*py));
        MITerms(~isfinite(MITerms)) = 0;
        MI(iEns,iModelSet) = sum(MITerms(:));
        
        % Run the Monte Carlo.
        xTrialNullStates = [ones([1,nObs/2]),2*ones([1,nObs/2])];
        for iMC = 1:nMC
            
            % Make a randomized version of one of the variables
            yTrialNullStates = [ones([1,nObs/2]),2*ones([1,nObs/2])];
            yTrialNullStates = yTrialNullStates(randperm(nObs));
            
            % Count the number of joint states.
            counts = accumarray({xTrialNullStates,yTrialNullStates},ones([1,nObs]),[2,2]);
            
            % Convert to a probability distribution.
            p = counts./nObs;
            
            % Calculate mutual information.
            px = sum(p,2);
            px = px(:,[1,1]);
            py = sum(p,1);
            py = py([1,1],:);
            MITerms = p.*log2(p./(px.*py));
            MITerms(~isfinite(MITerms)) = 0;
            nullMI(iEns,iMC,iModelSet) = sum(MITerms(:));
        end
        
        % Calculate the p-value for the real data
        MIp(iEns,iModelSet) = nnz(nullMI(iEns,:,iModelSet) >= MI(iEns,iModelSet))/nMC;
        
        % Calculate the p-value for the null data
        for iMC = 1:nMC
            nullMIp(iEns,iMC,iModelSet) = nnz(nullMI(iEns,:,iModelSet) >= nullMI(iEns,iMC,iModelSet))/nMC;
        end
        
    end
    
    % Correct zero p-values based on resolution from Monte Carlo
    MIp(MIp(:,iModelSet) == 0,iModelSet) = 1/(2*nMC);
    
    % Calculate the weights based on the p-values
    weights = -log10(MIp(:,iModelSet));
    if isequal(weights,zeros(size(weights)))
        weights = ones(size(weights));
    end
    weights = weights./sum(weights);
    
    % Calculate the weighted mean, standard error of the weighted mean,
    % and weighted standard deviation
    wmMI(iModelSet) = sum(MI(:,iModelSet).*weights);
    sewmMI(iModelSet) = std(MI(:,iModelSet))*sqrt(sum(weights.^2));
    wsdMI(iModelSet) = sqrt(sum(weights.*(MI(:,iModelSet) - sum(MI(:,iModelSet).*weights)).^2));
    
end

% Correct zero p-values based on resolution from Monte Carlo
nullMIp(nullMIp == 0) = 1/(2*nMC);

% Correct system resolution rounding errors
MI = round(MI./rFact).*rFact;
nullMI = round(nullMI./rFact).*rFact;

% Perform KS tests between real and null MI distributions
KSp = NaN([1,nModelSets]);
for iModelSet = 1:nModelSets
    nullTemp = nullMI(:,:,iModelSet);
    [~,KSp(iModelSet)] = kstest2(MI(:,iModelSet),nullTemp(:));
end

end