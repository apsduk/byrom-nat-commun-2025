%% ===== OUTPUT FUNCTION :: SIMULATEMODEL =================================
function cost = gfpcost(X, idx, loggy, circuit, ph, pc, externalvarIDs, mutrates, wXpercent, Y0)
% This function is used as a cost function for multiopjective optimisations
% which aim to simultaneously maximise P0, T+-10 and T50. It is used with
% the built-in matlab function gamultiobj.
% Inputs:
%   X        - Vector of parameters to be varied through the optimisation
%   idx      - IDs of the rows representing varying parameters in the
%              circuit parameters matrix pc
%   loggy    - Binary vector indicating whether elements of X should be
%              varied on a log scale
%   circuit  - String defining the name of the model to be optimised
%   ph       - Vector of host parameters
%   pc       - Matrix of circuit parameters - each column represents a
%              different subpopulation.
%   externalvarIDs - IDs of any circuit variables to be treated as
%              population-level rather than per-cell (e.g. external quorum sensing
%              molecule)
%   mutrates - Matrix of transition rates between subpopulations.
%   wXpercent- Matrix defining the percentage levels of function of
%              variables which mutate
%   Y0       - Initial conditions matrix for variables.
% Outputs:
%   Cost function - a vector of metrics to be minimised. (Set as negative
%                   because we want to maximise P0, T50, T+-10)
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

%no_of_mutating_variables
M = size(wXpercent,1);

% Assign parameter values from the vector X.
for i = 1:M
    pc(idx{i},:) = 10.^(X(i)).*wXpercent(i,:);
end
for i = M+1:numel(idx)
    if loggy(i) == 1
        pc(idx{i},:) = 10.^(X(i));
    else
        pc(idx{i},:) = X(i);
    end
end

% Batch setup
maxruntime = 1000*24*60;
feedtime = 24*60;
feedsize = Y0(1);
startpopsize = Y0(2);

n = size(mutrates,1);

% --- Simulate model ------------------------------------------------------
% Reach steady state without mutation over 3 days
    [~,~,Yy,~,~,pOutday0] = repeatbatch(feedtime*3,feedtime,feedsize,startpopsize,0,circuit,ph,pc,externalvarIDs,zeros(size(mutrates)),Y0);
    Y00 = reshape(Yy(end,:),[],n);
    pOutday0 = pOutday0(end);
    
    % Run the full simulation and calcuclate xlifes if the initial output is greater than 1e7
    if pOutday0 >1e7        
        [~,Tday,~,~,~,pOutday] = repeatbatch(maxruntime,feedtime,feedsize,startpopsize,0.49,circuit,ph,pc,externalvarIDs,mutrates,Y00);   
        life1 = xlifecalculator(0.5,Tday,pOutday)./(24*60);
        life2 = xlifecalculatoreitherside(0.9,Tday,pOutday)./(24*60);
    % Otherwise just set lifes to 0 so that they won't be selected for in
    % the optimisation
    else
        life1 = 0;
        life2 = 0;
    end
cost = -[log10(pOutday0) life1 life2];
end