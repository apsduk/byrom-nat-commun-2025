function [T,Tday,Y,Yday,pOut,pOutday] = repeatbatch(maxruntime,feedtime,feedsize,startpopsize,stoppercentage,circuit,phost,pcircuit,externalvarIDs,mutrates,Y0)
% This function runs a repeated batch simulation. An ODE model is simulated
% multiple times successively. Between each simulation, population size is
% reset and nutrient supply is reset.
% Inputs:
%   maxruntime:     Upper limit for simulation time
%   feedtime:       Time between successive batch simulations
%   feedsize:       Amount of nutrients xS to provide at the start of
%                   each batch simulation
%   startpopsize:   Population size to reset to at the start of each
%                   batch simulation
%   stoppercentage: Percentage of the population-wide output P after
%                   which the simulation ends
%   circuit:        String defining the circuit model being simulated
%   phost:          Vector of host parameters
%   pcircuit:       Matrix of circuit parameters - each column is a
%                   different subpopulation
%   externalvarIDs: IDs of variables which are population-level rather than
%                   per-cell (e.g. external quorum sensing molecule)
%   mutrates:       Matrix of transition rates between subpopulations
%   Y0:             Matrix of initial conditions - each column is a
%                   different subpopulation
% Outputs:
%   T:              Time series vector for the simulation
%   Tday:           Time series vector taking one point at the end of each
%                   simulation day
%   Y:              Matrix of variable outputs over time
%   Yday:           Matrix of variable outputs taken once at the end of
%                   each simulation day
%   pOut:           Population-wide circuit output P over time
%   pOutday:        Population-wide circuit output P taken once at the end
%                   of each simulation day
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

% Define circuit
circuitODE = str2func(circuit + "ODE");
circuitJac = str2func(circuit + "Jbatch");

% Retrieve number of variables and subpopulations
[vars,n] = size(Y0);

% Set up odeoptions
odeoptions = odeset('NonNegative',1:numel(Y0),'Jacobian',@(T,Y) circuitJac(Y,phost,pcircuit,mutrates));

% Track number of feedings
feedings = ceil(maxruntime/feedtime);
feedingscompleted = feedings;

% Set up cells for outputs of each repeated simulation
Tx = cell(feedings,1); Yx = cell(feedings,1); pOutx = cell(feedings,1);
Yday = zeros(feedings+1,vars,n); pOutday = zeros(feedings+1,1);
Yday(1,:,:) = Y0; pOutday(1) = sum(Y0(2,:).*Y0(21,:));

% For loop for repeated simulations
for i = 1:feedings
    % Set Y0 internal variables from previous cycle
    Y0(3:end,:) = Yday(i,3:end,:);

    % Replenish food to feedsize
    Y0(1,1) = feedsize;

    % Scale total population size to startpopsize, portioned according to
    % end of previous cycle, as well as other external variables
    Y0(2,:) = Yday(i,2,:).*startpopsize./sum(Yday(i,2,:));
    Y0(externalvarIDs,:) = Yday(i,externalvarIDs,:).*startpopsize./sum(Yday(i,2,:));

    % Simulate the system up to feedtime
    [Tx{i},Yy] = ode15s( @(T,Y) circuitODE(Y, phost, pcircuit, mutrates, "batch"),[feedtime*(i-1), feedtime*i-1],reshape(Y0,[],1),odeoptions);
    Yx{i} = reshape(Yy,numel(Tx{i}),[],n);
    pOutx{i} = sum(Yx{i}(:,2,:).*Yx{i}(:,21,:),3);

    % Calculate the variable outputs at the end of the cycle
    Yday(i+1,:,:) = Yx{i}(end,:,:);
    pOutday(i+1) = pOutx{i}(end);

    % If pOut is <stoppercentage% of initial, end the simulation
    if pOutday(i+1) < stoppercentage*pOutday(1)
        feedingscompleted = i;
        break
    end
end
% Produce final output vectors
T = cat(1,Tx{:}); Y = cat(1,Yx{:}); pOut = cat(1,pOutx{:});
Y = reshape(Y,numel(T),[],n);

Tday = [0; (feedtime-1:feedtime:T(end))'];
Yday = Yday(1:feedingscompleted+1,:,:);
pOutday = pOutday(1:feedingscompleted+1);

end