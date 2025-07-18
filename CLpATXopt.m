%%%%% CLpATXvsOL %%%%%
% This script performs a multi-objective optimisation for the CLpATX
% system, aiming to simultaneously maximise initial output P0, 'half-life'
% T50, and time taken to fall outside a window of +-10% of the original
% value T+-10
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

%% Host setup

% Clear and close everything
clc; close all; clear;

% Define culture
culture = "batch";

% Load parameters and guessed initial conditions for the host
[ph,Y0] = loadhostparameters();

% Set ODE options
odeoptions = odeset('NonNegative', 1:numel(Y0));

% Simulate single host in mid-exp phase to get steady-state initial conditions
[T0, Y0] = ode15s( @(T,Y) hostonlyODE(Y,ph,"midexp"), [0, 1e12],Y0,odeoptions);

%% Circuit setup

% Enter circuit filenames
circuit = "gfpCLpATX";
% Enter default base circuit parameterss 
wA = 100; oA = 4.38; nA = 300;
bA = 0.1; uA = 0.01; dymA = 0.1; dypA = 0;
nB = 300; bB = 0.001; uB = 0.01; dypB = 0;
kB = 1e5; hB = 2;

p1 = [wA; oA; nA; bA; uA; dymA; dypA; nB; bB; uB; dypB; kB; hB];

% Define variables and initial values
mA = 0; cA = 0; pA = 0;
        cB = 0; pB = 0;

% Update initial conditions vector
Y0 = [Y0(end,:)'; mA; cA; pA; cB; pB];

%% Mutation setup

% Define number of mutation states
M = 1;
N = 4;
n = N^M;

% Define global mutation rate
sigma = 1e-5;

% Define mutation rate matrix
mutrates = generatemutrates(M,N,sigma);

% Set up percentage values of mutatable parameters for each mutation state
wApercent = repelem(linspace(1,0,N),n/N);
wXpercent = wApercent;

% Define parameters for the different subpopulations
pc = repmat(p1,1,n);

% Set the wA values for each state
pc(1,:) = pc(1,:).*wApercent;

% We add a dummy parameter for the sake of the efast modelling
pc = [pc;zeros(1,size(pc,2))];

% Build initial conditions matrix
Y0 = repmat(Y0,1,n);

% External substrate is shared, so we remove all but one of these
% variables. We also set all mutated population sizes to zero.
Y0(1:2,2:end) = 0;

% Other external variables are also shared, such as quorum sensing
% molecules
externalvarIDs = [];
Y0(externalvarIDs,2:end) = 0;

%% ===== SET UP OPTIMISATION OPTIONS ======================================
% --- Set up unknown parameters with initial conditions from efast Pareto front


knames{ 1} = 'wA';   lb( 1) = 1e-1; ub( 1) = 1e3; loggy( 1) = 1; idx{ 1} =  1;
knames{ 2} = 'bB';   lb( 2) = 1e-3; ub( 2) = 1;   loggy( 2) = 1; idx{ 2} =  9;
knames{ 3} = 'kB';   lb( 3) = 1e3;  ub( 3) = 1e6; loggy( 3) = 1; idx{ 3} =  12;

for i = 1:numel(lb)
    if loggy(i) == 1
        lb(i) = log10(lb(i));
        ub(i) = log10(ub(i));
    end
end

popsize = 250;

objfun = @(X) gfpcost(X, idx, loggy, circuit, ph, pc, externalvarIDs, mutrates, wXpercent,Y0);

% --- Multi GA options ----------------------------------------------------
multigaoptions = optimoptions(@gamultiobj, 'PopulationSize', popsize, 'UseParallel', true, 'FunctionTolerance', 2e-5, 'ParetoFraction', 0.4, ...
    'MaxGenerations', 1000, 'MaxStallGenerations', 100, 'Display', 'iter', 'MaxTime', 60*60*24*5);

%% ===== CARRY OUT MULTI-OBJECTIVE OPTIMISATION ===========================
[bestX, fval, exitflag, outputstruct, population, scores] = gamultiobj(objfun, length(lb), [], [], [], [], lb, ub, multigaoptions);
for i = 1:numel(lb)
    if loggy(i) == 1
        bestX(:,i) = 10.^(bestX(:,i));
    end
end

%% ===== SAVE RESULTS =====================================================
fname = mknewdir('gfpCLpATXopt');
save(fname+"/results.mat","bestX","fval","exitflag","outputstruct","population","scores");
