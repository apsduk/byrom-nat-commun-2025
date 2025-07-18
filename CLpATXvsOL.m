%%%%% CLpATXvsOL %%%%%
% This script compares the time series output of two systems which each
% initially produce the same amount of output P0: an open-loop system (OL) and
% an intracircuit control system that acts using a transcription factor
% (CLpATX).
% The user can input the value of P0 they want to consider, and the program
% will look at the outputs from the multiobjective optimisation to find the
% system parameters for the closed-loop system that most closely align with
% that inputted value. The program then tweaks the maximal transcription rate wA
% so that P0 perfectly matches the user input, and similarly varies wA in
% the open-loop system to find a system of equal P0.
% Once the parameterisations are found, the full time series is simulated
% and plotted for the two systems.
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
circuitOL = "gfpOL";
circuitX = "gfpCLpATX";
% Enter default base circuit parameterss 
wA = 100; oA = 4.38; nA = 300;
bA = 0.1; uA = 0.01; dymA = 0.1; dypA = 0;
nB = 1; bB = 0.1; uB = 0.01; dypB = 0;
kB = 1e5; hB = 2;

pOL = [wA; oA; nA; bA; uA; dymA; dypA];

% Enter value of P0 to consider
P0 = 9.5;

% Find optimised parameters with initial output closest to P0.
load CLpATXOptimisations.mat
[~, I] = min(abs(fval(:,1)+P0));
P0XnB300 = -fval(I,1);
wA = bestX(I,1);
bB = bestX(I,2);
kB = bestX(I,3);
pXnB300 = [wA; oA; nA; bA; uA; dymA; dypA; 300; bB; uB; dypB; kB; hB];

% Define variables and initial values
mA = 0; cA = 0; pA = 0;
        cB = 0; pB = 0;

% Update initial conditions vector
Y0OL = [Y0(end,:)'; mA; cA; pA];
Y0X = [Y0(end,:)'; mA; cA; pA; cB; pB];

%% Mutation setup

% Define number of mutation states
MX = 1;
N = 4;
nX = N^MX;

% Define global mutation rate
sigma = 1e-5;

% Define mutation rate matrix
mutratesX = generatemutrates(MX,N,sigma);

% Set up percentage values of mutatable parameters for each mutation state
wApercentX = repelem(linspace(1,0,N),nX/N);

% Define parameters for the different subpopulations
pcOL = repmat(pOL,1,nX);
pcXnB300 = repmat(pXnB300,1,nX);

% Set the wA values for each state
pcOL(1,:) = pcOL(1,:).*wApercentX;
pcXnB300(1,:) = pcXnB300(1,:).*wApercentX;


% Build initial conditions matrix
Y0OL = repmat(Y0OL,1,nX);
Y0X = repmat(Y0X,1,nX);

% External substrate is shared, so we remove all but one of these
% variables. We also set all mutated population sizes to zero.
Y0OL(1:2,2:end) = 0;
Y0X(1:2,2:end) = 0;

% Other external variables are also shared, such as quorum sensing
% molecules
externalvarIDs = [];

%% Batch setup

% Set up parameters for repeated batch simulations
maxruntime = 1000*24*60;
feedtime = 24*60;
feedsize = Y0X(1);
startpopsize = Y0X(2);
stoppercentage = 0.005;

%% Run

% For the open-loop system, repeatedly modify the value of wA until
% initial output is close enough to P0.
[~,~,~,~,~,pOutday0] = repeatbatch(feedtime*3,feedtime,feedsize,startpopsize,stoppercentage,circuitOL,ph,pcOL,externalvarIDs,zeros(size(mutratesX)),Y0OL);
P0OL = log10(pOutday0(end));
while abs(P0OL - P0) > 0.001    
    pcOL(1,:) = pcOL(1,:).*exp(P0)./exp(P0OL);
    [~,~,Yy,~,~,pOutday0] = repeatbatch(feedtime*3,feedtime,feedsize,startpopsize,stoppercentage,circuitOL,ph,pcOL,externalvarIDs,zeros(size(mutratesX)),Y0OL);
    P0OL = log10(pOutday0(end));
end
Y0OL = reshape(Yy(end,:),[],nX);
% Then simulate full time series.
[TOL,TdayOL,YOL,YdayOL,pOutOL,pOutdayOL] = repeatbatch(maxruntime,feedtime,feedsize,startpopsize,stoppercentage,circuitOL,ph,pcOL,externalvarIDs,mutratesX,Y0OL);

% For the open-loop system, repeatedly modify the value of wA until
% initial output is close enough to P0.
while abs(P0XnB300 - P0) > 0.001    
    pcXnB300(1,:) = pcXnB300(1,:).*exp(P0)./exp(P0XnB300);
    [~,~,Yy,~,~,pOutday0] = repeatbatch(feedtime*3,feedtime,feedsize,startpopsize,stoppercentage,circuitX,ph,pcXnB300,externalvarIDs,zeros(size(mutratesX)),Y0X);
    P0XnB300 = log10(pOutday0(end));
end
Y0X = reshape(Yy(end,:),[],nX);
% Then simulate full time series.
[TXnB300,TdayXnB300,YXnB300,YdayXnB300,pOutXnB300,pOutdayXnB300] = repeatbatch(maxruntime,feedtime,feedsize,startpopsize,stoppercentage,circuitX,ph,pcXnB300,externalvarIDs,mutratesX,Y0X);


%% Output

set(0,'defaultaxesfontsize',20)
set(0,'defaultlinelinewidth',3)

blue = [70 150 220]/256;

tiledlayout(1,1,"TileSpacing","tight","Padding","compact")

nexttile
box
hold on
plot(TdayOL/1440,pOutdayOL,'k:')
plot(TdayXnB300/1440,pOutdayXnB300,'color',blue)
yline(1.1*10^9.5,'--','LineWidth',2)
yline(0.9*10^9.5,'--','LineWidth',2)
yline(0.55*10^9.5,'--','LineWidth',2)
xlim([0 50])
leg = legend('OL','CLp_ATX');
xlabel('Time (days)')
ylabel('P (molecules)')







