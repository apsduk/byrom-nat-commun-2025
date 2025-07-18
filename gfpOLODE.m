 
%% INFO %%
% This function generates a set of ODEs and other outputs for a host
% containing the open-loop circuit.
%
% Inputs:
%   Y        - Vector of initial conditions for the system
%   phost    - Vector containing all the fixed host parameters
%   pcircuit - Matrix containing circuit parameters for each subpopulation.
%              Each column represents one subpopulation.
%   mutrates - Matrix defining the rates of mutation between subpopulations
%   culture  - String determining the culture as one of batch,
%              chemostat, fed-batch or (population-free) mid-exponential
% Outputs:
%   dY      - Vector of differential equations, used for simulation with
%             ode15s.
%   lambda  - Matrix of growth rates showing the growth rates of each
%             subpopulation over time
%   mass    - Matrix of the masses of different protein types over time
%   gammaX  - Elongation rate
%   txnrates- Vector of transcription rates for different genes.
%   tlnrates- Vector of translation rates for different genes.
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

function [dY, lambda, mass, gammaX, txnrates, tlnrates] = gfpOLODE(Y, phost, pcircuit, mutrates, culture)

%% VARIABLE DEFINITION %%

% Retrieve number of competing subpopulations
n = size(mutrates,1);

% Split Y vector into matrix of subpopulations
YY = reshape(Y,[],n);

% Shared external substrate
xS = YY( 1,1);

% Define variable vectors
N  = YY( 2,:); iS = YY( 3,:); ee = YY( 4,:);
mT = YY( 5,:); cT = YY( 6,:); pT = YY( 7,:);
mE = YY( 8,:); cE = YY( 9,:); pE = YY(10,:);
mH = YY(11,:); cH = YY(12,:); pH = YY(13,:);
mR = YY(14,:); cR = YY(15,:); pR = YY(16,:); r = YY(17,:); R = YY(18,:);
% Circuit:
mA = YY(19,:); cA = YY(20,:); pA = YY(21,:);

% Define parameters
phie0   = phost(01); dyN  = phost(02); influx = phost(03); dilution = phost(04);
   vT   = phost(05); vE   = phost(06); kT   = phost(07); kE   = phost(08);
   wT   = phost(09); wE   = phost(10); wH   = phost(11); wR   = phost(12); wr = phost(13);
   oT   = phost(14); oE   = phost(15); oH   = phost(16); oR   = phost(17); or = phost(18);
   nT   = phost(19); nE   = phost(20); nH   = phost(21); nR   = phost(22);
   bT   = phost(23); bE   = phost(24); bH   = phost(25); bR   = phost(26);
   uT   = phost(27); uE   = phost(28); uH   = phost(29); uR   = phost(30); 
   dymT = phost(31); dymE = phost(32); dymH = phost(33); dymR = phost(34); dyrr = phost(35);
   brho = phost(36); urho = phost(37);
   maxG = phost(38); kG   = phost(39); M0   = phost(40); kH   = phost(41); hH = phost(42);

wA = pcircuit(1,:); oA = pcircuit(2,:); nA = pcircuit(3,:);
bA = pcircuit(4,:); uA = pcircuit(5,:); dymA = pcircuit(6,:); dypA = pcircuit(7,:);

%% RATE CALCULATION %%

% Elongation rate
gammaX = (maxG*ee)./(kG+ee);

% Growth rate
cX = cA;
lambda = (1/M0).*gammaX.*sum([cT;cE;cH;cR;cX]);

% Transcription rates
g2mT = wT*ee./(oT + ee);
g2mE = wE*ee./(oE + ee);
g2mH = (wH*ee./(oH + ee))./(1 + (pH./kH).^hH);
g2mR = wR*ee./(oR + ee);
g2rr = wr*ee./(or + ee);

% Translation rates
m2pT = (gammaX./nT).*cT;
m2pE = (gammaX./nE).*cE;
m2pH = (gammaX./nH).*cH;
m2pR = (gammaX./nR).*cR;

% Substrate import and catabolism
x2iS = pT.*vT.*xS./(kT + xS);
iS2e = pE.*vE.*iS./(kE + iS);

% Circuit rates
g2mA = (wA.*ee)./(oA + ee);
m2pA = (gammaX./nA).*cA;


%% ODEs %%

% Circuit ODEs
dmA = g2mA + m2pA - bA.*R.*mA + uA.*cA - (lambda + dymA).*mA;
dcA = bA.*R.*mA - uA.*cA - m2pA - lambda.*cA;
dpA = m2pA - (lambda + dypA).*pA;

% Circuit influence on host ODEs
iSX = 0;
eeX = - nA.*m2pA;
RX = m2pA - bA.*R.*mA + uA.*cA;

% Substrate import and population growth, based on culture type
% By default we assume midexp, unless a different culture is specified

switch culture
    case "midexp"
        dxS = 0;
        dN = zeros(1,n);  
    case "batch"
        dxS = -dot(x2iS,N);
        dN = ((mutrates' + diag(lambda' - sum(mutrates,2) - dyN))*N')';
    case "chemostat"
        dxS = -dot(x2iS,N) + influx - dilution*xS;
        dN = ((mutrates' + diag(lambda' - sum(mutrates,2) - dilution))*N')';
    case "fedbatch"
        dxS = -dot(x2iS,N) + influx;
        dN = ((mutrates' + diag(lambda' - sum(mutrates,2) - dyN))*N')';
end

% Internal substrate and energy
diS = x2iS - iS2e - lambda.*iS + iSX;
dee = phie0 * iS2e + eeX - lambda.*ee - (nT.*m2pT + nE.*m2pE + nH.*m2pH + nR.*m2pR);

% Transport
dmT = g2mT + m2pT - bT.*R.*mT + uT.*cT - mT.*(lambda + dymT);
dcT = bT.*R.*mT - uT.*cT - m2pT - lambda.*cT;
dpT = m2pT - lambda.*pT;

% Enzymes
dmE = g2mE + m2pE - bE.*R.*mE + uE.*cE - mE.*(lambda + dymE);
dcE = bE.*R.*mE - uE.*cE - m2pE - lambda.*cE;
dpE = m2pE - lambda.*pE;

% Housekeeping proteins
dmH = g2mH + m2pH - bH.*R.*mH + uH.*cH - mH.*(lambda + dymH);
dcH = bH.*R.*mH - uH.*cH - m2pH - lambda.*cH;
dpH = m2pH - lambda.*pH;

% Ribosomes
dmR = g2mR + m2pR - bR.*R.*mR + uR.*cR - (lambda + dymR).*mR;
dcR = bR.*R.*mR - uR.*cR - m2pR - lambda.*cR;
dpR = m2pR - lambda.*pR - brho.*pR.*r + urho.*R;
dr = g2rr - brho.*pR.*r + urho.*R - (lambda + dyrr).*r;
dR = brho.*pR.*r - urho.*R - lambda.*R ...
    + m2pT - bT.*R.*mT + uT.*cT ...
    + m2pE - bE.*R.*mE + uE.*cE ...
    + m2pH - bH.*R.*mH + uH.*cH ...
    + m2pR - bR.*R.*mR + uR.*cR ...
    + RX;

%% OUTPUTS %%

dYY = [dxS zeros(1,n-1); dN; diS; dee; dmT; dcT; dpT; dmE; dcE; dpE; dmH; dcH; dpH; dmR; dcR; dpR; dr; dR; dmA; dcA; dpA];
dY = reshape(dYY,[],1);
txnrates = [g2mT; g2mE; g2mH; g2mR; g2rr; g2mA];
tlnrates = [m2pT; m2pE; m2pH; m2pR; m2pA];
mass = [nT*pT; nE*pE; nH*pH; nR*(pR+R+cT+cE+cH+cR+cA); nA.*pA];
    
end





