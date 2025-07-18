 
%% INFO %%
% This function generates a set of ODEs and other outputs for a host
% containing no circuits.
%
% Inputs:
%   Y        - Vector of initial conditions for the system
%   phost    - Vector containing all the fixed host parameters
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

function [dY, lambda, mass, gammaX, txnrates, tlnrates] = hostonlyODE(Y, phost, culture)

%% VARIABLE DEFINITION %%

% Shared external substrate
xS = Y( 1);

% Define variable vectors
N  = Y( 2); iS = Y( 3); ee = Y( 4);
mT = Y( 5); cT = Y( 6); pT = Y( 7);
mE = Y( 8); cE = Y( 9); pE = Y(10);
mH = Y(11); cH = Y(12); pH = Y(13);
mR = Y(14); cR = Y(15); pR = Y(16); r = Y(17); R = Y(18);

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

%% RATE CALCULATION %%

% Elongation rate
gammaX = (maxG*ee)./(kG+ee);

% Growth rate
lambda = (1/M0).*gammaX.*sum([cT;cE;cH;cR]);

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


%% ODEs %%

switch culture
    case "midexp"
        dxS = 0;
        dN = 0;  
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
diS = x2iS - iS2e - lambda.*iS;
dee = phie0 * iS2e - lambda.*ee - (nT.*m2pT + nE.*m2pE + nH.*m2pH + nR.*m2pR);

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
    + m2pR - bR.*R.*mR + uR.*cR;

%% OUTPUTS %%

dY = [dxS; dN; diS; dee; dmT; dcT; dpT; dmE; dcE; dpE; dmH; dcH; dpH; dmR; dcR; dpR; dr; dR];
dY = reshape(dY,[],1);
txnrates = [g2mT; g2mE; g2mH; g2mR; g2rr];
tlnrates = [m2pT; m2pE; m2pH; m2pR];
mass = [nT*pT; nE*pE; nH*pH; nR*(pR+R+cT+cE+cH+cR)];
    
end





