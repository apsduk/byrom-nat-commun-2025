function lifes = xlifecalculatoreitherside(x,T,P)
% This function calculates the times taken lifes for a decaying output
% vector P to fall below a fraction x of its original value P(1), given a
% time series T, or to rise above 2-x of its original value P(1) -
% whichever comes first.
% (e.g. if x = 0.9, then we find the time taken to fall outside the window
% between 0.9 and 1.1 times the original value)
% x can be a vector, in which case lifes will be a vector of the same size.
% T and P must be the same size.
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

% First set the time variable to start from 0, and preallocate the lifes
% output vector.
T = T-T(1);
lifes = zeros(size(x));

% We use a for loop to calculate the lifes separately so we can evaluate
% whether they exist before assigning them.
% Lives are calculated by interpolating the two times steps either side of
% where P hits its desired level.
for i = 1:numel(x)
    if isempty(find(abs(P - P(1))>(1-x(i))*P(1),1))
        lifes(i) = NaN;
    else
        index = min(find(abs(P - P(1))>(1-x(i))*P(1),1),numel(T));
        Ti = T(index);
        Pi = P(index);
        Ti_1 = T(index-1);
        Pi_1 = P(index-1);
        if Pi_1 > Pi
            P0 = P(1)*x(i);
        else
            P0 = P(1)+P(1)*(1-x(i));
        end
        lifes(i) = Ti_1 + (P0-Pi_1)*(Ti-Ti_1)/(Pi-Pi_1);
    end
end
end