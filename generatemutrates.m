function mutrates = generatemutrates(M,N,sigma)
% M = no_of_mutatable_parameters
% N = no_of_dimensions_per_parameter
% sigma = global mutation rate
% mutrates = mutation matrix.
% This function defines the mutation matrix, determining the transition
% rates between different mutation states, as described in the paper.
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

baserates = zeros(N);
if N == 2
    powers = 1;
else
    powers = linspace(1,N-1,N-1);
end
for n = 1:N-1
    baserates = baserates + diag(sigma^powers(n)*ones(N-n,1),n);
end
mutrates = baserates;
if M ~= 1
    for m = 1:M-1
        mutrates = kron(eye(N),mutrates)+kron(baserates,eye(N^m));
    end
end

end
