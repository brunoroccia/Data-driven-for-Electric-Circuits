function [qtilde, vtilde] = phiESOperator(q, v, qdata, vdata, p)

% PHIESOPERATOR Feedback operator for DDCM (Data-Driven Computational Mechanics).
%
% This function projects a given state (q, v) onto the closest 
% data point (qtilde, vtilde) in the finite capacitor database. 
% It applies the pointwise operator phiES() for each component, 
% which finds the nearest admissible data pair in the 
% (charge-like, voltage-like) space.
%
% INPUT:
%   q      -> vector of charge-like variables at current step
%   v      -> vector of voltage-like variables at current step
%   qdata  -> database of admissible charge values 
%             (rows correspond to components, columns to data samples)
%   vdata  -> database of admissible voltage values 
%             (same structure as edata)
%   p      -> norm parameter (p = 2 usually, for Euclidean norm)
%
% OUTPUT:
%   qtilde -> projected charge-like variables (closest admissible data point)
%   vtilde -> projected voltage-like variables (corresponding voltage point)
%
% The projection step is key in DDEC:
%   - The solver generates a trial state (q, v) from equilibrium + compatibility.
%   - This function maps (q, v) to the closest constitutive data pair 
%     (qtilde, vtilde), enforcing material admissibility.
%   - The result is a feedback loop between electric and the data set.
%
% Notes:
%   - phiES() is called internally to handle the nearest-neighbor search
%     for each component.
%   - Transposition at the end ensures column-vector outputs.

Le = length(q);  % number of charge/voltage components

% Loop over all components and project onto closest data point
for i = 1:Le
    [qtilde(i), vtilde(i)] = phiES(q(i), v(i), qdata(i,:), vdata(i,:), p);
end

% Ensure outputs are column vectors
qtilde = qtilde';
vtilde = vtilde';

end