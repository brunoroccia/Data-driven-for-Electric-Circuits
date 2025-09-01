function [L, C, HR, IndxCap] = circuit_assembly(edges, vertices, Ls, Cs)
% CIRCUIT_ASSEMBLY assembles the constitutive and constraint matrices 
% for an electrical circuit represented as a graph.
%
% INPUTS:
%   edges    - structure array defining the circuit edges (elements). 
%              Each edge has fields:
%                 .type  -> 'L' for inductor, 'C' for capacitor
%                 .index -> index pointing to parameter value in Ls or Cs
%   vertices - structure array defining the circuit vertices (nodes).
%              Each vertex has a field:
%                 .adjacency -> list of incident edges with orientation
%                               (positive if leaving, negative if entering).
%   Ls       - vector of inductance values (indexed by edges(i).index)
%   Cs       - vector of capacitance values (indexed by edges(i).index)
%
% OUTPUTS:
%   L        - diagonal inductance matrix (ne x ne)
%   C        - diagonal inverse-capacitance matrix (ne x ne)
%   HR       - reduced incidence matrix (removes redundant constraints)
%   IndxCap  - indicator vector (1 if edge is a capacitor, 0 otherwise)

% Number of vertices (nodes) and edges (elements)
nv = length(vertices); 
ne = length(edges);

% Initialize matrices
L = zeros(ne);           % inductance matrix (diagonal)
C = zeros(ne);           % inverse capacitance matrix (diagonal)
IndxCap = zeros(1, ne);  % capacitor indicator vector

% ---------------------------------------------------------
% Loop over edges to assemble constitutive matrices (L, C)
% ---------------------------------------------------------
for i = 1:ne
    IndxCap(i) = 0; % default (not a capacitor)

    switch edges(i).type
        case 'L' % Inductor
            % Place inductance value in diagonal of L
            L(i,i) = Ls(edges(i).index);

        case 'C' % Capacitor
            % Place inverse capacitance in diagonal of C
            C(i,i) = 1 / Cs(edges(i).index);
            % Mark edge as capacitor
            IndxCap(i) = 1;
    end
end

% ---------------------------------------------------------
% Construct the incidence matrix H (nv x ne)
% Each row corresponds to a vertex (node)
% Each column corresponds to an edge (element)
%   - H(i,j) = +1 if edge j leaves node i
%   - H(i,j) = -1 if edge j enters node i
%   - H(i,j) = 0 otherwise
% ---------------------------------------------------------
H = zeros(nv, ne);
for i = 1:nv
    na = length(vertices(i).adjacency); % number of incident edges
    for j = 1:na
        k = abs(vertices(i).adjacency(j)); % edge index
        H(i, k) = sign(vertices(i).adjacency(j)); % orientation
    end
end

% ---------------------------------------------------------
% Reduced incidence matrix HR
% The standard incidence matrix H may contain redundant
% equations (due to Kirchhoff’s laws). We compute its
% reduced form using RREF to retain only independent
% constraints.
% ---------------------------------------------------------
HT = H.';                     % transpose incidence matrix
[~, RBH] = rref(HT);          % compute pivot columns via RREF
HR = HT(:, RBH).';            % reduced incidence matrix

end