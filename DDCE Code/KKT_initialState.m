function [HHi, Ai1, b0, KKT] = KKT_initialState(DATA)
% KKT_INITIALSTATE builds the KKT system at time t = 0 
% for the data-driven problem with nonholonomic constraints 
% and no holonomic constraints (nc1 = 0).
%
% INPUT:
%   DATA - structure containing circuit/system data:
%       .L       -> inductance matrix (ne x ne)
%       .C       -> capacitance matrix (ne x ne, inverted form)
%       .H       -> incidence matrix (nc2 x ne)
%       .Dt      -> time step
%       .IndxCap -> indicator vector (1 if edge is capacitor, 0 otherwise)
%       .Dq0     -> initial derivative of charges
%       .q0      -> initial charges
%       .Const   -> constraint constants (vector of length 2)
%
% OUTPUT:
%   HHi - block matrix related to capacitor constitutive laws
%   Ai1 - constraint matrix at t = 0
%   b0  - right-hand side vector at t = 0
%   KKT - full KKT matrix at t = 0
%
% The KKT system has the generic form:
%       [ HHi   Ai1ᵀ ] [Y] = [F]
%       [ Ai1    0  ] [λ]   [b0]
% where HHi encodes constitutive relations, Ai1 encodes constraints,
% and b0 encodes initial conditions. 

% ---------------------------------------------------------
% Extract data
% ---------------------------------------------------------
L  = DATA.L;       % inductance matrix
H  = DATA.H;       % incidence matrix (constraints)
Dt = DATA.Dt;      % time step

ne  = size(L,1);   % number of edges (elements)
nc1 = 0;           % number of holonomic constraints (set to zero here)
nc2 = size(H,1);   % number of nonholonomic constraints

r  = sum(DATA.IndxCap); % number of capacitors

% Dimensions of the extended system
n1 = ne + 2*r + nc1 + nc2;
m1 = ne + r + nc1 + nc2;

% ---------------------------------------------------------
% Assemble Ai1 matrix (constraint + dynamics)
% ---------------------------------------------------------
Ai1 = zeros(m1, n1);

% Extract capacitor indices
nonzero_indices = find(DATA.IndxCap ~= 0);
Vi              = diag(DATA.IndxCap);       % diagonal selector
Vi              = Vi(:, nonzero_indices);   % restrict to capacitor edges

% Fill Ai1 (time t+Dt contributions)
Ai1(1:ne, 1:ne)                    = 1/Dt * L;
Ai1(1:ne, (ne+r)+1:ne+2*r)         = Dt/2 * Vi;
Ai1(1:ne, ne+2*r+nc1+1:end)        = Dt * H.';
Ai1(ne+1:ne+r, 1:ne)               = -1/2 * Vi.';
Ai1(ne+1:ne+r, ne+1:ne+r)          = eye(r);
Ai1(ne+r+nc1+1:end, 1:ne)          = 1/Dt * H;

% ---------------------------------------------------------
% Assemble Ai0 matrix (time t contributions)
% ---------------------------------------------------------
Ai0 = zeros(m1, n1);

Ai0(1:ne, 1:ne)              = -1/Dt * L;
Ai0(ne+1:ne+r, 1:ne)         = -1/2 * Vi.';
Ai0(ne+r+nc1+1:end, 1:ne)    = -1/Dt * H;

% ---------------------------------------------------------
% Build right-hand side vector b0
% ---------------------------------------------------------
P0 = L * DATA.Dq0;           % initial momentum-like term

F0 = zeros(m1, 1);           % forcing vector
F0(1:ne,1) = P0;

Y0 = zeros(n1, 1);           % initial state vector
Y0(1:ne,1) = DATA.q0;        % initial charges

b0 = F0 - Ai0 * Y0;          % effective RHS

% -------------------------------------------------------------------------
% Assemble HHi matrix (DDCE cost function)
% CC encodes the constants used in the cost function: 
% k/2|q-q*|^2 + 1/(2k)|vi-vi*|^2
% -------------------------------------------------------------------------
CC = [DATA.Const(1), DATA.Const(2)];
CC = diag(CC);

HHi = zeros(n1, n1);
HHi(ne+1:ne+r, ne+1:ne+r) = CC * eye(r);      
HHi(ne+r+1:ne+2*r, ne+r+1:ne+2*r) = inv(CC) * eye(r);

% -------------------------------------------------------------------------
% Assemble full KKT system
% -------------------------------------------------------------------------
nH = size(HHi,1);
nA = size(Ai1.',2);
nn = nH + nA;

KKT = zeros(nn, nn);

KKT(1:nH, 1:nH)     = HHi;
KKT(1:nH, nH+1:end) = Ai1.';
KKT(nH+1:end, 1:nH) = Ai1;

end