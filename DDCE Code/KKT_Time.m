function [HH, A2, A1, A0, KKT] = KKT_Time(DATA)
% KKT_TIME builds the KKT system for a generic time step t > 0 
% in the data-driven problem with nonholonomic constraints 
% and no holonomic constraints (nc1 = 0).
%
% INPUT:
%   DATA - structure containing circuit/system data:
%       .L       -> inductance matrix (ne x ne)
%       .C       -> capacitance matrix (ne x ne, inverted form)
%       .R       -> resistance matrix (ne x ne)   [not explicitly used here]
%       .H       -> incidence matrix (nonholonomic constraints)
%       .G       -> (optional) matrix for resistors/damping
%       .Dt      -> time step
%       .IndxCap -> indicator vector (1 if edge is capacitor, 0 otherwise)
%       .Const   -> constitutive constants for capacitors
%
% OUTPUT:
%   HH  - constitutive block matrix (capacitor subproblem)
%   A2  - block matrix (coefficient for step n+1, "future")
%   A1  - block matrix (coefficient for step n, "current")
%   A0  - block matrix (coefficient for step n-1, "past")
%   KKT - full KKT matrix for step n+1
%
% The block matrices (A2, A1, A0) encode the time-discretized
% dynamics of inductors, capacitors, and constraints:
%
%   A2 * Y_{n+1} + A1 * Y_n + A0 * Y_{n-1} = F
%
% The KKT system has the form:
%       [ HH   A2ᵀ ] [Y] = [F]
%       [ A2    0 ] [λ]   [0]
% where HH encodes the constitutive relations of capacitors
% and A2 encodes the constraints and time-discretized dynamics.

% ---------------------------------------------------------
% Extract data
% ---------------------------------------------------------
L  = DATA.L;       % inductance matrix
H  = DATA.H;       % incidence matrix (constraints)
Dt = DATA.Dt;      % time step

ne  = size(L,1);   % number of edges (circuit elements)
nc1 = 0;           % number of holonomic constraints (set to zero)
nc2 = size(H,1);   % number of nonholonomic constraints

r  = sum(DATA.IndxCap); % number of capacitors

% Dimensions of the extended state
n1 = ne + 2*r + nc1 + nc2;
m1 = ne + r + nc1 + nc2;

% ---------------------------------------------------------
% Assemble A2 matrix (contribution from step n+1)
% ---------------------------------------------------------
A2 = zeros(m1, n1);

% Capacitor index selector
nonzero_indices = find(DATA.IndxCap ~= 0);
Vi              = diag(DATA.IndxCap);
Vi              = Vi(:, nonzero_indices);

% Fill A2
A2(1:ne, 1:ne)                    =  1/Dt * L;
A2(1:ne, (ne+r)+1:ne+2*r)         =  Dt/2 * Vi;
A2(1:ne, ne+2*r+nc1+1:end)        =  Dt * H.';
A2(ne+1:ne+r, 1:ne)               = -1/2 * Vi.';
A2(ne+1:ne+r, ne+1:ne+r)          =  eye(r);
A2(ne+r+nc1+1:end, 1:ne)          =  1/Dt * H;

% ---------------------------------------------------------
% Assemble A1 matrix (contribution from step n)
% ---------------------------------------------------------
A1 = zeros(m1, n1);

A1(1:ne, 1:ne)              = -2/Dt * L;
A1(1:ne, (ne+r)+1:ne+2*r)   =  Dt/2 * Vi;
A1(ne+1:ne+r, 1:ne)         = -1/2 * Vi.';
A1(ne+r+nc1+1:end, 1:ne)    = -1/Dt * H;

% ---------------------------------------------------------
% Assemble A0 matrix (contribution from step n-1)
% ---------------------------------------------------------
A0 = zeros(m1, n1);

A0(1:ne, 1:ne)              =  1/Dt * L;

% -------------------------------------------------------------------------
% Assemble HHi matrix (DDCE cost function)
% CC encodes the constants used in the cost function: 
% k/2|q-q*|^2 + 1/(2k)|vi-vi*|^2
% -------------------------------------------------------------------------
CC = [DATA.Const(1), DATA.Const(2)];
CC = diag(CC);

HH = zeros(n1, n1);
HH(ne+1:ne+r, ne+1:ne+r)         = CC * eye(r);     
HH(ne+r+1:ne+2*r, ne+r+1:ne+2*r) = inv(CC) * eye(r); 

% ---------------------------------------------------------
% Assemble full KKT matrix
% ---------------------------------------------------------
nH = size(HH,1);
nA = size(A2.',2);
nn = nH + nA;

KKT = zeros(nn, nn);

KKT(1:nH, 1:nH)     = HH;
KKT(1:nH, nH+1:end) = A2.';
KKT(nH+1:end, 1:nH) = A2;

end