function RESULTS = DDCM_ADM_Electric (DATA, qdata, vdata)
% DDCM_ADM_Electric
% -------------------------------------------------------------------------
% Data-Driven Computational Mechanics (DDCM) solver for electric circuits
% using the Alternating Direction Method (ADM).
%
% This function solves the time evolution of an electric network subject to
% non-holonomic constraints, where material behavior is not described by
% constitutive laws but by a finite material database (qdata, vdata).
% The solver alternates between:
%   (i)  enforcing equilibrium/compatibility (via KKT system),
%   (ii) projecting trial states (q, v) onto the closest data points
%        in the material database (via phiESOperator).
%
% INPUT:
%   DATA   -> structure containing system matrices and parameters:
%               L, C    : element matrices (inductors, capacitors)
%               H       : incidence/constraint matrices
%               Dt      : time step
%               Nsteps  : total number of steps
%               Iter    : max ADM iterations per step
%               TOL     : tolerance for residuum (equilibrium error)
%               TOLc    : tolerance for cost function convergence
%               q0      : initial charges (or generalized coords)
%               Const   : scaling constants for auxiliary variables
%               IndxCap : indicator vector for capacitor positions
%   qdata  -> admissible database of charge-like variables
%   vdata  -> admissible database of voltage-like variables
%
% OUTPUT (RESULTS struct):
%   q        : generalized charges (trajectories)
%   qdot     : time derivatives of charges
%   qmid     : midpoint charge variables
%   Vmid     : midpoint voltage variables
%   mu       : Lagrange multipliers (constraint forces)
%   qtilde   : projected charges from data
%   vtilde   : projected voltages from data
%   t, tm    : time grid and midpoints
%   Iter     : number of ADM iterations per time step
%   Costq    : projection cost in charge space
%   CostV    : projection cost in voltage space
%   ErrorCost: cost error evolution across iterations
%   Residuum : equilibrium residual at each time step
%
% METHOD:
%   1. Initial state solved using KKT system (steady state).
%   2. Time integration via variational scheme (KKT_Time).
%   3. At each step, the ADM loop iterates until both:
%        - equilibrium residual < TOL
%        - projection error < TOLc
%   4. Projection onto dataset done by phiESOperator (feedback).
%
% NOTE:
%   - The scheme is variational, structure-preserving, and
%     converges by alternating equilibrium + data admissibility.
%   - qdata and vdata encode the material response directly, without
%     fitting constitutive models.

% ----------------- INITIALIZATION -----------------
[HHi, Ai1, b0, KKT0] = KKT_initialState (DATA);

k = 1;
q0Index = find(DATA.IndxCap~=0);
qtilde0 = DATA.q0(q0Index); % initial projected charges
AUX     = DATA.C*DATA.q0;
vtilde0 = AUX(q0Index);       % initial projected voltages

Utilde  = zeros (size(HHi,1),1);
ne      = size(DATA.L,1);   % number of elements
r       = sum(DATA.IndxCap);% number of capacitors
nc1     = 0;                % optional constraints (unused)
nc2     = size(DATA.H,1);   % number of constraints

Residuum(1) = 1;
Cost0 = 0;
ErrorCost(1) = 1;

CC = diag([DATA.Const(1), DATA.Const(2)]);

% ----------------- INITIAL STEP -----------------
while ((Residuum(1) >= DATA.TOL)  || (ErrorCost(1) > DATA.TOLc)) && k <= DATA.Iter

    % Assemble RHS with projected states
    Utilde(ne+1:ne+r,1)     = - CC * qtilde0;
    Utilde(ne+r+1:ne+2*r,1) = - inv(CC) * vtilde0;
    RHS = [-Utilde; b0];

    % Solve KKT system
    X   = KKT0 \ RHS;

    % Extract midpoint trial states
    Vm  = X(ne+r+1:ne+2*r,1);
    qm  = X(ne+1:ne+r,1);

    % Project onto data
    [qtilde0, vtilde0] = phiESOperator(qm, Vm, qdata, vdata, 1);

    % Update cost and residuals
    Costq(1) = 0.5 * (qm - qtilde0).' * CC * (qm - qtilde0);
    Costv(1) = 0.5 * (Vm - vtilde0).' * inv(CC) * (Vm - vtilde0);
    Cost(1)  = Costq(1) + Costv(1);

    Residuum(1)  = norm (Ai1*X(1:ne+2*r+nc1+nc2,1) - b0, 1);
    ErrorCost(1) = norm (Cost(1)-Cost0, inf);
    Cost0 = Cost(1);

    k = k + 1;
end

% Store results for initial step
t    = [0 DATA.Dt];
tm   = DATA.Dt/2;
q      = [DATA.q0      X(1:ne,1)];
mu     = X(ne+2*r+nc1+1:ne+2*r+nc1+nc2,1);
qdot   = (q(:,2)-q(:,1))/DATA.Dt;
Vmid   = X(ne+r+1:ne+2*r,1);
qmid   = X(ne+1:ne+r,1);
vtilde = vtilde0;
qtilde = qtilde0;
Iter(1) = k;

% ----------------- TIME INTEGRATION LOOP -----------------
[HH, A2, A1, A0, KKT] = KKT_Time (DATA);
n1 = ne + 2*r + nc1 + nc2;
m1 = ne + r + nc1 + nc2;

for i = 2:DATA.Nsteps-1

    Residuum(i) = 1;
    k = 1;
    qtilde0 = qtilde(:,i-1);
    vtilde0 = vtilde(:,i-1);
    Cost0 = 0;
    ErrorCost(i) = 1;

    % ADM iterations per time step
    while ((Residuum(i) >= DATA.TOL)  || (ErrorCost(i) > DATA.TOLc)) && k <= DATA.Iter

        % Construct variational scheme (two-step)
        Y1 = zeros(n1,1);
        Y1(1:ne,1) = q(1:ne,i);
        Y1(ne+r+1:ne+2*r,1) = Vmid(:,i-1);

        Y2 = zeros(n1,1);
        Y2(1:ne,1) = q(1:ne,i-1);

        % Assemble RHS
        Utilde(ne+1:ne+r,1)     = - CC * qtilde0;
        Utilde(ne+r+1:ne+2*r,1) = - inv(CC) * vtilde0;
        b = -A1*Y1 - A0*Y2;
        RHS = [-Utilde; b];

        % Solve KKT system
        X    = KKT\RHS;
        Vm  = X(ne+r+1:ne+2*r,1);
        qm  = X(ne+1:ne+r,1);

        % Projection step
        [qtilde0, vtilde0] = phiESOperator(qm, Vm, qdata, vdata, 1);

        % Update cost and residual
        Costq(i) = 0.5 * (qm - qtilde0).' * CC * (qm - qtilde0);
        Costv(i) = 0.5 * (Vm - vtilde0).' * inv(CC) * (Vm - vtilde0);
        Cost(i)  = Costq(i) + Costv(i);

        Residuum(i)  = norm (A2*X(1:ne+2*r+nc1+nc2,1) - b, inf);
        ErrorCost(i) = norm(Cost(i)-Cost0,inf);
        Cost0 = Cost(i);

        k = k + 1;
    end

    % Update trajectories
    q      = [q        X(1:ne,1)];
    mu     = [mu       X(ne+2*r+nc1+1:ne+2*r+nc1+nc2,1)];
    qdotp  = (q(:,i+1)-q(:,i))/DATA.Dt;
    qdot   = [qdot     qdotp];
    Vmid   = [Vmid     X(ne+r+1:ne+2*r,1)];
    qmid   = [qmid     X(ne+1:ne+r,1)];
    vtilde = [vtilde   vtilde0];
    qtilde = [qtilde   qtilde0];

    % Update time grid
    t   = [t   i*DATA.Dt];
    tm  = [tm  (i-1/2)*DATA.Dt];
    Iter(i) = k;

    fprintf ('Time Step %i - Iterations %i\n', i, k)
end

% ----------------- OUTPUT -----------------
RESULTS.q = q;
RESULTS.t = t;
RESULTS.tm = tm;
RESULTS.mu = mu;
RESULTS.qmid = qmid;
RESULTS.Vmid = Vmid;
RESULTS.qtilde = qtilde;
RESULTS.vtilde = vtilde;
RESULTS.qdot = qdot;
RESULTS.Iter = Iter;
RESULTS.Costq = Costq;
RESULTS.CostV = Costv;
RESULTS.ErrorCost = ErrorCost;
RESULTS.Residuum = Residuum;

end
