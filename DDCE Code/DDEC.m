%% ========================================================================
%  MAIN SCRIPT: DDCM Solver for Electric Circuits
%  ------------------------------------------------------------------------
%  This script sets up, assembles, and solves an electrical circuit using
%  Data-Driven Computational Mechanics (DDCM) with the Alternating 
%  Direction Method (ADM). 
%
%  Steps:
%   1. Define inductances and capacitances
%   2. Specify circuit topology (vertices + edges)
%   3. Set initial conditions and integrator parameters
%   4. Generate synthetic data for capacitors
%   5. Assemble system matrices (via circuit_assembly.m)
%   6. Call DDCM solver (DDCM_ADM_Electric.m)
%   7. Plot results: phase space, time evolution, convergence, cost
%
%  HOW TO RUN:
%   - Place this script in the same folder as:
%       * circuit_assembly.m
%       * DDCM_ADM_Electric.m
%       * phiESOperator.m (and phiES.m if used internally)
%   - Run the script (F5 or "Run" button in MATLAB).
%   - Results will be plotted automatically.
%
%  Author:      Dr. Ing. Bruno A. Roccia
%  Year:        2025
%  Institution: Bergen Offshore Wind Centre (BOW), University of Bergen
%  Contact:     bruno.roccia@uib.no
%  ========================================================================
close all; clear all; clc;

%% ========================================================================
% 1. Inductances and Capacitances
% -------------------------------------------------------------------------
% Define inductance (L) and capacitance (C) values for the circuit elements.
% Ls: array of inductances for each inductor.
% Cs: array of capacitances for each capacitor.
% ========================================================================
Ls = [1 1];   % Inductor values
Cs = [1 10];  % Capacitor values


%% ========================================================================
% 2. Circuit Topology (Graph Representation)
% -------------------------------------------------------------------------
% Vertices: nodes of the circuit
% Edges: elements (L or C) connecting two vertices
%
% adjacency list convention:
%   positive sign -> outgoing edge
%   negative sign -> incoming edge
% connectivity: defines which two vertices are connected by an edge
% ========================================================================

% Define vertices (nodes) and their adjacency lists
vertices(1).adjacency = [+1 -4];
vertices(2).adjacency = [+4 -2 -3];
vertices(3).adjacency = [+3 +2 -1];

% Define edges (circuit elements)
edges(1).connectivity = [1 3];  % between vertex 1 and 3
edges(1).type = 'L';            % inductor
edges(1).index = 1;             % uses Ls(1)

edges(2).connectivity = [3 2];  
edges(2).type = 'L';
edges(2).index = 2;             % uses Ls(2)

edges(3).connectivity = [3 2];  
edges(3).type = 'C';            % capacitor
edges(3).index = 1;             % uses Cs(1)

edges(4).connectivity = [2 1];
edges(4).type = 'C';
edges(4).index = 2;             % uses Cs(2)


%% ========================================================================
% 3. Integrator Setup & Initial Conditions
% -------------------------------------------------------------------------
% q0: initial charges on each branch
% qdot0: initial rates of change of charge (currents)
% deltat: time step
% TOL, TOLc: tolerances for solver
% nsteps: total number of steps in simulation
% ========================================================================
q0     = [1 1 0 1]';   % Initial charges
qdot0  = [0 0 0 0]';   % Initial charge rates (currents)

deltat = 0.01;         % Time step
TOL    = 1e-10;        % Solver tolerance
TOLc   = 1e-6;         % Cost function tolerance
nsteps = 80001;        % Number of integration steps


%% ========================================================================
% 4. Synthetic Data Generation for Capacitors
% -------------------------------------------------------------------------
% Generates admissible charge-voltage pairs (q,v) for capacitors.
% Data includes nonlinear terms to mimic realistic capacitors.
% This database is used by the feedback operator to enforce 
% data-driven constraints.
% ========================================================================
nd = 101;  % number of database points per capacitor

% Range of admissible charge values
emin_1= -1.5; emax_1 = +1.5;
emin_2= -1.5; emax_2 = +1.5;

% Linear spacing of charge values
qdata1 = linspace(emin_1, emax_1, nd); 
qdata2 = linspace(emin_2, emax_2, nd); 

% Capacitor voltage-charge relations (linear + nonlinear terms)
vdata1 = 1/Cs(1)*qdata1 + 0.00*1/Cs(1)*qdata1.^3;
vdata2 = 1/Cs(2)*qdata2 + 3.25*1/Cs(2)*qdata1.^3;

% Assemble into data matrices
qdata  = [qdata1; qdata2];
vdata  = [vdata1; vdata2];

% Compute effective tangent constants (averaged slopes)
Indx1 = find(qdata1 ~= 0);  Indx2 = find(qdata2 ~= 0);
Auxq1 = qdata1(Indx1);      Auxv1 = vdata1(Indx1);
Auxq2 = qdata2(Indx2);      Auxv2 = vdata2(Indx2);

K1 = 1/length(Auxq1) * sum(Auxv1./Auxq1);
K2 = 1/length(Auxq2) * sum(Auxv2./Auxq2);


%% ========================================================================
% 5. Assembly of the System
% -------------------------------------------------------------------------
% Uses circuit_assembly.m to build:
%   L   - Inductance matrix
%   C   - Capacitance matrix
%   MR  - Constraint matrix
%   IndxCap - Indices of capacitors
% ========================================================================
[L, C, MR, IndxCap] = circuit_assembly(edges, vertices, Ls, Cs);


%% ========================================================================
% 6. Store Parameters in DATA Structure
% -------------------------------------------------------------------------
% Pack all parameters into a structure to pass into the solver.
% ========================================================================
DATA.L   = L;
DATA.C   = C;
DATA.H   = MR;
DATA.Dt  = deltat;
DATA.q0  = q0;
DATA.Dq0 = qdot0;
DATA.IndxCap = IndxCap;
DATA.TOL  = TOL;
DATA.TOLc = TOLc;
DATA.Nsteps = nsteps;
DATA.Iter = 50;

DATA.Const(1) = K1;
DATA.Const(2) = K2;


%% ========================================================================
% 7. Solve with DDCM + ADM
% -------------------------------------------------------------------------
% Call solver: DDCM_ADM_Electric.m
% Inputs: DATA (system parameters), qdata, vdata
% Output: RESULTS (simulation results)
% ========================================================================
RESULTS = DDCM_ADM_Electric(DATA, qdata, vdata);


%% ========================================================================
% 8. Postprocessing & Plots
% -------------------------------------------------------------------------
% Generate plots:
%   1. Phase-space (q-v) comparisons
%   2. Time evolution of charges
%   3. Residuum and Error cost evolution
%   4. Global cost function
% ========================================================================

% Plot q-v data vs projected solutions
for i = 1:length(Cs)
    figure(i)
    plot (qdata(i,:), vdata(i,:), 'o','MarkerSize',6, 'MarkerFaceColor','b')
    hold on
    plot (RESULTS.qtilde(i,:), RESULTS.vtilde(i,:), 'o','MarkerSize',14, 'MarkeredgeColor','k')
    plot (RESULTS.qmid(i,:), RESULTS.Vmid(i,:), 'x','MarkerSize',4, 'MarkeredgeColor','r')
    grid on
    xlabel ('Electric charge [C]')
    ylabel ('Capacitor voltage [V]')
    legend ('Synthetic data','Discrete pair y_n','Pair (q_{n+1/2},v_{n+1/2})')
end

% Plot time evolution of charges
figure (i+1)
subplot(2,2,1)
plot (RESULTS.t, RESULTS.q(1,:), 'color','k','LineWidth',2)
box on; grid on
xlabel('Time'); ylabel('Electric charge q_1')

subplot(2,2,2)
plot (RESULTS.t, RESULTS.q(2,:), 'color','k','LineWidth',2)
box on; grid on
xlabel('Time'); ylabel('Electric charge q_2')

subplot(2,2,3)
plot (RESULTS.t, RESULTS.q(3,:), 'color','k','LineWidth',2)
box on; grid on
xlabel('Time'); ylabel('Electric charge q_3')

subplot(2,2,4)
plot (RESULTS.t, RESULTS.q(4,:), 'color','k','LineWidth',2)
box on; grid on
xlabel('Time'); ylabel('Electric charge q_4')

% Plot convergence measures
figure (i+2)
subplot(1,2,1)
plot (RESULTS.Residuum, 'color','b','LineWidth',2)
box on; grid on
xlabel('Time step'); ylabel('Residuum')

subplot(1,2,2)
plot (RESULTS.ErrorCost, 'color','b','LineWidth',2)
box on; grid on
xlabel('Time step'); ylabel('Cost Function Error')

% Plot global cost function
figure (i+3)
plot (RESULTS.Costq + RESULTS.CostV, 'color','b','LineWidth',2)
box on; grid on
xlabel('Time step'); ylabel('Cost function')