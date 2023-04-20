%% Computational Lung Model (ver 2)
% BIOE-5011-001 S23 Systems Physiology for BioEngineers w/ Prof. Jeffery Jacot, PhD

%% About Script
% Authors:
%  - Neil A. Kumar
%  - Helio Sulbaran
% Dependencies: None
% Formatting Notes:
%  - Look good formating denoted "lgf"
%  - Best practices for overwriting plot data denoted "bppd"

%% Setup

% Generic Reset
format default;     % Reset command window formatting
close all;          % Close all open windows / plots
clear;              % Clear the workspace of any variables
clc;                % Clear the command line

% Print this File's Metadata to Command Line
fprintf("<strong># ~/Models/   m</strong>\n");                              %lgf
fprintf("<strong># File: Lung_2.m</strong>\n");                             %lgf
fprintf("<strong># Title: Computational Lung Model (ver 2)</strong>\n");    %lgf
fprintf("<strong># Authors:</strong>\n#   - Neil A. Kumar\n#   - " + ...    %lgf
    "Helio Sulbaran\n");                                                    %lgf
fprintf("<strong># Dependencies:</strong> None\n");                         %lgf

% Create A Figure Counter
figCount = 1;

%% Parameter Initialization
fprintf('\n--- Parameter Initialization ---\n'); %lgf

% Initialize Time Parameters
total_time = 60;        % [     s     ] Total simulation time
dt = 0.01;              % [     s     ] Time step 
t = 0:dt:total_time;    % [     s     ] Time vector
fprintf('Time parameters initialized\n'); %lgf

% Create Constants
V_A = 500e-6;           % [     L     ] Alveolar tidal volume (assume 500 mL)
R = 8.314;              % [ J/(mol*K) ] Molar gas constant
T = 310;                % [     K     ] Body temperature

%%%Helio New addition%%%
%permeabilit =BD
%B B= Conc. in/ Conc. Out
%D D= Diffusion coeff.

B_O2= 40/100;
B_CO2= 45/40;
D_O2= 0.024; %Diffusion coefficient of O2 in (cm^2/s)
D_CO2= 0.000016; %Diffusion coefficient of CO2 in (um^2/s)
fprintf('Model constants created\n'); %lgf

% Set Initial Conditions
P_O2_alv = 100;         % [   mmHg    ] Partial pressure of O2 in alveoli
P_CO2_alv = 40;         % [   mmHg    ] Partial pressure of CO2 in alveoli
fprintf('Initial conditions set\n'); %lgf

%% Main Computation Loop
fprintf('\n--- Main Computationn Loop ---\n'); %lgf

% Preallocate Memory for Variables
O2_exchange_rate = zeros(size(t));
CO2_exchange_rate = zeros(size(t));
fprintf('Space preallocated for variables\n'); %lgf

fprintf('Entering main loop...');
for i = 1:1:length(t)
    % Compute exchange rates (use simple linear relationships)
    O2_exchange_rate(i) = D_O2*B_O2 * (P_O2_alv - 40);
    CO2_exchange_rate(i) = D_CO2 *B_CO2 *(45 - P_CO2_alv);
    
    % Update alveolar partial pressures
    P_O2_alv = P_O2_alv + (O2_exchange_rate(i) - CO2_exchange_rate(i)) * dt;
    P_CO2_alv = P_CO2_alv + (CO2_exchange_rate(i) - O2_exchange_rate(i)) * dt;
    
    % Ensure alveolar partial pressures remain within physiologic limits
    P_O2_alv = max(0, min(P_O2_alv, 150));
    P_CO2_alv = max(0, min(P_CO2_alv, 80));
    %{
        This is going to have to be workshopped.  The linear model is clearly
        not working so we might want to jump straight into the more complex models
    %}
end
fprintf('Success! (Loop completed)\n'); %lgf

%% Result Plotting
fprintf('\n--- Result Plotting ---\n'); %lgf

% Plot the results
figure(figCount); clf; figCount=figCount+1; %bppd
subplot(2, 1, 1);
    plot(t, O2_exchange_rate);
    xlabel('Time (s)');
    ylabel('O_2 Exchange Rate (mmHg/s)');
    title('Oxygen Exchange Rate');
    fprintf('Results printed to figure %i, sublot 1\n',figCount);
subplot(2, 1, 2);
    plot(t, CO2_exchange_rate);
    xlabel('Time (s)');
    ylabel('CO_2 Exchange Rate (mmHg/s)');
    title('Carbon Dioxide Exchange Rate');
    fprintf('Results printed to figure %i, sublot 2\n',figCount);

%% Program End
fprintf("\n<strong>## End of Program</strong>\n"); %lgf