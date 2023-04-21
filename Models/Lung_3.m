%% Computational Lung Model (ver 3)
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
fprintf("<strong># File: Lung_3.m</strong>\n");                             %lgf
fprintf("<strong># Title: Computational Lung Model (ver 3)</strong>\n");    %lgf
fprintf("<strong># Authors:</strong>\n#   - Neil A. Kumar\n#   - " + ...    %lgf
    "Helio Sulbaran\n");                                                    %lgf
fprintf("<strong># Dependencies:</strong> None\n");                         %lgf

% Create A Figure Counter
figCount = 1;

%% Parameter Initialization
fprintf('\n--- Parameter Initialization ---\n'); %lgf

% Initialize Time Parameters
total_time = 60*1;      % [     s     ] Total simulation time
dt = 0.01;              % [     s     ] Time step 
t = 0:dt:total_time;    % [     s     ] Time vector
fprintf('Time parameters initialized\n'); %lgf

% Create Constants
V_A = 500e-6;           % [     L     ] Alveolar tidal volume (assume 500 mL)
SA = 1e-6;              % [   cm^2    ] Surface Area
R = 8.314;              % [ J/(mol*K) ] Molar gas constant
T = 310;                % [     K     ] Body temperature
B_O2 = 4.5;             % [ unit-less ] O2 Partition coefficient
B_CO2 = 4.5;            % [ unit-less ] CO2 Partition coefficient
D_O2 = 1;               % [  cm^2/s   ] O2 Diffusion coefficient (pressure)
D_CO2 = 1;              % [  cm^2/s   ] CO2 Diffusion coefficient (pressure)
l = 1e-4;               % [    cm     ] Width of lung membrane
fprintf('Model constants created\n'); %lgf

% Set Initial Conditions
P_O2_alv_init  = 160;   % [   mmHg    ] Partial pressure of O2 in alveoli
P_O2_cap_init  = 0.1;   % [   mmHg    ] Partial pressure of O2 in capillary
P_CO2_alv_init = 0.3;   % [   mmHg    ] Partial pressure of CO2 in alveoli
P_CO2_cap_init =  35;   % [   mmHg    ] Partial pressure of CO2 in capillary
fprintf('Initial conditions set\n'); %lgf

%% Main Computation Loop
fprintf('\n--- Main Computationn Loop ---\n'); %lgf

% Preallocate Memory for Variables
dPdt_O2_cap = zeros(size(t));
dPdt_CO2_cap = zeros(size(t));
dPdt_O2_alv = zeros(size(t));
dPdt_CO2_alv = zeros(size(t));
P_O2_alv = zeros(size(t));
P_CO2_alv = zeros(size(t));
P_O2_cap = zeros(size(t));
P_CO2_cap = zeros(size(t));
P_O2_alv(1) = P_O2_alv_init;
P_CO2_alv(1) = P_CO2_alv_init;
P_O2_cap(1) = P_O2_cap_init;
P_CO2_cap(1) = P_CO2_cap_init;
fprintf('Space preallocated for variables\n'); %lgf

fprintf('Entering main loop...');
for i = 1:1:length(t)
    % Compute exchange rates (use simple linear relationships)
    dPdt_O2_cap(i) = (D_O2*B_O2 * ( (P_O2_alv(i) - P_O2_cap(i))/l )) * SA;
    dPdt_CO2_cap(i) = (D_CO2*B_CO2 * ( (P_CO2_alv(i) - P_CO2_cap(i))/l )) * SA;
    dPdt_O2_alv(i) = -dPdt_O2_cap(i);
    dPdt_CO2_alv(i) = -dPdt_CO2_cap(i);
    
    % Update alveolar partial pressures
    P_O2_alv(i+1) = P_O2_alv(i) + (dPdt_O2_alv(i) * dt);
    P_CO2_alv(i+1) = P_CO2_alv(i) + (dPdt_CO2_alv(i) * dt);
    P_O2_cap(i+1) = P_O2_cap(i) + (dPdt_O2_cap(i) * dt);
    P_CO2_cap(i+1) = P_CO2_cap(i) + (dPdt_CO2_cap(i) * dt);
    
    % Ensure alveolar partial pressures remain within physiologic limits
    % P_O2_alv(i+1) = max(0, min(P_O2_alv(i+1), 150));
    % P_CO2_alv(i+1) = max(0, min(P_CO2_cap(i+1), 80));
    % P_O2_cap(i+1) = max(0, min(P_O2_alv(i+1), 150));
    % P_CO2_cap(i+1) = max(0, min(P_CO2_cap(i+1), 80));
    %{
        This is going to have to be workshopped.  The linear model is clearly
        not working so we might want to jump straight into the more complex models
    %}
end
fprintf('Success! (Loop completed)\n'); %lgf

%% Result Plotting
fprintf('\n--- Result Plotting ---\n'); %lgf

% Plot the results (Partial Pressures)
figure(figCount); clf; figCount=figCount+1; %bppd
subplot(2, 2, 1);
    plot(t, P_O2_alv(2:end));
    xlabel('Time (s)');
    ylabel('Alveolus O_2 Partial Pressure (mmHg)');
    title('Alveolus Oxygen Partial Pressure');
    fprintf('Results printed to figure %i, sublot 1\n',figCount-1);
subplot(2, 2, 2);
    plot(t, P_CO2_alv(2:end));
    xlabel('Time (s)');
    ylabel('Alveolus CO_2 Partial Pressure (mmHg)');
    title('Alveolus Carbon Dioxide Partial Pressure');
    fprintf('Results printed to figure %i, sublot 2\n',figCount-1);
subplot(2, 2, 3);
    plot(t, P_O2_cap(2:end));
    xlabel('Time (s)');
    ylabel('Capillary O_2 Partial Pressure (mmHg)');
    title('Capillary Oxygen Partial Pressure');
    fprintf('Results printed to figure %i, sublot 3\n',figCount-1);
subplot(2, 2, 4);
    plot(t, P_CO2_cap(2:end));
    xlabel('Time (s)');
    ylabel('Capillary CO_2 Partial Pressure (mmHg)');
    title('Capillary Carbon Dioxide Partial Pressure');
    fprintf('Results printed to figure %i, sublot 4\n',figCount-1);

% Plot the results (Exchange Rate)
figure(figCount); clf; figCount=figCount+1; %bppd
subplot(2, 2, 1);
    plot(t, dPdt_O2_alv);
    xlabel('Time (s)');
    ylabel('Alveolus O_2 Exchange Rate (mmHg/s)');
    title('Alveolus Oxygen Exchange Rate');
    fprintf('Results printed to figure %i, sublot 1\n',figCount-1);
subplot(2, 2, 2);
    plot(t, dPdt_CO2_alv);
    xlabel('Time (s)');
    ylabel('Alveolus CO_2 Exchange Rate (mmHg/s)');
    title('Alveolus Carbon Dioxide Exchange Rate');
    fprintf('Results printed to figure %i, sublot 2\n',figCount-1);
subplot(2, 2, 3);
    plot(t, dPdt_O2_cap);
    xlabel('Time (s)');
    ylabel('Capillary O_2 Exchange Rate (mmHg/s)');
    title('Capillary Oxygen Exchange Rate');
    fprintf('Results printed to figure %i, sublot 3\n',figCount-1);
subplot(2, 2, 4);
    plot(t, dPdt_CO2_cap);
    xlabel('Time (s)');
    ylabel('Capillary CO_2 Exchange Rate (mmHg/s)');
    title('Capillary Carbon Dioxide Exchange Rate');
    fprintf('Results printed to figure %i, sublot 4\n',figCount-1);

%% Program End
fprintf("\n<strong>## End of Program</strong>\n"); %lgf
