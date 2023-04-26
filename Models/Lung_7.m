%% Computational Lung Model (ver 7)
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
fprintf("<strong># File: Lung_7.m</strong>\n");                             %lgf
fprintf("<strong># Title: Computational Lung Model (ver 7)</strong>\n");    %lgf
fprintf("<strong># Authors:</strong>\n#   - Neil A. Kumar\n#   - " + ...    %lgf
    "Helio Sulbaran\n");                                                    %lgf
fprintf("<strong># Dependencies:</strong> None\n");                         %lgf

% Create A Figure Counter
figureCounter = 1;

%% Parameter Initialization
fprintf('\n--- Parameter Initialization ---\n'); %lgf

% Initialize Time Parameters
total_time = 60*1;      % [     s     ] Total simulation time
dt = 0.01;              % [     s     ] Time step 
time = 0:dt:total_time; % [     s     ] Time vector
fprintf('Time parameters initialized\n'); %lgf

% Initialize Molecule Parameters
beta__O2 = 4.5;         % [ unit-less ] O2 oil-water partition coefficient
beta_CO2 = 4.5;         % [ unit-less ] CO2 oil-water partition coefficient
D__O2 = 1;              % [  cm^2/s   ] O2 Diffusion coefficient (pressure)
D_CO2 = 1;              % [  cm^2/s   ] CO2 Diffusion coefficient (pressure)
fprintf('Molecule parameters initialized\n'); %lgf

% Initialize Environmental Parameters
P__O2_Air = 160;        % [   mmHg    ] Room air O2 partial pressure
P_CO2_Air = 0.3;        % [   mmHg    ] Room air CO2 partial pressure
fprintf('Environmental parameters initialized\n'); %lgf

% Create Model Constants
V_Tidal = 500e-3;       % [     L     ] Alveolar tidal volume (assume 500 mL)
BR = 16;                % [breaths/min] Breath Rate
SA = 1e-6;              % [   cm^2    ] Surface Area
l = 1e-4;               % [    cm     ] Width of lung membrane
V_Cap = 1;              % [     L     ] Capilary volume (current random)
Q__Cap = V_Cap/10;      % [    L/s    ] Capilary flow rate
P__O2_Art = 0.1;        % [   mmHg    ] Arterial blood O2 partial pressure
P_CO2_Art = 35;         % [   mmHg    ] Arterial blood CO2 partial pressure
V_Dead = 150e-3; % [L] Deadspace volume (assume 150 mL)
fprintf('Model constants created\n'); %lgf

% Volume functions 
% (note: dVdt_Func should be an analytic derivative of V_Func)
V_Func = @(t) (V_Tidal/2) * sin(2*pi*BR*t/60) + 5;
dVdt_Func = @(t) (2*pi*BR*V_Tidal/120) * cos(2*pi*BR*t/60);

% Set Initial Conditions
initCond = struct;
P_O2_alv_init  = 160;   % [   mmHg    ] Partial pressure of O2 in alveoli
P_O2_cap_init  = 0.1;   % [   mmHg    ] Partial pressure of O2 in capillary
P_CO2_alv_init = 0.3;   % [   mmHg    ] Partial pressure of CO2 in alveoli
P_CO2_cap_init =  35;   % [   mmHg    ] Partial pressure of CO2 in capillary
fprintf('Initial conditions set\n'); %lgf

%% Main Computation Loop
fprintf('\n--- Main Computationn Loop ---\n'); %lgf

% Preallocate Memory for Variables
dPdt__O2_cap = zeros(size(time));
dPdt_CO2_cap = zeros(size(time));
dPdt__O2_alv = zeros(size(time));
dPdt_CO2_alv = zeros(size(time));
P__O2_alv = zeros(size(time));
P_CO2_alv = zeros(size(time));
P__O2_cap = zeros(size(time));
P_CO2_cap = zeros(size(time));
V_Alv = zeros(size(time));
fprintf('Space preallocated for variables\n'); %lgf

% Store Initial Conditions
P__O2_alv(1) = P_O2_alv_init;
P_CO2_alv(1) = P_CO2_alv_init;
P__O2_cap(1) = P_O2_cap_init;
P_CO2_cap(1) = P_CO2_cap_init;
fprintf('Initial conditions stored\n'); %lgf

fprintf('Entering main loop...');
for currentTimeStep = 1:1:length(time)
    
    % Find the current volume and flow rates (alveolus)
    V_Alv(currentTimeStep) = V_Func(time(currentTimeStep));
    dVdt = dVdt_Func(time(currentTimeStep));
    if dVdt > 0
        % Flow into the lungs
        Q__AIn = dVdt;
        Q_AOut = 0;
    elseif dVdt < 0
        % Flow out of the lungs
        Q__AIn = 0;
        Q_AOut = -dVdt;
    else
        % No flow
        Q__AIn = 0;
        Q_AOut = 0;
    end

    % Calculate partial pressure changes due to diffusion (Ficks 1st Law)
    diffRate__O2Alv = SA*D__O2*beta__O2 * ( (P__O2_cap(currentTimeStep) - P__O2_alv(currentTimeStep))/l );
    diffRate_CO2Alv = SA*D_CO2*beta_CO2 * ( (P_CO2_cap(currentTimeStep) - P_CO2_alv(currentTimeStep))/l );
    diffRate__O2Cap = SA*D__O2*beta__O2 * ( (P__O2_alv(currentTimeStep) - P__O2_cap(currentTimeStep))/l );
    diffRate_CO2Cap = SA*D_CO2*beta_CO2 * ( (P_CO2_alv(currentTimeStep) - P_CO2_cap(currentTimeStep))/l );

    % Calculate partial pressure changes due to convection
    flowRate__O2Alv = ( Q__AIn*P__O2_Air - Q_AOut*P__O2_alv(currentTimeStep) ) / V_Alv(currentTimeStep);
    flowRate_CO2Alv = ( Q__AIn*P_CO2_Air - Q_AOut*P_CO2_alv(currentTimeStep) ) / V_Alv(currentTimeStep);
    flowRate__O2Cap = ( Q__Cap*P__O2_Art - Q__Cap*P__O2_cap(currentTimeStep) ) / V_Cap;
    flowRate_CO2Cap = ( Q__Cap*P_CO2_Art - Q__Cap*P_CO2_cap(currentTimeStep) ) / V_Cap;

    % Compute total partial pressure change (diffusion + convective flow)
    dPdt__O2_alv(currentTimeStep) = diffRate__O2Alv + flowRate__O2Alv;
    dPdt_CO2_alv(currentTimeStep) = diffRate_CO2Alv + flowRate_CO2Alv;
    dPdt__O2_cap(currentTimeStep) = diffRate__O2Cap + flowRate__O2Cap;
    dPdt_CO2_cap(currentTimeStep) = diffRate_CO2Cap + flowRate_CO2Cap;
    
    % Update alveolar partial pressures
    P__O2_alv(currentTimeStep+1) = P__O2_alv(currentTimeStep) + (dPdt__O2_alv(currentTimeStep) * dt);
    P_CO2_alv(currentTimeStep+1) = P_CO2_alv(currentTimeStep) + (dPdt_CO2_alv(currentTimeStep) * dt);
    P__O2_cap(currentTimeStep+1) = P__O2_cap(currentTimeStep) + (dPdt__O2_cap(currentTimeStep) * dt);
    P_CO2_cap(currentTimeStep+1) = P_CO2_cap(currentTimeStep) + (dPdt_CO2_cap(currentTimeStep) * dt);

end
fprintf('Success! (Loop completed)\n'); %lgf

%% Result Plotting
fprintf('\n--- Result Plotting ---\n'); %lgf

% Plot the volume of the lungs
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
plot(time, V_Alv);
xlabel('Time (s)');
ylabel('Volume (L)');
title('Alveolar Volume Over Time');
fprintf('Results printed to figure %i\n',figureCounter-1);
ylim([0 6]);

% Plot the results (Partial Pressures)
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
subplot(2, 2, 1);
    plot(time, P__O2_alv(2:end));
    xlabel('Time (s)');
    ylabel('Alveolus O_2 Partial Pressure (mmHg)');
    title('Alveolus Oxygen Partial Pressure');
    fprintf('Results printed to figure %i, sublot 1\n',figureCounter-1);
subplot(2, 2, 2);
    plot(time, P_CO2_alv(2:end));
    xlabel('Time (s)');
    ylabel('Alveolus CO_2 Partial Pressure (mmHg)');
    title('Alveolus Carbon Dioxide Partial Pressure');
    fprintf('Results printed to figure %i, sublot 2\n',figureCounter-1);
subplot(2, 2, 3);
    plot(time, P__O2_cap(2:end));
    xlabel('Time (s)');
    ylabel('Capillary O_2 Partial Pressure (mmHg)');
    title('Capillary Oxygen Partial Pressure');
    fprintf('Results printed to figure %i, sublot 3\n',figureCounter-1);
subplot(2, 2, 4);
    plot(time, P_CO2_cap(2:end));
    xlabel('Time (s)');
    ylabel('Capillary CO_2 Partial Pressure (mmHg)');
    title('Capillary Carbon Dioxide Partial Pressure');
    fprintf('Results printed to figure %i, sublot 4\n',figureCounter-1);
sgtitle('Partial Pressure');

% Plot the results (Change)
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
subplot(2, 2, 1);
    plot(time, dPdt__O2_alv);
    xlabel('Time (s)');
    ylabel('Change in Alveolus O_2 (mmHg/s)');
    title('Change in Alveolus Oxygen');
    fprintf('Results printed to figure %i, sublot 1\n',figureCounter-1);
subplot(2, 2, 2);
    plot(time, dPdt_CO2_alv);
    xlabel('Time (s)');
    ylabel('Change in Alveolus CO_2 (mmHg/s)');
    title('Change in Alveolus Carbon Dioxide');
    fprintf('Results printed to figure %i, sublot 2\n',figureCounter-1);
subplot(2, 2, 3);
    plot(time, dPdt__O2_cap);
    xlabel('Time (s)');
    ylabel('Change in Capillary O_2 (mmHg/s)');
    title('Change in Capillary Oxygen');
    fprintf('Results printed to figure %i, sublot 3\n',figureCounter-1);
subplot(2, 2, 4);
    plot(time, dPdt_CO2_cap);
    xlabel('Time (s)');
    ylabel('Change in Capillary CO_2 (mmHg/s)');
    title('Change in Capillary Carbon Dioxide');
    fprintf('Results printed to figure %i, sublot 4\n',figureCounter-1);
sgtitle('Change in Partial Pressure');

%% Program End
fprintf("\n<strong>## End of Program</strong>\n"); %lgf
