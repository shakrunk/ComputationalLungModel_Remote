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

% Initialize Molecule Parameters (O2 and CO2)
molPara.O2.beta = 4.5;  % [ unit-less ] O2 oil-water partition coefficient
molPara.O2.D = 1;       % [  cm^2/s   ] O2 Diffusion coefficient (pressure)
molPara.CO2.beta = 4.5; % [ unit-less ] CO2 oil-water partition coefficient
molPara.CO2.D = 1;      % [  cm^2/s   ] CO2 Diffusion coefficient (pressure)
fprintf('Molecule parameters initialized\n'); %lgf

% Initialize Environmental Parameters
air.O2.P = 160;         % [   mmHg    ] Room air O2 partial pressure
air.CO2.P = 0.3;        % [   mmHg    ] Room air CO2 partial pressure
fprintf('Environmental parameters initialized\n'); %lgf

% Create Model Constants
V_Tidal = 500e-3;       % [     L     ] Alveolar tidal volume (assume 500 mL)
BR = 16;                % [breaths/min] Breath Rate
SA = 1e-6;              % [   cm^2    ] Surface Area
l = 1e-4;               % [    cm     ] Width of lung membrane
cap.V = 1;              % [     L     ] Capillary volume (current random)
cap.flow = cap.V/10;    % [    L/s    ] Capillary flow rate
art.O2.P = 0.1;         % [   mmHg    ] Arterial blood O2 partial pressure
art.CO2.P = 35;         % [   mmHg    ] Arterial blood CO2 partial pressure
dead.V = 150e-3;        % [     L     ] Deadspace volume (assume 150 mL)
fprintf('Model constants created\n'); %lgf

% Volume functions 
% (note: dVdt_Func should be an analytic derivative of V_Func)
V_Func = @(t) (V_Tidal/2) * sin(2*pi*BR*t/60) + 5;
dVdt_Func = @(t) (2*pi*BR*V_Tidal/120) * cos(2*pi*BR*t/60);

% Set Initial Conditions
alv.O2.init = 160;      % [   mmHg    ] Partial pressure of O2 in alveoli
alv.CO2.init = 0.3;     % [   mmHg    ] Partial pressure of CO2 in alveoli
cap.O2.init = 0.1;      % [   mmHg    ] Partial pressure of O2 in capillary
cap.CO2.init = 35;      % [   mmHg    ] Partial pressure of CO2 in capillary
dead.O2.init = 160;     % [   mmHg    ] Partial pressure of O2 in deadspace
dead.CO2.init = 0.3;    % [   mmHg    ] Partial pressure of CO2 in deadspace
fprintf('Initial conditions set\n'); %lgf

%% Main Computation Loop
fprintf('\n--- Main Computationn Loop ---\n'); %lgf

% Preallocate Memory for Variables    
alv.O2.P = zeros(size(time));       % Alveolus   | O2
alv.O2.dPdt = zeros(size(time));    % Alveolus   | O2
alv.CO2.P = zeros(size(time));      % Alveolus   | CO2
alv.CO2.dPdt = zeros(size(time));   % Alveolus   | CO2
alv.V = zeros(size(time));          % Alveolus   | Volume
cap.O2.P = zeros(size(time));       % Capillary  | O2
cap.O2.dPdt = zeros(size(time));    % Capillary  | O2
cap.CO2.P = zeros(size(time));      % Capillary  | CO2
cap.CO2.dPdt = zeros(size(time));   % Capillary  | CO2
dead.O2.P = zeros(size(time));      % Dead Space | O2
dead.O2.dPdt = zeros(size(time));   % Dead Space | O2
dead.CO2.P = zeros(size(time));     % Dead Space | CO2
dead.CO2.dPdt = zeros(size(time));  % Dead Space | CO2

fprintf('Space preallocated for variables\n'); %lgf

% Store Initial Conditions
alv.O2.P(1) = alv.O2.init;          % Alveolus   | O2
alv.CO2.P(1) = alv.CO2.init;        % Alveolus   | CO2
cap.O2.P(1) = cap.O2.init;          % Capillary  | O2
cap.CO2.P(1) = cap.CO2.init;        % Capillary  | CO2
dead.O2.P(1) = dead.O2.init;        % Dead Space | O2
dead.CO2.P(1) = dead.CO2.init;      % Dead Space | CO2
fprintf('Initial conditions stored\n'); %lgf

fprintf('Entering main loop...');
for currentTimeStep = 1:1:length(time)
    
    % Find the current volume and flow rates (alveolus)
    alv.V(currentTimeStep) = V_Func(time(currentTimeStep));
    dVdt = dVdt_Func(time(currentTimeStep));
    if dVdt > 0
        % Breathing in
        alv.flow.in = dVdt;
        alv.flow.out = 0;
        dead.flow.in.air = dVdt;
        dead.flow.in.alv = 0;
        dead.flow.out = -dVdt;
    elseif dVdt < 0
        % Breathing out
        alv.flow.in = 0;
        alv.flow.out = -dVdt;
        dead.flow.in.air = 0;
        dead.flow.in.alv = dVdt;
        dead.flow.out = -dVdt;
    else
        % No flow
        alv.flow.in = 0;
        alv.flow.out = 0;
        dead.flow.in.air = 0;
        dead.flow.in.alv = 0;
        dead.flow.out = 0;
    end

    % Calculate partial pressure changes due to diffusion (Ficks 1st Law)
    diffusion.O2Alv = SA*molPara.O2.D*molPara.O2.beta * ( (cap.O2.P(currentTimeStep) - alv.O2.P(currentTimeStep))/l );
    diffusion.CO2Alv = SA*molPara.CO2.D*molPara.CO2.beta * ( (cap.CO2.P(currentTimeStep) - alv.CO2.P(currentTimeStep))/l );
    diffusion.O2Cap = -diffusion.O2Alv;
    diffusion.CO2Cap = -diffusion.CO2Alv;

    % Calculate partial pressure changes due to convection
    convection.O2Dead = ( dead.flow.in.air*air.O2.P + dead.flow.in.alv*alv.O2.P(currentTimeStep) - dead.flow.out*dead.O2.P(currentTimeStep) ) / dead.V;
    convection.CO2Dead = ( dead.flow.in.air*air.CO2.P + dead.flow.in.alv*alv.CO2.P(currentTimeStep) - dead.flow.out*dead.CO2.P(currentTimeStep) ) / dead.V;
    convection.O2Alv = ( alv.flow.in*dead.O2.P(currentTimeStep) - alv.flow.out*alv.O2.P(currentTimeStep) ) / alv.V(currentTimeStep);
    convection.CO2Alv = ( alv.flow.in*dead.CO2.P(currentTimeStep)- alv.flow.out*alv.CO2.P(currentTimeStep) ) / alv.V(currentTimeStep);
    convection.O2Cap = ( cap.flow*art.O2.P- cap.flow*cap.O2.P(currentTimeStep) ) / cap.V;
    convection.CO2Cap = ( cap.flow*art.CO2.P- cap.flow*cap.CO2.P(currentTimeStep) ) / cap.V;

    % Compute total partial pressure change (diffusion + convective flow)
    dead.O2.dPdt(currentTimeStep) = convection.O2Dead;
    dead.CO2.dPdt(currentTimeStep) = convection.CO2Dead;
    alv.O2.dPdt(currentTimeStep) = diffusion.O2Alv + convection.O2Alv;
    alv.CO2.dPdt(currentTimeStep) = diffusion.CO2Alv + convection.CO2Alv;
    cap.O2.dPdt(currentTimeStep) = diffusion.O2Cap + convection.O2Cap;
    cap.CO2.dPdt(currentTimeStep) = diffusion.CO2Cap + convection.CO2Cap;
    
    % Update alveolar partial pressures
    alv.O2.P(currentTimeStep+1) = alv.O2.P(currentTimeStep) + (alv.O2.dPdt(currentTimeStep) * dt);
    alv.CO2.P(currentTimeStep+1) = alv.CO2.P(currentTimeStep) + (alv.CO2.dPdt(currentTimeStep) * dt);
    cap.O2.P(currentTimeStep+1) = cap.O2.P(currentTimeStep) + (cap.O2.dPdt(currentTimeStep) * dt);
    cap.CO2.P(currentTimeStep+1) = cap.CO2.P(currentTimeStep) + (cap.CO2.dPdt(currentTimeStep) * dt);

end; clear diffusion convection;
fprintf('Success! (Loop completed)\n'); %lgf

%% Result Plotting
fprintf('\n--- Result Plotting ---\n'); %lgf

% Plot the volume of the lungs
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
plot(time, alv.V);
xlabel('Time (s)');
ylabel('Volume (L)');
title('Alveolar Volume Over Time');
fprintf('Results printed to figure %i\n',figureCounter-1);
ylim([0 6]);

% Plot the results (Dead Space)
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
subplot(2, 1, 1);
    plot(time, dead.O2.P); hold on;
    plot(time, dead.CO2.P);
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    title('Dead Space Partial Pressure');
    fprintf('Results printed to figure %i, sublot 1\n',figureCounter-1);
subplot(2, 1, 2);
    plot(time, dead.O2.dPdt); hold on;
    plot(time, dead.CO2.dPdt);
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure Change (mmHg/s)');
    title('Change in Dead Space Partial Pressure');
    fprintf('Results printed to figure %i, sublot 2\n',figureCounter-1);
sgtitle("Dead Space");

% Plot the results (Alveolus)
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
subplot(2, 1, 1);
    plot(time, alv.O2.P(2:end)); hold on;
    plot(time, alv.CO2.P(2:end));
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    title('Alveolus Partial Pressure');
    fprintf('Results printed to figure %i, sublot 1\n',figureCounter-1);
subplot(2, 1, 2);
    plot(time, alv.O2.dPdt); hold on;
    plot(time, alv.CO2.dPdt);
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure Change (mmHg/s)');
    title('Change in Alveolus Partial Pressure');
    fprintf('Results printed to figure %i, sublot 2\n',figureCounter-1);
sgtitle("Alveolus");

% Plot the results (Capillary)
figure(figureCounter); clf; figureCounter=figureCounter+1; %bppd
subplot(2, 1, 1);
    plot(time, cap.O2.P(2:end)); hold on;
    plot(time, cap.CO2.P(2:end));
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure (mmHg)');
    title('Capillary Partial Pressure');
    fprintf('Results printed to figure %i, sublot 1\n',figureCounter-1);
subplot(2, 1, 2);
    plot(time, cap.O2.dPdt); hold on;
    plot(time, cap.CO2.dPdt);
    legend("O_2","CO_2");
    xlabel('Time (s)');
    ylabel('Pressure Change (mmHg/s)');
    title('Change in Capillary Partial Pressure');
    fprintf('Results printed to figure %i, sublot 2\n',figureCounter-1);
sgtitle("Capillary");

%% Program End
fprintf("\n<strong>## End of Program</strong>\n"); %lgf
