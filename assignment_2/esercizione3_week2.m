% =========================================================================
% Oleksandra Golub 856706
% =========================================================================
% Questo script crea un modello Simulink con un solo Chart Stateflow
% per simulare la dinamica delle cellule staminali ematopoietiche (stem)
% secondo la logica che è stata discussa a lezione con pdf "DCB 2025 04 Deterministic Sim"
% =========================================================================

clear all; 
close all;
clc;

% popolazioni iniziali
INITIAL_STEM = 100;
INITIAL_PROG = 50;

% tassi di eventi 
RATE_2S = 0.30;   % S -> 2S (auto-rinnovamento simmetrico)
RATE_ASYM = 0.50;   % S -> S + P (divisione asimmetrica)
RATE_2P = 0.15;   % S -> 2P (differenziazione simmetrica)
RATE_DEATH = 0.05;   % S -> morte (apoptosi)

% tassi di differenziazione progenitori
RATE_MACRO = 0.40;   % P -> Macrophage
RATE_GRAN  = 0.60;   % P -> Granulocyte

% passo temporale e durata simulazione
DT        = 0.1;     % passo di integrazione
SIM_TIME  = 100;     % tempo totale di simulazione

% creazione del modello simulink + stateflow
modelName = 'StemCellStateflow_min';
modelFile = [modelName '.slx'];

% se esiste già, apri
if exist(modelFile, 'file')
    % se il modello esiste già, aprilo
    fprintf('Il modello "%s" esiste già. Lo apro...\n', modelFile);
    open_system(modelName);
else
    % se non esiste, crealo da zero
    fprintf('Il modello "%s" non esiste. Lo creo ora...\n', modelFile);
    new_system(modelName);
    % aggiungi un Chart Stateflow
    add_block('sflib/Chart', [modelName '/Chart'], ...
               'Position', [100 100 380 320]);
end 

% imposta stop time e solver
set_param(modelName, 'StopTime', num2str(SIM_TIME));
set_param(modelName, 'Solver', 'ode4');          % Runge-Kutta 4
set_param(modelName, 'FixedStep', num2str(DT));        % passo base

%% salvaggio dei paramentri su workspace
assignin('base','INITIAL_STEM',INITIAL_STEM);
assignin('base','INITIAL_PROG',INITIAL_PROG);
assignin('base','RATE_2S',RATE_2S);
assignin('base','RATE_ASYM',RATE_ASYM);
assignin('base','RATE_2P',RATE_2P);
assignin('base','RATE_DEATH',RATE_DEATH);
assignin('base','RATE_MACRO',RATE_MACRO);
assignin('base','RATE_GRAN',RATE_GRAN);
assignin('base','DT',DT);
