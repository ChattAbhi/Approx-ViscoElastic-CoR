% This file can be used to generate all figures in the paper. Some figures
% generated figures need to be modified mannually with annotations. The
% code pertaining to different results are seperated into sections, and
% can be run independently in the matlab editor by selecting a particular
% section and then clicking Run Section in matlab editor.

error('Comment the error line to run all! Note this will overwrite previous results!')

clear all; close all; clc;

%% Generate Figures to show parameter dependence
addpath ../
ShowGamGtil
rmpath ../
clear all; close all; clc; 

%% Generate Figures to show the evolution of displacement, speed, and COR (Needs Manual Modification)
addpath ../
numericCOR(3/2,3/2,0.006,0.1,'ShowPlots',true,'figsave','Figures/matfig/NumIntResults/numeric1.fig',...
    'pdfsave','Figures/pdf/NumIntResults/numeric1.pdf')

numericCOR(3/2,3/2,0.06,4,'ShowPlots',true,'figsave','Figures/matfig/NumIntResults/numeric2.fig',...
    'pdfsave','Figures/pdf/NumIntResults/numeric2.pdf')
rmpath ../
clear all; close all; clc;

%% Generate Figures for the first-order COR comparison (Needs Manual Modification, and does not automatically save image files)
addpath ../
FirstOrderCompare(0.0001)
FirstOrderCompare(0.001)
FirstOrderCompare(0.01)
rmpath ../
clear all; close all; clc;

%% Generate Figures to show the Non-rebounding condition (Needs Manual Modification)
addpath ../

NoReboundCondition(3/2,1,0.0001,1,0.01,'case_1/Stopg','case_1/Stopg_err');
NoReboundCondition(3/2,5/4,0.0001,1,0.01,'case_2/Stopg','case_2/Stopg_err');
NoReboundCondition(3/2,3/2,0.0001,1,0.01,'case_3/Stopg','case_3/Stopg_err');
NoReboundCondition(3/2,5/2,0.0001,1,0.01,'case_4/Stopg','case_4/Stopg_err');

% NoReboundCondition(3/2,1,0.0001,0.01,'case_1/Stopg','case_1/Stopg_err');
% NoReboundCondition(3/2,1,0.0001,0.1,'case_1/Stopg2','case_1/Stopg_err2');
% NoReboundCondition(3/2,5/4,0.0001,0.01,'case_2/Stopg','case_2/Stopg_err');
% NoReboundCondition(3/2,5/4,0.0001,0.1,'case_2/Stopg2','case_2/Stopg_err2');
% NoReboundCondition(3/2,3/2,0.0001,0.01,'case_3/Stopg','case_3/Stopg_err');
% NoReboundCondition(3/2,3/2,0.0001,0.1,'case_3/Stopg2','case_3/Stopg_err2');
% NoReboundCondition(3/2,5/2,0.0001,0.01,'case_4/Stopg','case_4/Stopg_err');
% NoReboundCondition(3/2,5/2,0.0001,0.1,'case_4/Stopg2','case_4/Stopg_err2');
rmpath ../
clear all; close all; clc; 

%% Generate Figures for Schwager-Poschel comparison (Needs Manual Modification)
addpath ../
ConstCoeffFormCompare
rmpath ../
clear all; close all; clc;

%% Generate Figures to compare linearized and e^2 formulations (Needs Manual Modification)
addpath ../
NumVarCoeffFormCompare
rmpath ../
clear all; close all; clc; 

%% Generate Component Integral Approximation Comparisons for specific beta values
addpath ../
ComponentIntegrals(1)
ComponentIntegrals(2)
ComponentIntegrals(3)
ComponentIntegrals(4)
rmpath ../
clear all; close all; clc;

%% Generate Component Integral Approximation Comparisons for varying beta values
addpath ../
ComponentIntegralErrors
clear all; close all; clc;
% ComponentIntegralErrors_Alt
% clear all; close all; clc;
rmpath ../


%% Generate the Final Comparison Results Figures (Needs Manual Modification)
addpath ../

AnalVarCoeffCompare(1)
AnalVarCoeffCompare(2)
AnalVarCoeffCompare(3)
AnalVarCoeffCompare(4)
AnalVarCoeffCompare(5)
AnalVarCoeffCompare(6)
rmpath ../
clear all; close all; 
