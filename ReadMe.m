%%  The code and data herein distributed reproduces the results published in
%  the paper 
%
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and Jose M. Bioucas-Dias,
% "Regularization Parameter Selection in Minimum Volume Hyperspectral Unmixing",
% TGRS, 2019.
%  
%  URL:http://www.lx.it.pt/~bioucas/publications.html 
%      or https://sites.google.com/hkbu.edu.hk/linazhuang/home
%
%%  Notes:
%
%  
% 1) demo_simulatedDataset1.m <-- This demo illustrates the NMF_QMV hyperspectral 
%    unmixing algorithm operating in simulated Dataset1 (SCENARIO:  non pure 
%    pixels, various number of endmembers and noise level)
% 
% 2) demo_unmixing_Rcuprite_TERRAIN.m  <-- This demo illustrates the NMF_QMV 
%    hyperspectral unmixing  algorithm operating in Rcuprite image and TERRAIN image.
% 
% 3) GenerateHSIFromTerrain.m <-- This script illustrates the procudure to generates 
%    a simulated image from TERRAIN data.
% 
% 4) NMF_QMV.m <-- hyperspectral unmixing method based on Nonnegative Matrix 
%    Factorization via Quadractic Minimum Volume (NMF_QMV)
% 
% 5) nmf_qmv_subspace.m <-- NMF-QMV working in a dimension-reduced space.
%      
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Lina Zhuang and Jose M. Bioucas Dias, Oct. 2018
%

