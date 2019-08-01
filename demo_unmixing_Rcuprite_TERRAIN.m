%% Name: demo_unmixing_Rcuprite_TERRAIN
%
%-----------------------------------------------------------------------
%
%  This demo illustrates the NMF_QMV hyperspectral unmixing  algorithm
%  operating in Rcuprite image and TERRAIN image.
%
%  Results are given in the paper:
%
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and Jose M. Bioucas-Dias,
% "Regularization Parameter Selection in Minimum Volume Hyperspectral Unmixing",
% TGRS, 2019.
%
%  URL:http://www.lx.it.pt/~bioucas/publications.html
%      or https://sites.google.com/hkbu.edu.hk/linazhuang/home
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Lina Zhuang (lina.zhuang@lx.it.pt)
%          Oct., 2018
%%


%% Input
clear;clc;
close all;
rng('default');
rng(1);
addpath('NMF-QMV');
%%
%--------------------------------------------------------------------------
%  Load Rcuprite data
%--------------------------------------------------------------------------
%  % Rcuprite, of size 50 lines by 90 columns by 188 spectral bands,
% % is a subset of the well-known AVIRIS Cuprite data set.

 
load Rcuprite50x90x188; %a sub-image of Cuprite image
img = Rcuprite50x90x188;
img = img./10000;
p=3;

load USGS_lib_DesertVanish_Montmorillonite_Alunite_188b; %spectral signatures from USGS library
endmember_true = UGSG_lib_DesertVanish_Montmorillonite_Alunite_188b;
clear UGSG_lib_DesertVanish_Montmorillonite_Alunite_188b;


%%
%--------------------------------------------------------------------------
%   Load TERRAIN data
%--------------------------------------------------------------------------


%% TERRAIN
% % we simulate a semi-real image base on the publicly available TERRAIN hyperspectral
% % image acquired by the HYDICE sensor (see more details in 'GenerateHSIFromTerrain.m').

%   load SimulatedTerrainDataset;
% img = img_syn_noisy; %clear img_syn_noisy;
% abundance_true = S_low_pass;
% endmember_true = M;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_candidates =   10.^(-5:5);
results_per_term=[];
iterm = 0;
for  item = 1 : 3
    iterm = iterm +1;
    switch item
        case 1
            term =  'boundary';
        case 2
            term = 'center';
        case 3
            term ='totalVar';
    end
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------------------------------------------------------------------
    %  Unmixing based on NMF-QMV
    %--------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if ~exist('endmember_true','var')
        [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, term);
    else
        if  ~exist('abundance_true','var')
            [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, ...
                term, 'ENDMEMBER_TRUE',endmember_true);
        else
            [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, ...
                term, 'ENDMEMBER_TRUE',endmember_true, 'ABUNDANCE_TRUE', ...
                abundance_true);
            %  [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, ...
            %                 term, 'ENDMEMBER_TRUE',endmember_true, 'ABUNDANCE_TRUE', ...
            %                 abundance_true);
        end
    end
    %
    
    results_per_term{item} = results_save;
    results_A_output{item} = A_output;
    results_S_output{item} = S_output;
    results_beta_best{item} = beta_best;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------------------------------------------------------------------------
    %  Figures showing results
    %--------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if  strcmp(term,'totalVar')
        term = 'TV';
    end
    
    
    fig = figure(9);
    left_color = [0 0 1];
    right_color = [0.3 0.3 0.3];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    subplot(2,3, iterm );
    
    if exist('endmember_true','var') %if no endmember_true, results_save(:,2) = 0.
        yyaxis left;
        semilogx(results_save(:,1),results_save(:,2),'-o','LineWidth',2 ,'MarkerSize',10);
        hold on;
        [~,tmp] = min(results_save(:,2));
        semilogx(results_save(tmp,1),results_save(tmp,2),'o','MarkerFaceColor','r','MarkerSize',9);
        ylabel('NMSE_A','FontSize',14);
    end
    
    yyaxis right;
    semilogx(results_save(:,1),results_save(:,4),'-o','LineWidth',2 ,'MarkerSize',10);
    hold on;
    [~,tmp] = min(results_save(:,4));
    semilogx(results_save(tmp,1),results_save(tmp,4),'o','MarkerFaceColor','r','MarkerSize',9);
    ylabel('$D(\mathcal{G},\widehat{\bf M})$','Interpret','latex','FontSize',14);
    
    
    
    xlabel('$\beta$','Interpret','latex','FontSize',14);
    
    
    
    
    
    
    
    title({['MV = ''', term,'''' ] },'FontSize',14);
    set(gca,'XTick',10.^(-5:5)) ;
    set(gca,'FontSize',13.5);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if  exist('abundance_true','var') %if no abundance_true, results_save(:,3) = 0.
        left_color = [0 .6 0];
        right_color = [0.3 0.3 0.3];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        
        subplot(2,3, iterm+3 );
        yyaxis left;
        semilogx(results_save(:,1),results_save(:,3),'-o','LineWidth',2 ,'MarkerSize',10);
        hold on;
        [~,tmp] = min(results_save(:,3));
        semilogx(results_save(tmp,1),results_save(tmp,3),'o','MarkerFaceColor','r','MarkerSize',9);
        ylabel('NMSE_S','FontSize',14);
        
        
        yyaxis right;
        semilogx(results_save(:,1),results_save(:,4),'-o','LineWidth',2 ,'MarkerSize',10);
        hold on;
        [~,tmp] = min(results_save(:,4));
        semilogx(results_save(tmp,1),results_save(tmp,4),'o','MarkerFaceColor','r','MarkerSize',9);
        ylabel('$D(\mathcal{G},\widehat{\bf M})$','Interpret','latex','FontSize',14);
        
        
        
        xlabel('$\beta$','Interpret','latex','FontSize',14);
        set(gca,'XTick',10.^(-5:5)) ;
        set(gca,'FontSize',13.5);
        title({['MV = ''', term,'''' ] },'FontSize',14);
    end
    
end
