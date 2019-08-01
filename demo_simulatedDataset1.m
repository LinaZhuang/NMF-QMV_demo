%% Name: demo_simulatedDataset1
%
%  Generate the unmixing results of NMF-QMV reported in Fig. 6-7
%  of paper:
%
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and Jose M. Bioucas-Dias,
% "Regularization Parameter Selection in Minimum Volume Hyperspectral Unmixing",
% TGRS, 2019.
%
%  URL:http://www.lx.it.pt/~bioucas/publications.html 
%      or https://sites.google.com/hkbu.edu.hk/linazhuang/home
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors: Lina Zhuang (lina.zhuang@lx.it.pt)
%         &
%         Jose M. Bioucas-Dias (bioucas@lx.it.pt)
%         Oct., 2018
%%

clear all
clc
close all
addpath('NMF-QMV');
% define random states for reproducible result
rng(100);  
%


%--------------------------------------------------------------------------
% Load Library (matrix A)
%--------------------------------------------------------------------------
% 1 - USGS           (L = 224; m = 498)
% 2 - USGS - pruned  (L = 224; m = 342) (3 deg)
% 3 - USGS - pruned  (L = 224; m = 62)  (10 deg)
% 4 - USGS - pruned  (L = 224; m = 12)  (20 deg)
% 5 - USGS - pruned  (L = 224; m = 6)   (30 deg)
% 6 - iid Gaussian (0,1)
% 7 - iid Uniform [0,1]

% size of the mixing matrix [Lxm] (only applies to libraries 6 and 7)
L = 200;

% set library
library = 3;

% number of pixels
np = 4000;




i_trial = 0;

for  SNR =  20:10:50 %SNR = 30;  
    p =   7; %6:3:15   % number of endmembers
    i_trial = i_trial + 1;
    
    % maximum number of endmembers per pixel
    p_pix = min(p,4);  % see description below
    
    
    
    % abundance related parameters
    MAX_PURIRY = 0.8;         % do not theshold abundances
    OUTLIERS   = 0;           % no outliers
    PURE_PIXELS = 'no';       % no pure pixels
    SHAPE_PARAMETER = 1;      %(Dirichlet parameter) abundances uniformely
    %distribted  on each pidel
    
    %% -------------------------  abundance generation ------------------
    %
    % Assumptions:
    %
    %   1)  the dataset contains p endmenbers
    %   2)  number of endmembers per pixel is p_pix
    %   3)  SHAPE_PARAMETER is the Dirichlet parameter applied to the sets of
    %       p_pix endmembers active at each pixel
    %   4)  the groups of active pixels are sets  with p_mix components. The
    %       sets are mod([i,1+1,...,i+p_mix-1], p) for i=1,...,p
    %    5) the groups have equal probability of being active
    %
    
    % Dirchlet parameters
    pdf_pars = zeros(p, p);
    aux = [ones(1,p_pix) zeros(SHAPE_PARAMETER,p-p_pix)];
    for i=1:p
        pdf_pars(i,:) = circshift(aux',[i-1])';
    end
    % add weights
    pdf_pars = [1/p*ones(p,1) pdf_pars];
    
    
    
    %% --------------------------------------------------------------------------
    %       Start simulation
    %--------------------------------------------------------------------------
    
    switch library
        case 1  % A_1
            load USGS_1995_Library.mat
            wavelengths = datalib(:,1);
            [dummy, indexes] = sort(wavelengths);
            A = datalib(indexes,4:end);
            names = names(4:end,:);
            clear datalib;
        case 2
            load USGS_pruned_3_deg.mat
            A = B;
            clear B;
        case 3
            load USGS_pruned_10_deg.mat
            %A = B(1:6:end,:);
            A = B;
            clear B;
        case 4
            load USGS_pruned_20_deg.mat
            %A = B(1:30:end,:);
            A = B;
            clear B;
        case 5
            load USGS_pruned_30_deg.mat
            A = B;
            clear B;
        case 6
            A = randn(L,2*L);
        case 7
            A = rand(L,2*L);
        otherwise
            disp('Unknown library')
    end
    
    [L,m] = size(A);
    % L = number of bands; m = total number of materials in the library
    % normalize A
    % nA = sqrt(sum(A.^2));
    % A = A./repmat(nA,L,1);
    
    
    %%
    % -------------------------------------------------------------------
    %            Generate data
    % -------------------------------------------------------------------
    
    % mixing matrix
    index = randperm(m);
    M = A(:,index(1:p)); %true endmembers
    
    % generate the data
    n_aux  = 0;
    Y = [];
    Xaux = [];
    N = [];
    while n_aux < np
        [Ya,Xaaux,Na] = spectMixGen(M,np, ...
            'Source_pdf', 'Diri_mix', ...
            'pdf_pars',pdf_pars,...
            'max_purity',ones(1,p), ...
            'no_outliers',OUTLIERS, ...
            'pure_pixels', PURE_PIXELS, ...
            'violation_extremes',[1,1.2], ....
            'snr', SNR, ...
            'noise_shape','uniform');
        
        
        mask =   sum(Xaaux > MAX_PURIRY) ==  0;
        Y = [Y Ya(:,mask)];
        Xaux = [Xaux Xaaux(:,mask)];
        N = [N Na(:,mask)];
        n_aux = length(Y);
        
    end
    
    Y = Y(:,1:np);
    X = Xaux(:,1:np); %true abundance
    N = N(:,1:np);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    endmember_true = M;
    abundance_true = X;
    img = reshape(Y', [np 1 L]); %a image of size np rows x 1 column x L bands
    beta_candidates =  10.^(-5:5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iterm = 0;
    time_selection  = [];
    for  item = 1:3
        iterm = iterm +1;
        switch item
            case 1
                term =  'boundary';
            case 2
                term = 'center';
            case 3
                term ='totalVar';
        end
        
        
        
        %%  
        %--------------------------------------------------------------------------
        %  Unmixing based on NMF-QMV
        %--------------------------------------------------------------------------
           t_start=clock;
        
        [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, ...
            term, 'ENDMEMBER_TRUE',endmember_true, 'ABUNDANCE_TRUE', ...
            abundance_true,'DRAWFIGS','no');
        
         t_end=clock;
        time_selection = [time_selection;etime(t_end,t_start)];
        
        
        %%  
        %--------------------------------------------------------------------------
        %  Figures showing results
        %--------------------------------------------------------------------------
         
        
        if  strcmp(term,'totalVar')
            term = 'TV';
        end
        
        
        fig = figure(item);
        left_color = [0 0 1];
        right_color = [0.3 0.3 0.3];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        subplot(2,4, i_trial );
        
        yyaxis left;
        semilogx(results_save(:,1),results_save(:,2),'-o','LineWidth',2 ,'MarkerSize',10);
        hold on;
        [~,tmp] = min(results_save(:,2));
        semilogx(results_save(tmp,1),results_save(tmp,2),'o','MarkerFaceColor','r','MarkerSize',9);
        ylabel('NMSE_A','FontSize',14);
        
        
        yyaxis right;
        semilogx(results_save(:,1),results_save(:,4),'-o','LineWidth',2 ,'MarkerSize',10);
        hold on;
        [~,tmp] = min(results_save(:,4));
        semilogx(results_save(tmp,1),results_save(tmp,4),'o','MarkerFaceColor','r','MarkerSize',9);
        ylabel('$D(\mathcal{G},\widehat{\bf M})$','Interpret','latex','FontSize',14);
        
        
        
        xlabel('$\beta$','Interpret','latex','FontSize',14);
        
        set(gca,'XTick',10.^(-5:5)) ;
        set(gca,'FontSize',13.5);
        
        title({[' SNR = ', num2str(SNR), ' dB, p = ',num2str(p)],...
            ['MV = ''', term,''', ', num2str(round(time_selection(end))), 's' ] },'FontSize',14);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        left_color = [0 .6 0];
        right_color = [0.3 0.3 0.3];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        
        subplot(2,4, i_trial+4 );
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
         title({[' SNR = ', num2str(SNR), ' dB, p = ',num2str(p)],...
            ['MV = ''', term,''', ', num2str(round(time_selection(end))), 's' ] },'FontSize',14);
        
        
    end
    
    
    
end

   
        
        

    
 