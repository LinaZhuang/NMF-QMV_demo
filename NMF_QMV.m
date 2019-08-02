function [beta_best, A_output, S_output, results_save] = NMF_QMV(img, p, beta_candidates, term, varargin)
%
%
% --------------- Description ---------------------------------------------
%
%  Nonnegative Matrix Factorization via Quadractic MinimumVoule (NMF_QMV):
%
%  Let  -> Y  [bands,np]  matrix containing the observed spectral vectors in its
%                    columns
%
%       -> A  [bands,p] mixing matrix with p spectral signatures (usually bands >> p)
%
%       -> S  [p,np] abundance matrix
%
%       -> N  {bands,np] additive Gaussian Noise
%
%
%  Linear mixing observation model
%
%     Y = AS + N    with   S>= 0,   and   sum(S) = ones(1,np)
%
% -----------------------------------------------------------------------
%  Optimization problem:
%
%   min  (1/2)  ||AS-Y||^2_F + beta ||AB-O||^2_F
%    S >= 0
%    A \in (p-1)affine set defined by the data
%
%   OPTIONAL CONSTRAINTS:
%
%       Sum-To-One sum(S) = ones(1,np)
%
%   where  ||AB-O||^2_F is a term with minimum volume flavor.  The exact
%   meaning of this term depends in the couple(B,O):
%
%  'MIN_VOLUME' - {'boundary', 'center','totalVar'}  define the minimum volume
%                  term ||AB-O||^2_F
%
%                 'boundary'        B = I, O =  extremes given by VCA (or other pure pixel algorithm)
%                 'center'          B = I,  O = center of mass
%                 'totalVar'        B = eye(p) - ones(p)./p, O = 0
%
% ---------------------------- -------------------------------------------
%
% See more details in papers:
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and Jose M. Bioucas-Dias,
% "Regularization Parameter Selection in Minimum Volume Hyperspectral Unmixing",
% TGRS, 2019.
%
%  URL:http://www.lx.it.pt/~bioucas/publications.html 
%      or https://sites.google.com/hkbu.edu.hk/linazhuang/home
% -----------------------------------------------------------------------
%
% Input:
% img                 hyperspectral data set with (lines x samples x bands).
% p                   number of endmember
% beta_candidates     A set of regularization parameter candidates
% term              = {'boundary', 'center', 'totalVar' }  define the minimum volume
%                     term ||AB-O||^2_F
%                    'boundary'        B = I, O =  extremes given by VCA (or other pure pixel algorithm)
%                    'center'          B = I,  O = center of mass
%                    'totalVar'        B = eye(p) - ones(p)./p, O = 0
%
% Optional inputs:
% endmember_true      true endmember matrix of size bands x p
% abundance_true      true abundance matrix of size p x np, where np = lines x samples.
% drawFigs           = {'yes', 'no'}
%                     'no' - work silently
%                     'yes' -  plot figures
%                     Default 'yes'
%
% ---------------------------- -------------------------------------------
%
% Output:
% beta_best             selected regularization parameter
% A_output              endmember matrix of size bands x p
% S_output              abundance matrix of size p x np, where np = lines x samples.
% results_save          a matrix containing parameter candidates (in first column),
%                       NMSE_A (in second column), NMSE_S (in third column),
%                       proposed criterion (in fourth column)
%
%
%% -------------------------------------------------------------------------
%
% Copyright (Oct. 2018):
%             Lina Zhuang (lina.zhuang@lx.it.pt)
%             &
%             José Bioucas-Dias (bioucas@lx.it.pt)
%
%
% NMF-QMV is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------





%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------

if (nargin-length(varargin)) ~= 4
    error('Wrong number of required parameters');
end

drawFigs  = 'yes';
%--------------------------------------------------------------
% Read the optional input parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'ENDMEMBER_TRUE'
                endmember_true = varargin{i+1};
            case 'ABUNDANCE_TRUE'
                abundance_true = varargin{i+1};
            case 'DRAWFIGS'
                drawFigs = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end





[lines,samples,bands]=size(img);
np = lines*samples;
Y=[];
Y = reshape(img, np, bands)';



%%
%--------------------------------------------------------------------------
% Optimization
%--------------------------------------------------------------------------


%%
%--------------------------------------------------------------------------
% Dimensionality reduction
%--------------------------------------------------------------------------
my = mean(Y,2); % mean of data set
Yp = Y - repmat(my,1,np);
global svd_D;
[Up,svd_D] = svds(Yp*Yp'/np,p-1);
% represent yp in the subspace R^(p-1)
G = Up'*Yp; %G: Dimension-reduced observed image matrix
if exist('endmember_true','var')
    endmember_truep = endmember_true - repmat(my, 1, p);
    endmember_truep = Up'*endmember_truep;
end







%%
%--------------------------------------------------------------------------
% Boundary Pixel Identification: FastBPI
%--------------------------------------------------------------------------
idx_boundaryP = [];
idx_boundaryP = FastBPI(G);






if strcmp(drawFigs,'yes')
    figure;
    subplot(1,2,1);
    hold on;
    plot(G(1,:),G(2,:),'.');
    plot(G(1,idx_boundaryP'),G(2,idx_boundaryP'),'o' );
    legend('Pixels','Boundary pixels');
    hold off;
end




fprintf('\n\n\n\n');
fprintf('========================\n');
disp(['Regularization term = ',term]);


result_save_Mhat = [];
result_save_Shat = [];
results_save = [];
i_beta = 0;
for  beta_set = beta_candidates
    i_beta = i_beta+1;
    fprintf('\n\n');
    disp(['beta = ',num2str(beta_set)]);
    
    
    
    
    [Mhat,Shat,rerr,L] = ...
        nmf_qmv_subspace(G,p, ...
        'POSITIVITY','yes', ...
        'BETA', beta_set ,...    % minimum volume regularization parameter
        'ADDONE','yes', ...
        'AO_ITERS', 100, ...
        'DELTA', 1e-8,  ...                   % (STOP) relative reconstruction error
        'SUNSAL_ITERS',100, ...
        'MU_A', 1e-4*(np*p), ...       %proximity weight for A
        'MU_X', 1e-1, ...              %proximity weight for X
        'SPHERIZE', 'cov', ...         %{'no','cov', 'M'}
        'MIN_VOLUME', term, ...  %{'boundary', 'center', 'totalVar'}
        'VERBOSE','no', ...
        'PLOT_A', 'no');
    
    
    
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    % Proposed selection criterion: Dist_criterion
    %--------------------------------------------------------------------------
    
    D = [];
    for ip = 1:p %p-th facet
        point_facet = Mhat;
        point_facet(:,ip) = [];
        c = (inv(point_facet))'* (ones(p-1,1)*(-1));
        %plane function: c'*x + 1 = 0
        D(ip,:) = abs(c'* G(1:p-1,idx_boundaryP) + 1)./ norm(c,2);
    end
    D_mm =   min(D) ;
    Dist_criterion = mean(D_mm);
    fprintf('Dist_criterion = %2.2f\n',Dist_criterion);
    
    
    
    
    
    
    if exist('endmember_true','var') %true endmember existing
        
        % find alignment
        Ahat = Up*Mhat;
        Ahat = Ahat + repmat(my, 1, p);
        P = align_matrices(Ahat,endmember_true,'angle');
        
        Ahat = Ahat*P;
        Shat = P'*Shat;
        Mhat = Mhat*P; %update Mhat
        
        
        % compute SAD
        %     aux = sum((Ahat).*endmember_true)./ sqrt(sum((Ahat).^2).*sum(endmember_true.^2));
        %     fprintf('\n\nSAD (deg) = %2.2f\n', mean(abs(acos(aux)))*180/pi);
        
        %normalized mean square error of endmembers
        NMSE_M = (norm((Ahat) - endmember_true,'fro')^2)/(norm(endmember_true,'fro')^2);
        fprintf('NMSE_M = %2.5f\n', NMSE_M);
        
        results_save = [results_save;...
            beta_set,NMSE_M , 0,Dist_criterion];
        
        if exist('abundance_true','var')
            NMSE_S = (norm(Shat - abundance_true,'fro')^2)/(norm(abundance_true,'fro')^2);
            fprintf('NMSE_S = %2.5f\n', NMSE_S);
            
            results_save(end,3) = NMSE_S;
        end
        
        
        
        
        
    else
        results_save = [results_save;...
            beta_set,0 , 0,Dist_criterion];
        
        
        
    end
    result_save_Mhat{i_beta} = Mhat;
    result_save_Shat{i_beta} = Shat;
end

[~,idx_selected] = min(results_save(:,4));

beta_best = beta_candidates(idx_selected);
S_output = result_save_Shat{idx_selected};
M_output = result_save_Mhat{idx_selected};
A_output = Up*M_output + repmat(my, 1, p);

if strcmp(drawFigs,'yes')
    
    subplot(1,2,2);
    hold on;
    plot(G(1,:),G(2,:),'.');
    plot(M_output(1,:),M_output(2,:),'^' );
    str_aux = ['Estimated Endmembers, \beta = ', num2str(beta_best),' , MV=', term];
    legend('Pixels',str_aux);
    hold off;
end

end
