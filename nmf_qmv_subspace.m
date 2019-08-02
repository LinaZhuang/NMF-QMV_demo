function [M,S,rerr,L] = nmf_qmv_subspace(G,p,varargin)

%
%  NMF-QMV working in a dimension-reduced space.
%
%% --------------- Description ---------------------------------------------
%
%
%  Let  -> G [d,n]  matrix containing the dimension-reduced spectral vectors in its
%                    columns
%
%       -> M  [d,p] mixing matrix with p spectral signatures  
%
%       -> S  [p,n] abundance matrix
%
%       -> N  {L,n] additive Gaussian Noise
%
%
%  Linear mixing model
%
%     G = MS + N    with   S>= 0,   and   sum(S) = ones(1,np)
%
% -----------------------------------------------------------------------
%  Optimization problem:
%
%   min  (1/2)  ||G-MS||^2_F + beta ||MB-O||^2_F
%    S >= 0
%
%    OPTIONAL CONSTRAINTS:
%
%       Sum-To-One sum(S) = ones(1,np)
%
%   where the couple (B,O) define a minimum volume-type  regularizer:
%   
%   'MIN_VOLUME' - {'boundary', 'center', 'total variance'}
%   
%      'boundary' - > B=I, O = (VCA solution)
%    
%      'center' - >   B=I, O = (sample mean value)
%
%      'totalVar'     B = B = eye(p) - ones(p)./p, O = 0
%           
%       In the the three above cases,  the term beta ||MB-O||^2_F
%       pushes the columns of M  towards the simplex thus having a minimum 
%       volume flavor
%
%
%% -------------------- Line of Attack  -----------------------------------
%
%  NMF_QMV implemets proximal alternating sptimization  (PAO). 
%
% ------------------------------------------------------------------------
% 
% See more details in papers:
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and Jose M. Bioucas-Dias,
% "Regularization Parameter Selection in Minimum Volume Hyperspectral Unmixing",
% TGRS, 2018.
%
%%  ======================  Required inputs ====================== 
%
% G - matrix with  n(channels) x d(pixels).
%     each dimension-reduced pixel is a linear mixture of p dimension-reduced endmembers
%     signatures g = M*s + noise,
%
% p -  the number of endmembers
%
%
%%  ====================== Optional inputs =============================
%
%  'POSITIVITY'  = {'yes', 'no'}; Enforces the positivity constraint:
%                   S >= 0
%                   Default 'yes'
%
%  'ADDONE'  = {'yes', 'no'}; Enforces the positivity constraint: S >= 0
%              Default 'no'
%
%  'BETA'  - regularization parameter for the minimum volume term
%
%  'MU_A' -  weight for A  proximal term
%
%  'MU_S' -  weight for S  proximal term
%
%  'AO_ITERS' -  number of iterations
%
%  'DELTA' -  relative reconstruction error. Used as STOP criterion
%
%  'THETA' -  thresold to detect the active rows of S
%             default = 1e-3*n;
%
%  'SUNSAL_ITERS' - CSUNSAL number of iterations
%                    Default: 100;
%
%  'SPHERIZE' -  {'no','cov', 'M'} sperize the observed data Y
%                'cov'  based on the covariance matrix
%                'M'    based on an estimated mixing matrix
%   
%  'MIN_VOLUME' - {'boundary', 'center', }  define the minimum volume 
%                  term ||AB-O||^2_F
%
%                 'boundary'    B= I, O =  extremes given by VCA
%                 'center'      B=I,  O = center of mass
%                 'totalVar'    B = eye(p) - ones(p)./p, O = 0
%
%  'A0' -  initial A
%
%  'TRUE_M' - accept the true mixing matrix (only for displaying purposes)
%
%  'VERBOSE'   = {'yes', 'no'};
%                 'no' - work silently
%                 'yes' - display warnings and plot figures
%                  Default 'no'
%
%
%
%%  =========================== Outputs ==================================
%
% M  =  [d x p] estimated mixing matrix
%
% S  =  [p x n] estimated abundances
%
% rerr ->  reconstruction error along the iterations
%
% L   ->   objective fuction along the iterations
%
%% -------------------------------------------------------------------------
%
% Copyright (Oct. 2018):     José Bioucas-Dias (bioucas@lx.it.pt)
%                            %
%                            Lina Zhuang (lina.zhuang@lx.it.pt)
%
% CoNMF is distributed under the terms of
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
% ----------------------------------------------------------------------
%
% -------------------------------------------------------------------------

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
 
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% data set size
[d,n] = size(G);
 
%%
%--------------------------------------------------------------
% Set the default for the optional input parameters
%--------------------------------------------------------------
%positivity
positivity = 'yes';
%addone
addone = 'no';
%regularization parameters
lam_a = 5e-1*(n/1000);      % minumum volume regularization parameter
% parameter for the quadratic proximal term for A
mu_a = 1e0*(n/1000);
% parameter for the quadratic proximal term for S
mu_s = 1e0;
% maximum number of AO iterations
AOiters = 100;
% relative reconstruction error (RECONS_ERR) for the termination test
delta = 1e-4;
% thresold to detect the active rows of S
theta = 1e-3*n;
%  CSUNSAL number of iterations
sunsal_iters = 100;
% true A
M_true = [];
% display only CoNMF warnings
verbose = 'no';
%spherize
spherize = 'M';
% minimum volume tetrm
min_volume = 'boundary';
% no initial simplex
A = [];
% use csunsal_hinge instead of csunsal
hinge = 0;
% plot evolution of A(t)
plot_a_t = 'no';

%--------------------------------------------------------------
% Read the optional input parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'POSITIVITY'
                positivity = varargin{i+1};
            case 'HINGE'
                hinge = varargin{i+1};
            case 'ADDONE'
                addone = varargin{i+1};
            case 'BETA'
                lam_a = varargin{i+1};
            case 'MU_A'
                mu_a = varargin{i+1};
            case 'MU_X'
                mu_s = varargin{i+1};
            case 'AO_ITERS'
                AOiters = varargin{i+1};
            case 'DELTA'
                delta = varargin{i+1};
            case 'THETA'
                theta = varargin{i+1};
            case 'SUNSAL_ITERS'
                sunsal_iters = varargin{i+1};
            case 'TRUE_M'
                M_true = varargin{i+1};
            case 'SPHERIZE'
                spherize = varargin{i+1};
            case 'MIN_VOLUME'
                min_volume = varargin{i+1};
            case 'A0'
                A = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};            
            case 'PLOT_A'
                 plot_a_t = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%%
%--------------------------------------------------------------
% take care of output variable
%--------------------------------------------------------------

% compute objective function
comp_obj_func = 0;

if nargout == 4
    comp_obj_func = 1;
end

%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------

% regulatization parameter for spherization
lam_sphe_reg  = 1e-6;
% VCA number of VCA runs
VCA_runs = 30;

%%
  

%%
%------------------------------------------
% spherize if requested
%------------------------------------------
G_aux = [G; repmat(1,1,size(G,2))];
if strcmp(spherize,'cov')
    global svd_D;
    sing_values = diag(svd_D);
    C = diag([1./(sqrt(sing_values(1:p-1))+lam_sphe_reg)]);
elseif strcmp(spherize,'M')
    % estimate M with VCA
    Aux = VCA(G_aux,'Endmembers',p,'SNR',1,'verbose',verbose);
    C = inv(Aux + lam_sphe_reg*eye(p));
else
    C = eye(p-1);
end

% spherize data
G = C*G;
% represent A in the new coordinate system
if ~isempty(A)
   M = C*A;
end
% represent M_true in the new coordinate system
if ~ isempty(M_true)
   M_true = C*M_true;
end





%%
%---------------------------------------------
%  Initialization
%---------------------------------------------
small = 1e-4;
vol_ant = -inf;

% Initialize with VCA
G_aux = [G; repmat(1,1,size(G,2))];
for i=1:VCA_runs
    Aux = VCA(G_aux,'Endmembers',p,'SNR',1,'verbose',verbose);
      vol = sum(log(svd(Aux) + small ));
    if vol > vol_ant
        Avca = Aux;
        vol_ant = vol;
    end
end
 Mvca =Avca(1:p-1,:);
 clear Avca;


if isempty(A)
    M = Mvca;
end



% initial abundances
S = sunsal(M,G,'POSITIVITY',positivity,'VERBOSE',verbose,'ADDONE',addone, ...
    'lambda', 1e-3,'AL_ITERS',1000, 'TOL', 1e-8, 'verbose',verbose);

 

% A is represented in the affine set as A = Yc + F*W
% Yc = repmat(yc,1,p);

%set the couple(B,O) for the minimum volume term
if strcmp(min_volume, 'boundary')
    %use VCA extremes  
 O = Mvca;
 B = eye(p);
elseif strcmp(min_volume, 'center')
    % mass-center
 O = zeros(p-1, p);
    B = eye(p);
else  % total variance
 O = zeros(p-1, p);
    B = eye(p) - ones(p)./p;
end



% use vca or mass-center as points close to vertices to define the 
% minimum volume term

 
if strcmp(plot_a_t, 'yes')
    figure(1000);
    hold off;
    plot(G(1,:), G(2,:),'.');
    hold on
    plot(Mvca(1,:), Mvca(2,:),'go','LineWidth',3);
    title(['beta = ',num2str(lam_a)]);
    if ~ isempty(M_true)
        plot(M_true(1,:),M_true(2,:)','ro','LineWidth',3);
    end
end


%---------------------------------------------
%   Constants
%---------------------------------------------

const1 = lam_a*O*B';
const2 =  lam_a*B*B'+mu_a*eye(p);


% inilialize rerr and L
% norm of Y
norm_y = norm(G,'fro');
% initial reconstruction errors
rerr(1) = norm(G-M*S,'fro');

L(1)=0.5*norm(G-M*S,'fro')^2  + 0.5*lam_a*norm(M*B-O,'fro')^2;
             
%---------------------------------------------
%   main body
%---------------------------------------------       
t=1;
k=1;

while (t <= AOiters) && (rerr(t)/norm_y > delta)



    % update M (dimensional reduced endmembers)
  M_aux1 = G*S'+const1+mu_a*M;
  M_aux2 = S*S'+const2;
  M_aux2 = inv(M_aux2);
  M = M_aux1*M_aux2 ;
   
    
    rerr(t+1) = norm(G-M*S,'fro');
    
    
    if comp_obj_func
        L(t) = 0.5*norm(G-M*S,'fro')^2  + 0.5*lam_a*norm(M*B-O,'fro')^2;
    end
    
    % update S
    G_aux = [G; sqrt(mu_s)*S];
    M_aux = [M; sqrt(mu_s)*eye(p)];
    
    S = sunsal(M_aux,G_aux,'POSITIVITY',positivity,'VERBOSE',verbose,'ADDONE',addone, ...
        'lambda', 0,'AL_ITERS',sunsal_iters, 'TOL', 1e-6, ...
        'X0',S,'verbose', verbose);    
 
    
    if strcmp(plot_a_t, 'yes')
        figure(1000);
        plot(M(1,:), M(2,:),'k.','LineWidth',3);
        drawnow;
    end
    t = t+1;
    
end % main body



% undo spherization and back to the original representation
M = inv(C)*M;




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
