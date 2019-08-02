
function [idx_boundaryP] = FastBPI(X, varargin) 
%
%  Fast boundary pixel learning (FastBPI)
% 
%   ====================== Required input ======================
%
%   X [p-1,np] :  Dimension-reduced observed image matrix containing 
%                      the spectral vectors in its columns
%
%
%    ====================== Optional inputs =============================
%
%   K: number of iterations
%   n:  dimension of subspace
%
%   =========================== Outputs ==================================
%
%   idx_boundaryP : indexes of detected boundary  pixels
%
%
%  ===================See more details in the paper:======================
%  Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and 
%  Jose M. Bioucas-Dias, "Regularization Parameter Selection in Minimum 
%  Volume Hyperspectral Unmixing", TGRS, 2019.
%  -------------------------------------------------------------------------
%
% Copyright (Oct., 2018):     Chia-Hsiang Lin (chiahsiang.steven.lin@gmail.com)
%
% CoNMF is distributed under the terms of
% the GNU General Public License 2.0.

%% 
rng(1);  
%--------------------------------------------------------------
% Set the default for the optional input parameters
%--------------------------------------------------------------

K = 100;
n = 5;

%--------------------------------------------------------------
% Read the optional input parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional input parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'n'
                n = varargin{i+1};
            case 'K'
                K = varargin{i+1};
        end
    end
end



idx_boundaryP = [];
p = size(X,1) + 1; %number of endmembers
if (p-1)<=n
    Xidx = convhulln(X');
else
    
    
    vec = @(x)   x(:);
    %  boundary learning
    
    Xidx=[];
    for j=1:K
        
        [Q,R] = qr(randn(p-1,p-1));
        R = diag(diag(R)./abs(diag(R)));
        U = Q*R;
        
        
        pro=U(:,1:n);
        Xp=pro'*X;
        Xi=convhulln(Xp');
        Xidx=union(Xidx,vec(Xi));
        
    end
end

idx_boundaryP  = union(Xidx,Xidx);

