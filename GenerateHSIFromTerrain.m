%%
% In order to have a quantitative unmixing evaluation for NMF-QMV algorithm , 
% we simulate a semi-real image base on the publicly available
% TERRAIN hyperspectral image acquired by the HYDICE sensor.  
%
% This script generates a simulated image from TERRAIN data
%
%     y = M*S+noise
%
% where M is the mixing matrix containing the endmembers, S is the
% the abundance of each enmember at each pixel, and noise is a
% Gaussian independent (band and pixel - wise) additive perturbation.
% 
%
%Simulated TERRAIN image is used to evaluate non-pure-pixel based unmixing
%method NMF-QMV in the paper:
%
% Lina Zhuang, Chia-Hsiang Lin, Mario A.T. Figueiredo, and 
% Jose M. Bioucas-Dias, "Regularization Parameter Selection in Minimum 
% Volume Hyperspectral Unmixing", TGRS, 2019.
%
%
% If you use this code, please cite the above paper. 
% For the latest version of the code and papers, please check at
% http://www.lx.it.pt/~bioucas/publications.html
%
%
%
% 
% Author: Lina Zhuang, Oct., 2018.
% -----------------------------------------------------------------------
% Copyright (2018): Lina Zhuang
% 
% 'GenerateHSIFromTerrain' is distributed under the terms of 
% the GNU General Public License 2.0.
% ----------------------------------------------------------------------
% 
%  =====================================================================
%  ====== important variables: ==========================================
% M <-- endmember matrix
% S_low_pass <-- simulated abundance matrix
% img_syn <-- simulated clean image
% noise <-- nosie matrix
% img_syn_noisy <-- simulated noisy image
 
 

clear;clc;
close all;
rng(100);


load terrain_envi_original_b210;
img = terrain_envi; clear terrain_envi;
%img = img./1000;
[r, c, b] = size(img);
N = r*c;
Y = reshape(img, N, b)';clear img;

%remove 44 low SNR bands:
[w Rw] = estNoise(Y,'additive','on');
%figure;plot(diag(Rw));
sigma_noise_oriImg = sqrt(diag(Rw));
xlabel('Band');ylabel('noise variance');title('Noise level in original TERRAIN image');
[band_lowSNR,band_idx] = sort(diag(Rw),'descend');
Y(band_idx(1:44),:) = [];
sigma_noise_oriImg(band_idx(1:44,1),:) = [];
band_idx = sort(band_idx(45:end),'ascend'); %band_idx records indexes of remaining bands


p = 5; %number of endmebers
b = size(Y,1); %update the number of bands



%% endmember extraction:  NFINDR
%Setp (a): five endmembers (namely, Soil 1, Soil2, Tree, Shadow, Grass) are 
% extracted by the N-FINDR algorithm;  
[Endm_nfindr,location_mask] = EIA_NFINDR(Y,p, 1000);
M = Endm_nfindr;


%show endmembers
figure;
plot(band_idx,M,'-o');
legend('1','2','3','4','5','6');
xlabel('Band'); title('Endmembers estimated by VCA');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%unmixing: FCLS
%Setp (b): the original TERRAIN image is unmixed by FCLS, obtaining five
%abundance maps with respective to estimated endmembers;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The unmixing problem is solved by sparse unmixing by variable splitting
% and augmented Lagrangian (SUnSAL) algorithm  by disabling the sparsity
% regularization (i.e., 'lambda'=0).

S = sunsal(M,Y,'POSITIVITY','yes','VERBOSE','no','ADDONE','yes', ...
    'lambda', 0,'AL_ITERS',100, 'TOL', 1e-6, ...
    'verbose', 'no');

% for ip=1:p
%     max_S(ip) = max(S(ip,:));
%     min_S(ip) = min(S(ip,:));
% end


S_img = reshape(S', r,c,p);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% low pass filter (no pure pixels)
% max(s)<0 <- no pure pixels
% min(s)=0 <- we have pixels lying on all facets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step (c) in order to simulate
% scenarios without pure pixels, we convolve each abundance
% map using a low-pass filter with a 2-D Gaussian smoothing
% kernel with standard deviation of 2 and finally maximum of
% each convolved abundance map is less than one (implying
% no pure pixel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

for ip=1:p
    sigma_imgaussfilt =  2;
%     if (ip == 3) || (ip ==4)
%         sigma_imgaussfilt =  2;
%     end
    S_img_low_pass(:,:,ip) = imgaussfilt(S_img(:,:,ip),sigma_imgaussfilt);
    
    subplot(2,6,ip)
    imagesc(S_img(:,:,ip),[0,1]);
    title(['Original abundances  ',num2str(ip)]);
    subplot(2,6,ip+6)
    imagesc(S_img_low_pass(:,:,ip),[0,1]);
    title(['Convolved abundances  ',num2str(ip)]);
end



% check abundance distrubtion
for ip=1:p
    max_S(ip) = max(max(S_img_low_pass(:,:,ip)));
    min_S(ip) = min(min(S_img_low_pass(:,:,ip)));
end
% max_S
% min_S



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step (d): in order to simulate sparse abundances,
% we sort the elements of each abundance map in ascending
% order, followed by setting the first 1% of elements as 0;
%
%  min(s)=0 <- we have pixels lying on all facets of endmeber simplex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% for ip = 1:p
%  
%         s_aux = S_img_low_pass(:,:,ip);
%         s_aux = s_aux(:);
%       s_aux(find(s_aux<0.01)) = 0;
%       S_img_low_pass(:,:,ip) = reshape(s_aux, r, c);
% end

for ip = 1:p
   % if min_S(ip)>0
        s_aux = S_img_low_pass(:,:,ip);
        s_aux = s_aux(:);
        [~, idx_sort] = sort(s_aux, 'ascend');
        idx_selected = idx_sort(1:fix(N*0.01));
        s_aux(idx_selected) = 0;
        S_img_low_pass(:,:,ip) = reshape(s_aux, r, c);
    %end
end


% %%  max(s)<0 <- no pure pixels
% for ip=1:p
%     S_aux = S_img_low_pass(:,:,ip);
%    if max_S(ip)>0.7
%        S_aux( find(S_aux>0.7) ) = 0.7;
%        S_img_low_pass(:,:,ip) =  S_aux;
%    end
% end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step (e):
% abundance vector of each pixel is normalized to satisfy sum-toone
% property;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalization:
S_low_pass = reshape(S_img_low_pass, N, p)';
S_aux = sum(S_low_pass);
S_aux = repmat(S_aux, p, 1);
S_low_pass = S_low_pass./S_aux;


%  figure;
% % check abundance distrubtion
% for ip=1:p
%     max_S(ip) = max(S_low_pass(ip,:));
%     min_S(ip) = min(S_low_pass(ip,:));
%     subplot(1,p,ip);
%     histogram(S_low_pass(ip,:));
% end
% max_S
% min_S





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step (f): a clean HSI is synthesized based on the linear
% mixing model, i.e. M*S_low_pass, where M is endmembers extracted in
% step (a) and S_low_pass is abundances generated in step (e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% synthesize
Y_syn = M*S_low_pass;
img_syn = reshape(Y_syn', r, c, b);





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step (g):
% Gaussian i.i.d. noise is added to simulate a noisy TERRAIN HSI, i.e.
% img_syn_noisy = M*S_low_pass +noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adding Gaussian noise with variances equal to those estimated from
% each band of the original TERRAIN using HySime
% img_syn_noisy=img_syn+noise;


for ib = 1:size(img_syn,3)
     noise(:,:,ib) = sigma_noise_oriImg(ib)*randn(size(img_syn(:,:,1)));

%       noise(:,:,ib) = sigma_noise_oriImg(ib)*10*randn(size(img_syn(:,:,1)));
%  sigma_noise_oriImg(ib) = rand(1); 
% noise(:,:,ib) = 20* sigma_noise_oriImg(ib)*randn(size(img_syn(:,:,1)));
end
img_syn_noisy=img_syn+noise;
disp(['SNR = ',num2str(snr( img_syn(:), noise(:))), ' dB']);

save SimulatedTerrainDataset.mat img_syn noise img_syn_noisy  M S_low_pass p;