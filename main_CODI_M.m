%%%% Counting Objects by Diffused Index:geometry-free and training-free approach %%%%
%%% This script is about counting the number of objects in a given image by
%%% using CODI-M

%% PSEUDO-CODE
%{

INPUTS
------
DownSampleRate      : The percentage of pixels kept from an original image
fig_name            : Directary of a fig 
half_seed_len       : half of the length of a squared seed
dist_seed_center    : distance between the center of consecutive seeds
eps                 : PARA eps in DBSCAN
MP                  : PARA m in DBSCAN
mu       = 0.00005  : PARA mu in diffusion algorithm
theta    = 1;       : PARA theta in diffusion algorithm
eta      = 0.0001;  : PARA eta in diffusion algorithm
R_n      = 0.01;    : Stoppting criteria of diffusion algorithm
seedsdim = 4        : Total dimension of seeds
treshvalue          : Threshold in preprocessing step to get mask image


OUTPUTS
-------
NumOfClusters       : Total Number of clusters found


FUNCTIONS:
-------
normalize(x)
    -A helper method to return the normalized data of input in [0,255]
    -inputs:
     x            : Input matrix of size m, n
    -output :
     u            : Output matrix of size m, n, values in [0,255]

add_frame(phi,k, value)
    -A helper method to add a frame of some value and width k into a matrix
    -inputs:
     phi          : a gray-scale image of size M * N
     k            : the width of the frame
     value        : the value of the boundary

    -output:
     framed_phi   : a grapy-scale image of size (M+k) * (N+k)

Treshold2BW_up (x,t)
    -A helper method to return a binary matrix of value 0 and 255 of the 
     input matrix
    -inputs:
     x            : a matrix of size M * N
     t            : threshold
    -output:
     u            : a matrix of size M * N, with value 0 where x(i,j) < t
                                            with value 255 o.w.

PlotClusterinResult(D, IDX, dir1, dir2)
    -A helper method to plot the clustering results from 4 dimensional data
     into a subspace contructed by dir1 and dir2 dimension
    -inputs:
     D            : a matrix of size m*n, m data points in n dimensions
     IDX          : a list of cluster index of each point in D
     dir1, dir2   : the index of the subspace, default values are 1 and 2

    -output:
     Figure       : projection of clustering results from high dimension to
                    low dimension

AlgIEdge(SeedIm,EdgeIm,X_D,mu,theta, eta, perc)
    -A helper function about the Edge-weighted harmonic variational optimization model
     proposed in the paper to produce a diffused image in a single
     dimension

    -inputs:
     SeedIm       : Seed Image
     EdgeIm       : Mask Image of the image that $g$ applied onto
     mu           : multiplier on the constraint V−U=0.
     lambda       : λ, the Lagrange multiplier associated with the linear constraint V−U=0.
     eta          : η, fidelity parameter

    -output :
     DiffiusedIm  : Diffused image
     objval       : Energy Function
     Itercount    : Iteration Number
     CPUcount     : CPU time

DBSCAN(D,eps,MP)
    -A helper method to cluster data points in a high dimension cloud by
    using a density-based clustering algorithm DBSCAN
    -inputs:
     D            : a matrix of size m*n, m data points in n dimensions
     eps          : epsilon, PARAMETER in DBSCAN
     MP           : MinPoints, PARAMETER in DBSCAN

    -output:
     IDX          : A list of m*1 data, ith element represent the index of
                    the cluster ith row in D

CODI_M_visualization.m
    Visualize the segmentation results of CODI-S on original image
    pixels that is clustered into a single objects are assigned to a unidex
    color. The noise pixels identified in DBSCAN has the same color as the
    background

%}

clc
clear
close all

%% INITIALIZE
DownSampleRate     = 0.5;     % The percentage of pixels kept from an original image
fig_name           = 'sample_img.jpg';
half_seed_len      = 1;       % half of the length of a squared seed
dist_seed_center   = 6;       % distance between the center of consecutive seeds
eps                = 2;       % PARA eps in DBSCAN
MP                 = 6;       % PARA m in DBSCAN
mu                 = 0.00005; % PARA mu in diffusion algorithm
theta              = 1;       % PARA theta in diffusion algorithm
eta                = 0.0001;  % PARA eta in diffusion algorithm
R_n                = 0.01;    % Stoppting criteria of diffusion algorithm
seedsdim           = 4;       % Total dimension of seeds
plot_gt            = 1;       % Whether to plot ground truth
treshvalue         = 205;      % Threshold in preprocessing step to get mask image

% read an image
Image = imread(fig_name);

disp('Finished reading');

tic
disp('Start timing')
figure;imshow(Image);title('original image')

% Convert image from unit8 into double
Im = double(Image);

% Dowmsample the origninal image
% DownSampleRate : # pixels in the downsampled image / # pixels in the
%                  original image
Im    = imresize(Im,DownSampleRate);

% transfer the image from color scale into gray scale
Im    = sqrt(sum(Im.^2,3));
phi   = Im;

% normalize phi to have value between 0 and 255
phi = normalize(phi);

% add a frame
phi   = add_frame(phi, 3, 0);
[n,m] = size(phi);



%% Preprocess the gray_scale image

% identify background by using threshvalue
phi = Treshhold2BW_up(phi,treshvalue); % phi = 255 if phi < threshvalue;
% figure;imshow(phi,'border','loose'); title('Thresholding the gray scale image')

% generate the descretized g function/ a mask image
EdgeMask = phi>0;
%%   generate seeds in nD

seed_image         = zeros(n,m,seedsdim);
seed_binary_domain = zeros(n,m);

seed_loc_i = half_seed_len+1:dist_seed_center:n; % x - seeds centers
seed_loc_j = half_seed_len+1:dist_seed_center:m; % y - seeds centers

% generate seeds in the 1st component
% increse j first then increase i
k = 1; % start with value 1 and every time increase k by 1
for i= seed_loc_i
    for j= seed_loc_j
        seed_image(i-half_seed_len:i,j-half_seed_len:j,1) = k;
        seed_binary_domain(i-half_seed_len:i,j-half_seed_len:j) = 1;
        k=k+1;
    end
end


% generate seeds in the 2nd component
% increse i first then decrease j
k = k - 1;
for  j= flip(seed_loc_j)
    for i= seed_loc_i
        seed_image(i-half_seed_len:i,j-half_seed_len:j,2)=k;
        k=k-1;
    end
end

% generate permutated seeds in 3rd, 4th components

% generate random permutated uniformly distributed values of seeds
% between [1,255] and store them in random_seeds_val.
% each column is a permutation

random_seeds_val = [];
for p = 1:2
    rng(p);
    random_seeds_val (:,p) = randperm(length(seed_loc_i)*length(seed_loc_j));
end

l = 1;
for i= seed_loc_i
    for j= seed_loc_j
        seed_image(i-half_seed_len:i,j-half_seed_len:j,3) = random_seeds_val(l,1);
        seed_image(i-half_seed_len:i,j-half_seed_len:j,4) = random_seeds_val(l,2);
        l=l+1;
    end
end

% Ensure that the values of seeds are located in [1,255] % to distinct the
% seed values from the background value 0
seed_image = seed_image./max(seed_image,[],'all').*255;
fprintf('Finished generating seeds \n')
%%  Edge-weighted harmonic variational optimization model

% perform diffusion in each dimension 
for p = 1: seedsdim
    [AlgIdiffIm(:,:,p),AlgIobj1,AlgIiter1,AlgICPU1] =  AlgIEdge(seed_image(:,:,p),EdgeMask,seed_binary_domain,mu,theta,eta,R_n);
    % normalize the diffused image and only preserve the pixels within
    % the mask image
    AlgIdiffIm(:,:,p) = normalize(AlgIdiffIm(:,:,p)).*EdgeMask;
end

fprintf('Finished diffusing seeds \n')
%% Apply DBSCAN to diffused image
% Vectorize the diffused pixels in all objects.
D = [];
for p = 1:seedsdim
    DD =  AlgIdiffIm(:,:,p);
    D(:,p) = DD(DD >0);
    % Each row in D is associated with a pixel in the mask in multiple
    % dimension
    % Each column of D is associated with all pixels in the mask in one
    % dimension
end

% Apply a high dimensional density based clustering algorithm to D with
% predeterminated paramters eps(\epsilon) and MP(MinPoints)

IDX = DBSCAN(D,eps,MP);
NumOfClusters = max(IDX);

% End counting CPU time
t = toc;
fprintf('Finished clustering \n')
fprintf(['CPU                       : ',num2str(t),'\n'])
fprintf(['Number of Clusters found  : ', num2str(NumOfClusters),'\n'])
fprintf(['Ground Truth              : ', num2str(30),'\n'])
% Projecting counting results onto original image
run CODI_M_visualization.m

