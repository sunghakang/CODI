%%%% Counting Objects by Diffused Index:geometry-free and training-free approach %%%%
%%% This script is about counting the number of objects in a given image by
%%% using CODI-S
%% PSEUDO-CODE
%{
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


%}



clc
clear
close all

%% INITIALIZE
DownSampleRate     = 0.5;     % The percentage of pixels kept from an original image
fig_name           = 'sample_img.jpg'; % The name of the figure
half_seed_len      = 1;       % half of the length of a squared seed
dist_seed_center   = 6;       % distance between the center of consecutive seeds
mu                 = 0.00005; % PARA mu in diffusion algorithm
theta              = 1;       % PARA theta in diffusion algorithm
eta                = 0.0001;  % PARA eta in diffusion algorithm
R_n                = 0.2;    % Stoppting criteria of diffusion algorithm
treshvalue         = 205;     % Threshold in preprocessing step to get mask image
r                  = 8;       % PARA in Gaussian convolution curve
sigma              = 0.75;    % PARA in Gaussian convolution curve



% read an image
Image = imread(fig_name);
disp('Finished reading');

tic
disp('Start timing')
% figure;imshow(Image);title('original image')

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
figure;imshow(phi,'border','loose'); title('Thresholding the gray scale image')


% generate the descretized g function/ a mask image
EdgeMask = phi>0;

%%   generate seeds in 1D

seed_image = zeros(n,m);
seed_binary_domain=zeros(n,m);


index_i = half_seed_len+1:dist_seed_center:n;
index_j = half_seed_len+1:dist_seed_center:m;


% generate vertically increasing seeds in 1st component
k = 1;
for i = index_i
    for j = index_j
        seed_image(i-half_seed_len:i,j-half_seed_len:j,1) = k;
        seed_binary_domain(i-half_seed_len:i,j-half_seed_len:j) = 1;
        k=k+1;
    end
end

seed_image = seed_image./max(seed_image,[],'all').*255;
% Ensure that the values of seeds are located in [1,255] % to distinct the
% seed values from the background value 0

fprintf('Finished generating seeds \n')
%%  Edge-weighted harmonic variational optimization model

[AlgIdiffIm,AlgIobj1,AlgIiter1,AlgICPU1] = AlgIEdge(seed_image,EdgeMask,seed_binary_domain,mu,theta,eta,R_n);
AlgIdiffIm = normalize(AlgIdiffIm).* EdgeMask;

% test the result from bwconncomp
% a =  bwconncomp(phi./255);
% mablabobjectnumber = a.NumObjects

U = AlgIdiffIm;

fprintf('Finished diffusing seeds \n')

%% Apply Gaussian Fitting curve to Diffused image
U0 = U(U>0);

binwidth= 0.1; edges = 0:binwidth:256;
fig=figure; set(fig,'visible','off');
h = histogram (U0,edges,'Normalization','pdf');
f = h.Values; f(f<0.005)=0;

figure(1);
hb= bar((binwidth:binwidth:256),f); hb.FaceColor = 'r';
hold on

g  = normpdf(-r:r,0,sigma); p1=conv(f,g);
pg = plot(-(r-1)*binwidth: binwidth:256+r*binwidth,p1,'k--', 'Marker','.');
xlabel('Pixel Intensity'); xlim([53 75])
legend ([hb;pg], 'histogram', 'Gaussian Fitted')
hold off
[pks, locspks] = findpeaks(p1);
NumOfClusters = length(pks);
title([num2str(NumOfClusters),' objects are found'])

t = toc;
fprintf('Finished clustering \n')
fprintf(['CPU                       : ',num2str(t),'\n'])
fprintf(['Number of Clusters found  : ', num2str(NumOfClusters),'\n'])
fprintf(['Ground truth              : ', num2str(30),'\n'])

