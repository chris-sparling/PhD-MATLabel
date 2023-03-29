function [Data_Abel, I, b_n, C] = pBASEXReconstruct(Data,Basis_Set)

%% pBASEX Reconstruction Function Without Interpolation - Introduction

%{

This function is based on the pBASEX approach, originally described by
Garcia et. al. in their pBASEX publication. The inverse Abel transform of a VMI image is
reconstructed by expanding the image as sum of basis functions with known
Abel inverses. The coefficients of this expansion can be used in
conjunction with these Abel inverted functions to reconstruct the central
slice of the original 3D distribution.

Also using these coefficients, the integrated radial distribution and
anisotropy parameters can also be recovered.

Unlike the original pBASEX code, however, this version does not interpolate
the image onto a polar grid. Rather, it samples the distribution in a
similar manner to that described my Roberts et. al. in their POP routine.

To use this function, a basis set needs to be first constructed using the
newpBASEX_WriteBasisSet.m function.

%}

%% Unpack Basis Set

%{

Note there are additional objects included in the basis set, but not all of
them are strictly needed to reconstruct an image.

%}

G_inv = Basis_Set.G_inv; % Inverse of G-matrix, for which each row is a transformed basis function

res = Basis_Set.res; % Resolution of basis set (each basis functions is [res x res] pixels)

N = Basis_Set.N; % Photon order of the basis set

width = Basis_Set.width; % Width of basis function

spacing = Basis_Set.spacing; % Pixel spacing between adjacent basis functions

F_set = Basis_Set.F_set; % Matrix with each row being an original basis function

G_set = Basis_Set.G_set;

Data_Abel = zeros(size(Data));

I = zeros(1,0.5*res,size(Data,3));

b_n = zeros(N+1,0.5*res,size(Data,3));

%% Transform into polar coordinates, find basis coefficients and create reconstruction

for n = 1:size(Data,3) % for each image in the Data array

Datapol = pol_sample(Data(:,:,n)); % to polar coordinates

Datapol = Datapol(:); % vector rep of polar image

Data_vec = Datapol(~isnan(Datapol)); % remove NaNs from triangular array

C = G_inv*Data_vec; % project image vector onto basis set

Abel = F_set*C; % create inverted image from coefficients

Data_Abel(:,:,n) = reshape(Abel,[res res]); % reshape into Cartesian image

%% Calculate radial distribution and angular parameters

R = linspace(0,0.5*res,0.5*res); % create radial coordinate

b = zeros(N+1,0.5*res); % create space for b-parameters

for k = 1:((0.5*res)/spacing) % for each radial basis
    
    for l = 1:N+1 % for each angular basis
        
        idx = k + (l-1)*((0.5*res)/spacing); % find coefficient index
        
        R_k = (k-1)*spacing; % peak of basis function
            
        b(l,:) = b(l,:) + C(idx,1)*exp(-((R-R_k)/width).^2); % weight Gaussian by coefficient
        
    end
    
end

I(1,:,n) = b(1,:).*R.^2; % intensity is simply P0 basis values weighted by R^2

b_n(:,:,n) = b.*R.^2;

% IF YOU WANT BETAS, YOU NEED TO DIVIDE BY I TO NORMALISE

end

end

