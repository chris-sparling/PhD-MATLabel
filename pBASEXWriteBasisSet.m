function [Basis_Set] = pBASEXWriteBasisSet(res,N,width,spacing)

%% Writing pBASEX Basis Sets - Introduction

%{

This function is based on the pBASEX approach, originally described by
Garcia et. al. in their pBASEX publication. The inverse Abel transform of a VMI image is
reconstructed by expanding the image as sum of basis functions with known
Abel inverses. The coefficients of this expansion can be used in
conjunction with these Abel inverted functions to reconstruct the central
slice of the original 3D distribution.

Gaussian-like radial basis sets are used here, but this is ultimatley a
users choice. Sech-squared functions, Lorentzians etc. should all behave
similarly. If modifying the form of radial function here, be sure to also
modify the radial functions in the pBASEXReconstruct.m function, too.

Rather than numerically integrating to find the Abel transform of each
basis function, we use the matrix method of Livingstone et. al. to make use of simple 
matrix multiplication to minimise computation time.

%}

%% Create basis space and variables

t0 = tic;

F = zeros(res,res,(0.5*res)/spacing,N+1); % space for original functions

G = F; % space for Abel transformed functions

% Create variables space

x = linspace(-0.5*res,0.5*res,res); 

[X, Y] = meshgrid(x,x); % create Cartesian grid

R = sqrt(X.^2 + Y.^2); % define corresponding polar grid

T = atan2(X,Y);

% Write polynomials

leg_P = zeros(res,res,N+1); % create space for polynomials

leg_P(:,:,1) = ones; % P0

leg_P(:,:,2) = cos(T); % P1

for l = 2:N % for P2 to PN
    
    n = l - 1;
    
    leg_P(:,:,l+1) = (1/(n+1))*((2*n+1)*leg_P(:,:,2).*leg_P(:,:,l) - n*leg_P(:,:,l-1)); % use recursion formula
    
end

% Write basis functions

for k = 1:((0.5*res)/spacing)
    
    R_k = (k-1)*spacing;
    
    R_basis = exp(-((R - R_k)/width).^2); % create radial basis
    
    for l = 1:N+1
        
        F_working = R_basis.*leg_P(:,:,l); % modify with angular part
        
        F(:,:,k,l) = F_working/max(max(F_working)); % store in F basis set
        
        G_working = Abel(F_working);
        
%         if mod(l,2)~=0
%             
%             G_working = Abel(F_working')';
%             
%         end
%         

        G(:,:,k,l) = G_working; %store in G basis set
        
    end
    
end

% Functions written!

Basis_Set_Written = toc(t0)

%% Transform each G basis function to polar coordinates and save result in G matrix

t0 = tic;

r = [0:0.5*res];

n = floor(pi*(r));

G_set = zeros([sum(n), ((0.5*res)/spacing)*(N+1)]); % matrix of G functions

F_set = zeros([res^2, ((0.5*res)/spacing)*(N+1)]); % matrix of F functions (note different sizes)

for k = 1:((0.5*res)/spacing)
    
    for l = 1:N+1
        
        idx = k + (l-1)*((0.5*res)/spacing);
        
        G_sample = pol_sample(G(:,:,k,l)); % to polar using POP method
        
        G_sample = G_sample(:); % to vector
        
        G_sample = G_sample(~isnan(G_sample)); % remove NaNs
        
        G_set(:,idx) = G_sample; % store in G matrix
        
        F_working = F(:,:,k,l); % lookup corresponding Cartesian F function
        
        F_set(:,idx) = F_working(:); % store in F matrix
        
    end
    
end

G_inv = pinv(G_set); % Moore-Penrose inverse via SVD - very fast!

% Inverted!

G_inverted = toc(t0)
    
% Write everything to output

Basis_Set.F = F; % unprojected basis functions

Basis_Set.G = G; % projected basis functions

Basis_Set.F_set = F_set; % matrix of unprojected basis functions

Basis_Set.G_set = G_set; % matrix of projected basis functions

Basis_Set.G_inv = G_inv; % inverse of above

Basis_Set.res = res; % square resoltuion

Basis_Set.N = N; % max legendre polynomial

Basis_Set.width = width; 

Basis_Set.spacing = spacing;

Basis_Set.Basis_Set_Written_t = Basis_Set_Written;

Basis_Set.G_inverted_t = G_inverted;

end

