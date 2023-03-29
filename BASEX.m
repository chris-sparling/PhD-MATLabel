function [Data_Abel_BASEX] = BASEX(Data, sym, tol)

%% BASEX Function - Introduction

%{

The BASEX method expands an image as a sum of well-behaved functions which
have a known Abel inverse. The calculated expansion coefficients can then
be used to reconstruct the Abel inverse of the original image by instead
expanding with the Abel inverted basis functions.

This method was first introduced by Dribinski et. al. 
See https://doi.org/10.1063/1.1482156 for details.

The basis functions used here are Gaussians, and there Abel transform is
numerically calculated using the matrix Abel transform.


%}

%% Write basis functions

%{

This entire section can be prelocated and saved ahead of time for optimal
speed. This section creates basis functions in the projected and Abel
inverted regimes.

%}

tic

w = length(Data)/2;

A = zeros(w,w);

for r=1:length(A)          %for each radius
    
    for y=1:length(A)       %for each slice of y
        
        if y==r     %if your touching the y axis
           A(r,y)=-0.5*sqrt(-1+2*r)*r+0.5*sqrt(-1+2*r)-0.5*r^2*asin((r-1)/r)+0.25*r^2*pi;
           
        elseif y<r  %if below diagonal do nothing 
            
        else        %if above diagonal
            A(r,y)=(-0.5*y^2+y-0.5)*asin(r/(y-1)) + (0.5*y^2-y+0.5)*asin((r-1)/(y-1))...
                + 0.5*sqrt(y^2-2*y+2*r-r^2)*r - 0.5*sqrt(y^2+2*r-1-r^2)*r - 0.5*sqrt(y^2-2*y+2*r-r^2)...
                + 0.5*sqrt(y^2+2*r-1-r^2) + 0.5*sqrt(y^2-r^2)*r - 0.5*y^2*asin((r-1)/y)...
                + 0.5*y^2*asin(r/y) - 0.5*sqrt(y^2-2*y+1-r^2)*r;
        end
        
    end
    
end

Data_working = Data(w+1:2*w,w+1:2*w)';

x = linspace(0,w,w);

F = zeros(w,w);

G = F;

for r = 0:w-1

    F(:,r+1) = exp(-((abs(x)-r)/sqrt(2)).^2);

    G(:,r+1) = A*F(:,r+1);

end

G_inv = pinv(G,tol); % invert the basis matrix to a given tolerance

toc

%% Reconstruct image

tic

%{

This section is the inversion step. For multiple images of the same size,
only this section needs to be repeated. Depending on the symmetry of the
image, either one half of the image or a quarter of the image need be
calculated.

%}

C = zeros(w,w);

Data_Abel_BASEX = zeros(size(Data_working));

for n = 1:w
    
    C(:,n) = G_inv*Data_working(:,n);
    
    Data_Abel_BASEX(:,n) = F*C(:,n);
    
end

Data_Abel_BASEX = [flipud(Data_Abel_BASEX); Data_Abel_BASEX];

if sym == '4-fold'
    
    Data_Abel_BASEX = [fliplr(Data_Abel_BASEX), Data_Abel_BASEX];
    
elseif sym == '2-fold'
    
    Data_working = Data(1:w,w+1:2*w)';
    
    C = zeros(w,w);

    Data_Abel_BASEX_2 = zeros(size(Data_working));

    for n = 1:w
    
        C(:,n) = G_inv*Data_working(:,n);
    
        Data_Abel_BASEX_2(:,n) = F*C(:,n);
    
    end
    
    Data_Abel_BASEX_2 = [flipud(Data_Abel_BASEX_2); Data_Abel_BASEX_2];
    
    Data_Abel_BASEX = [(Data_Abel_BASEX_2), Data_Abel_BASEX];
    
else
    
    disp('Invalid symmetry input: 2-fold and 4-fold only')
    
    Data_Abel_BASEX = zeros(size(Data));
    
    return
    
end

Data_Abel_BASEX = Data_Abel_BASEX';
    
end