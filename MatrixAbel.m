function [Data_Abel] = MatrixAbel(Data,tol)

%% Matrix Abel Transform Function - Introduction

%{

This fast Abel transform function calculates the inverse of the Abel
transform matrix and applies it directly to the image data.

This method is already implemented in the Townsend group for the analysis of
time-resolved photoelectron imaging data.

See https://doi.org/10.1063/1.4765104 for details.

This method can be viewed as a limiting case of the BASEX method, if
Kronecker delta type basis functions are used without any additional
regularisation.

%}

%% Define Abel transform matrix

%{

This entire section can be prelocated and saved ahead of time for optimal
speed. Images of different size will require a different Abel transform
matrix.

%}

Data_Abel = zeros(size(Data));

Data = Data';

% Write matrix Abel transform

w =((size(Data,2))/2);

A = zeros(w,w);

for r=1:length(A) % For each radius
    
    for y=1:length(A) % For each slice of y
        
        if y==r % If your touching the y axis
           A(r,y)=-0.5*sqrt(-1+2*r)*r+0.5*sqrt(-1+2*r)-0.5*r^2*asin((r-1)/r)+0.25*r^2*pi;
           
        elseif y<r  %If below diagonal do nothing 
            
        else        %If above diagonal
            A(r,y)=(-0.5*y^2+y-0.5)*asin(r/(y-1)) + (0.5*y^2-y+0.5)*asin((r-1)/(y-1))...
                + 0.5*sqrt(y^2-2*y+2*r-r^2)*r - 0.5*sqrt(y^2+2*r-1-r^2)*r - 0.5*sqrt(y^2-2*y+2*r-r^2)...
                + 0.5*sqrt(y^2+2*r-1-r^2) + 0.5*sqrt(y^2-r^2)*r - 0.5*y^2*asin((r-1)/y)...
                + 0.5*y^2*asin(r/y) - 0.5*sqrt(y^2-2*y+1-r^2)*r;
        end
        
    end
    
end

A_inv = pinv(A,tol); % invert A

%% Abel transform image

%{

This section is the inversion step. For multiple images of the same size,
only this section needs to be repeated.

%}

tic

Data_Abel(w:-1:1,:) = 0.5*A_inv*Data(w:-1:1,:); % invert half image
Data_Abel(w+1:2*w,:) = 0.5*A_inv*Data(w+1:2*w,:); % invert other half


Data_Abel = Data_Abel'; % transpose output

end

