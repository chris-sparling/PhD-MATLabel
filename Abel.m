function [Image_Abel] = Abel(Image)

%% Matrix Abel Transform Function - Introduction

%{


%}

%% Define Abel transform matrix

%{

This entire section can be prelocated and saved ahead of time for optimal
speed. Images of different size will require a different Abel transform
matrix.

%}

tic

Image_Abel = zeros(size(Image));

Image = Image';

% Write matrix Abel transform

w =((size(Image,2))/2);

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

%% Abel transform image

%{

This section is the inversion step. For multiple images of the same size,
only this section needs to be repeated.

%}

Image_Abel(w:-1:1,:) = 2*A*Image(w:-1:1,:);
Image_Abel(w+1:2*w,:) = 2*A*Image(w+1:2*w,:);

Image_Abel = Image_Abel';

end

