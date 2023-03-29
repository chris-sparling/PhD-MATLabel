function [Data_FinA] = C_FinA(Data)

%% FinA Fit - Introduction

%% Generate spherical basis functions

n_rad = length(Data)/2;

% make basis

f = zeros([n_rad n_rad]);

for i = 1:n_rad
    
    for j = 1:i
        
        f(i,j) = i*(acos((j-1)/i) - acos((j)/i));
        
    end
    
end

Data_FinA = zeros(size(Data));

I_working = Data(:,n_rad+1:end);

I_recon = zeros(size(I_working));

for i = n_rad:-1:1
    
    I_recon(:,i) = I_working(:,i)/f(i,i);
    
    slice = I_recon(:,i).*f(i,:);
    
    I_working = I_working - slice;
    
    I_working(I_working<0) = 0;
    
end

Data_FinA(:,n_rad+1:end) = (I_recon);

I_working = Data(:,n_rad+1:end);

I_recon = zeros(size(I_working));

for i = n_rad:-1:1
    
    I_recon(:,i) = I_working(:,i)/f(i,i);
    
    slice = I_recon(:,i).*f(i,:);
    
    I_working = I_working - slice;
    
    I_working(I_working<0) = 0;
    
end

Data_FinA(:,1:n_rad) = fliplr(I_recon);



end

