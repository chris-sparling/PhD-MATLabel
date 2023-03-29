function [Data_FinA, I_r] = FinA(Data, n_rad)

%% FinA Fit - Introduction

%% Transform data into polar coordinates

n_ang = 361;

I_pol = cart_2_polar(Data', [n_rad, n_ang]);

%I_pol = I_pol(1:181,:);

%% Generate spherical basis functions

r_pol = linspace(0,size(Data,2)/2,n_rad).*ones(size(I_pol));

theta_pol = (pi/180)*[0:360]'.*ones(size(I_pol));


% make basis

f = zeros([n_rad n_rad]);

for i = 1:n_rad
    
    for j = 1:i
        
        f(i,j) = i*(acos((j-1)/i) - acos((j)/i));
        
    end
    
end


% Process image

I_recon = zeros(size(I_pol));

W = zeros([n_rad,n_rad]);

for i = n_rad:-1:1 % for each radii
    
    I_recon(:,i) = I_pol(:,i)/f(i,i);
    
    for t = 1:361 % for each angular bin
        
        if I_pol(t,i) ~= 0 
        
            for j = i:-1:1

                th_idx = round((180/pi)*acos((j/i)*cos(theta_pol(t,1)))) + 1;

                W(i,j) = I_recon(th_idx,i);
                
            end

            f_i = f(i,:).*W(i,:);
    
            % subtract scaled and shaped basis function
    
            I_pol(t,:) = I_pol(t,:) - f_i;
            
        else
            
            % do nothing
        
        end
    
    end
    
%I_recon(I_recon<0) = 0;   
   
end



I_r = sum(I_recon.*r_pol.^2.*abs(sin(theta_pol)));

%I_recon = [I_recon; flipud(I_recon)];

Data_FinA = polar_2_cart(I_recon,size(Data,1));

end

