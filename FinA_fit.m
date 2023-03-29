function [Data_FinA, I_r, B] = FinA_fit(Data, n_rad, N)

%% FinA Fit - Introduction


B = zeros([N+1, n_rad, size(Data,3)]);

for n = 1:size(Data,3)

%% Transform data into polar coordinates

n_ang = 361;

I_pol = cart_2_polar(Data(:,:,n)', [n_rad, n_ang]);

I_pol = I_pol(1:181,:);


%% Generate spherical basis functions

r_pol = linspace(0,size(Data,2)/2,n_rad).*ones(size(I_pol));

theta_pol = (pi/180).*[0:180]'.*ones(size(I_pol));

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
    
    [b, I_fit] = dLt1D(I_pol(:,i)',N,'half');
    
    B(:,i,n) = b';
    
    I_recon(:,i) = I_fit'/f(i,i);
    
    for t = 1:181 % for each angular bin
        
        if I_pol(t,i) ~= 0 
        
            for j = i:-1:1

                th_idx = round((180/pi)*acos((j/i)*cos(theta_pol(t,1)))) + 1;
                
                %[val, th_idx] = min(abs(th-theta_pol(:,1)));

                W(i,j) = I_recon(th_idx,i);
                
            end

            f_i = f(i,:).*W(i,:);
    
            % subtract scaled and shaped basis function
    
            I_pol(t,:) = I_pol(t,:) - f_i;
            
        else
            
            % do nothing
        
        end
    
    end
    I_recon(I_recon<0) = 0;
    
    
end

B(:,:,n) = B(:,:,n).*linspace(0,size(Data,2)/2,n_rad).^2;

I_r(:,:,n) = sum(I_recon.*r_pol.^2.*sin(theta_pol));

I_recon = [I_recon; flipud(I_recon)];

Data_FinA(:,:,n) = flipud(polar_2_cart(I_recon,size(Data,1)));

end

end

