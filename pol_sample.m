function polar_image = pol_sample(image)

R = [0:0.5*length(image)];

n = floor(pi*(R+1));

polar_image = NaN([length(R), max(n)]);

tic

for i = 1:length(R)-1
    
    r = R(1,i);
    
    theta = linspace(0,pi,n(1,i));
    
    for a = 1:length(theta)
        
        x_p = r*sin(theta(1,a)) + 0.5*length(image);
        y_p = r*cos(theta(1,a)) + 0.5*length(image);
        
        x_c = round(x_p);
        y_c = round(y_p);
        
        x_dif = x_c - x_p;
        y_dif = y_c - y_p;
        
        sign_x = sign(x_dif);
        sign_y = sign(y_dif);
        
        x_dif = 1 - abs(x_dif);
        y_dif = 1 - abs(y_dif);
        
        polar_image(i,a) = image(y_c+1,x_c+1)*x_dif*y_dif + image(y_c+1,x_c+1-sign_x)*(1-x_dif)*y_dif + image(y_c+1-sign_y,x_c+1)*x_dif*(1-y_dif) + image(y_c+1-sign_y,x_c+1-sign_x)*(1-x_dif)*(1-y_dif);
        
    end
    
end

end

