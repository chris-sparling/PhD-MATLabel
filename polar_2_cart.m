function [I_cart, Int] = polar_2_cart(I_polar,res)

% Define source (polar) image coordinate system

r_res = size(I_polar,2);
t_res = size(I_polar,1);

r = linspace(0,r_res,r_res);
theta = linspace(-pi,pi,t_res);



[R, T] = meshgrid(r,theta);

% Define destination (Cartesian) image coordinate system

x = linspace(-r_res,r_res,res);

[X_cart, Y_cart] = meshgrid(x,x);

R_cart = sqrt(X_cart.^2 + Y_cart.^2);

T_cart = atan2(Y_cart,X_cart)';

I_cart = interp2(R,T,I_polar,R_cart,T_cart);

I_cart(isnan(I_cart)) = 0;