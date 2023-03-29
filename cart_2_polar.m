function I_polar = cart_2_polar(I_cart,res)

%%

% Takes input I_cart and bins it into a polar coordinate system using a
% natural interpolation method. The grid size and resolution can be
% chosen for the codes application.

%%

x_res = size(I_cart,2);

y_res = size(I_cart,1);

x = linspace(-0.5*x_res,0.5*x_res,x_res);

y = linspace(-0.5*y_res,0.5*y_res,y_res);

[X_cart, Y_cart] = meshgrid(x,y);

%%

r = linspace(0,0.5*(x_res),res(1,1));
theta = linspace(0,2*pi,res(1,2));

[R,T] = meshgrid(r, theta);

X_pol = R.*cos(T);
Y_pol = R.*sin(T);

I_polar = interp2(X_cart,Y_cart,I_cart,X_pol,Y_pol);
