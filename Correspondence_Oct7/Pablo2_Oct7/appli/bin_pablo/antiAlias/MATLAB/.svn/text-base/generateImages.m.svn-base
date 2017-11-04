%% initialize
dims = [65 65 65];
img = false(dims(1), dims(2), dims(3));
dims = dims - 1;
x = -dims(1)/2:dims(1)/2;
y = -dims(2)/2:dims(2)/2;
z = -dims(3)/2:dims(3)/2;
dims = dims + 1;
[x,y,z] = ndgrid(x,y,z);

%% a thin slab
img(abs(x) <= 20 & abs(y) <= 20 & -1 <= z & z < 0.5) = true;

%% circle
img( sqrt(x.^2 + y.^2 + z.^2) <= 20) = true;
clear x y z;

%% ellipsoid
img( (x/(dims(1)/3)).^2 + (y/(dims(1)/6)).^2 + (z/(dims(1)/12)).^2 <= 1) = true;
clear x y z;

%% some quartic surface (surface has concavities)

a = 0; b = -7; c = 0;
x = x/8;
y = y/7;
z = z/6;
img( x.^4 + y.^4 + z.^4 + a*(x.^2 + y.^2 + z.^2).^2 + b*(x.^2 + y.^2 + z.^2) + c <= 0) = true;
img( 31, 31, 31 ) = 1;
clear x y z a b c;

%%
for z = 31
	clf;
% 	imagesc(phi(:,:,z)*10, [-10 10]);
	imagesc(img(:,:,z), [0 1]);
	axis image;
	colorbar;
	hold on;
	contour(phi(:,:,z),1e-100,'k');
	pause(0.5)
end
%plot( 21+10*cos(0:pi/50:2*pi),21+10*sin(0:pi/50:2*pi),'w')
