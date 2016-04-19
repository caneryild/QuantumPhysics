% 
% This script plots the probability density function of the 1s orbital
% of a hydrogen-like atom with charge number Z and Bohr radius a0.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;
Z_a0 = Z/a0;

% The normalization constant for the wave function.
N = 1/sqrt(pi) * (Z/a0) ^ 1.5;

% Compute contours of squared wave function values.
% psi = N * exp(-sigma);

theta = 0:0.05:pi;
phi = 0:0.05:2*pi;
[theta,phi]=meshgrid(theta,phi);

C = 0.0001;
c = sprintf('%g *exp(-2*r)- %g',N^2,C);
r = fzero(inline(c,'r'),10);
 
% convert to cartesian coordinates
x = r*sin(theta) .* cos(phi);
y = r*sin(theta) .* sin(phi);
z = r*cos(theta);

% Generate contour plot
figure(1);
hold off;
contour3(x,y,z,25);
c = sprintf('|\\Psi^2| = %g Level Surface',C);
title(c);

% Generate surface plot
figure(2);
hold off;
surf(x,y,z,ones(size(z)));
title(c);

% Generate the radial probability distribution

% Values for r-coordinate.  There is no phi or theta dependence
% for the 1s orbitals.
r = [0:.01:6];
sigma = Z/a0 * r;

% Compute radial probability density
psi_r = 4 * pi * r.^2 .* (N * exp(-sigma)).^2;

% Plot radial distibution function
figure(3);
hold off;
plot(r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');

%%2s--------------------------
% 
% This script plots the radial probability density function of 
% the 2s orbital of a hydrogen-like atom with charge number Z = 1
% and Bohr radius a0 = 1.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;

% The normalization constant for the wave function.
N = 1/4/sqrt(2*pi) * (Z/a0) ^ 1.5;

% Values for r-coordinate.  There is no phi or theta dependence
% for 2s orbitals.
r = Z/a0 * [0 : .01 : 17];
phi = [0 : 0.01 : 2*pi];
theta = [0 : 0.01 : pi];


% The wave function values.
psi = N * (2 - r) .* exp(-r/2);

% Radial wave function
psi_r = 4 * pi * r.^2 .* psi.^2;

figure(1);
plot(r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');
%%2p--------------------------
% 
% This script plots the probability density function of the 2p_z orbital
% of a hydrogen-like atom with charge number Z and Bohr radius a0.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;
Z_a0 = Z/a0;

% The normalization constant for the wave function.
N = 1/4/sqrt(2*pi) * (Z/a0) ^ 1.5;

% Compute contours of the squared wave function values.
% psi = N*sigma*exp(-sigma/2)*cos(theta) 

theta = 0:0.05:pi;
phi = 0:0.05:2*pi+0.1;

[theta,phi]=meshgrid(theta,phi);

siz = size(theta);
r = zeros(1,siz(2));
C = 0.0001;
for i=1:siz(2)
  c = sprintf('%g*r*r*exp(-r)*cos(%g)^2-%g',N^2,theta(1,i),C);
  r(i) = fzero(inline(c,'r'),10);
end

% eliminate NaNs and negative r values.
index = find(r<0);
r(index) = 0;
index = find(r==NaN);
r(index) = 0;

% convert to cartesian coordinates
x = sin(theta).*cos(phi);
x = x*diag(r);
y = sin(theta).*sin(phi);
y = y*diag(r);
z = cos(theta);
z = z*diag(r);

% Generate contour plot
figure(1);
hold off;
contour3(x,y,z,40);
c = sprintf('|\\Psi|^2 = %g Level Surface',C);
title(c);


% Generate surface plot
figure(2);
hold off;
Y = cos(theta);
surf(x,y,z,Y);
title(c);


% Generate the radial probability distribution

% Values for r-coordinate.  There is no phi or theta dependence
% for the 1s orbitals.
r_r = [0:15/1000:15];
sigma_r = Z/a0 * r_r;
theta_r = [-pi:2*pi/1000:pi];

% Compute radial probability density
psi_r = (4*pi/3) * r_r.^2 .* (N * sigma_r .* exp(-sigma_r/2)).^2;

figure(3);
hold off;
plot(r_r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');
%%3d2----------------------
% 
% This script plots the probability density function of the 2p_z orbital
% of a hydrogen-like atom with charge number Z and Bohr radius a0.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;
Z_a0 = Z/a0;

% The normalization constant for the wave function.
N = 1/4/sqrt(2*pi) * (Z/a0) ^ 1.5;

% Compute contours of the squared wave function values.
% psi = N*sigma*exp(-sigma/2)*cos(theta) 

theta = 0:0.05:pi;
phi = 0:0.05:2*pi+0.1;

[theta,phi]=meshgrid(theta,phi);

siz = size(theta);
r = zeros(1,siz(2));
C = 0.0001;
for i=1:siz(2)
  c = sprintf('%g*r*r*exp(-r)*cos(%g)^2-%g',N^2,theta(1,i),C);
  r(i) = fzero(inline(c,'r'),10);
end

% eliminate NaNs and negative r values.
index = find(r<0);
r(index) = 0;
index = find(r==NaN);
r(index) = 0;

% convert to cartesian coordinates
x = sin(theta).*cos(phi);
x = x*diag(r);
y = sin(theta).*sin(phi);
y = y*diag(r);
z = cos(theta);
z = z*diag(r);

% Generate contour plot
figure(1);
hold off;
contour3(x,y,z,40);
c = sprintf('|\\Psi|^2 = %g Level Surface',C);
title(c);


% Generate surface plot
figure(2);
hold off;
Y = cos(theta);
surf(x,y,z,Y);
title(c);


% Generate the radial probability distribution

% Values for r-coordinate.  There is no phi or theta dependence
% for the 1s orbitals.
r_r = [0:15/1000:15];
sigma_r = Z/a0 * r_r;
theta_r = [-pi:2*pi/1000:pi];

% Compute radial probability density
psi_r = (4*pi/3) * r_r.^2 .* (N * sigma_r .* exp(-sigma_r/2)).^2;

figure(3);
hold off;
plot(r_r,psi_r);
grid on;
title('|\Psi(r)|^2');
xlabel('r/a_0');
ylabel('Probability density');
%%3d xy---------------------
% 
% This script plots the probability density function of the 3d_xy orbital
% of a hydrogen-like atom with charge number Z and Bohr radius a0.
%
% Kevin Chu
% Autumn 2001
%

% The charge number and Bohr radius.
Z = 1;
a0 = 1;
Z_a0 = Z/a0;

% The normalization constant for the wave function.
N = 1/81/sqrt(2*pi) * (Z/a0)^1.5;

% Compute contours of the squared wave function values.
% psi = N*sigma^2*exp(-sigma/3)*sin(theta)*sin(theta)*sin(2*phi)

theta = 0:0.1:pi;
phi = 0:0.05:2*pi;

[theta,phi]=meshgrid(theta,phi);

siz_theta = size(theta);
siz_phi = size(phi);
r = zeros(siz_phi(1),siz_theta(2));
C = 0.00001;
for i=1:siz_phi(1)
  for j=1:siz_theta(2)
    c = sprintf('%g*r^4*exp(-2*r/3)*(%g)^4*(%g)^2-%g',N^2,              sin(theta(1,j)), sin(2*phi(i,1)),C);
    r(i,j) = fzero(inline(c,'r'),20);
  end
end

% eliminate NaNs and negative r values.
index = find(r<0);
r(index)=0;
index = find(r==NaN);
r(index)=0;

% convert to cartesian coordinates
x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);
for i=1:siz_phi(1)
  for j=1:siz_theta(2)
    x(i,j) = r(i,j)*x(i,j);
    y(i,j) = r(i,j)*y(i,j);
    z(i,j) = r(i,j)*z(i,j);
  end
end

% Generate a contour plot
figure(1);
hold off;
contour3(x,y,z,40);
c = sprintf('|\\Psi^2| = %g Level Surface',C);
title(c);


% Generate a surface plot
figure(2);
hold off;
Y = sin(theta).^2 .* sin(2*phi);
surf(x,y,z,Y);
title(c);
