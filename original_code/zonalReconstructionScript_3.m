% using original matrix
%simple script which shows zonal reconstruction as well as
% the corresponding slopes
%
% copyright by:
% Steffen Mauch, (c) 10/2014
% email: steffen.mauch (at) gmail.com
%
% You can redistribute it and/or modify it under the terms of the GNU 
% General Public License as published by the 
% Free Software Foundation, version 2.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


clc
close all
clear all;

n = 32;             % nb lenses in x/y-direction
f = 4796;           % in units [mu]
pixelSize = 7.4;    % in units [um]

nbPixels = 512;
nbPixelsPerLens = nbPixels/n;
defaultGrid = nbPixelsPerLens/2:nbPixelsPerLens:nbPixels;


xDefault = ones(n,1)*defaultGrid;
yDefault = defaultGrid'*ones(1,n);



%%%%%
[tempX,tempY] = meshgrid(-1:1/n:1-1/n); 
%z = exp(-(tempX.^2+tempY.^2) * 1.25 )*10;
A=randn(1,10);
%z = A(1)*tempY+A(2)*tempX+A(3)*(tempX.^2+tempY.^2)+A(4)*tempX.*tempY+A(5)*(tempX.^2+3*tempY.^2);
z = A(1)*...
   ( A(2)*(1-tempX).^2.*exp(-(tempX.^2) - (tempY+1).^2) ... 
   - A(3)*(tempX/5 - tempX.^3 - tempY.^5).*exp(-tempX.^2-tempY.^2) ... 
   - A(4)*exp(-(tempX+1).^2 - tempY.^2) ...
   );
intens = tempX;
[dx,dy,intens] = createSpots(z, intens, tempX, tempY, n, n );
%%%%%

Sx = dx*pixelSize/f*nbPixelsPerLens*pixelSize;
Sy = dy*pixelSize/f*nbPixelsPerLens*pixelSize;

W = zonalReconstruction_hudgin(Sx, Sy, 1);

%W = intgrad2( Sx, Sy );
subplot(2,1,1);
minima = min(min(W));
maxima = max(max(W));
W = resizem(W,2,'bilinear');
norm_W=(W-minima)/(maxima-minima);
surf(norm_W);
colorbar
title('3D plot')
set(gca,'XTickLabel',[0:nbPixelsPerLens:nbPixels]*7.4);
set(gca,'YTickLabel',[0:nbPixelsPerLens:nbPixels]*7.4);
set(gca,'xtick',[0:1:nbPixelsPerLens]);
set(gca,'ytick',[0:1:nbPixelsPerLens]);
xlabel('x-position in [um]')
ylabel('y-position in [um]')
zlabel('z-height in [um]')

subplot(2,1,2);
minima = min(min(z));
maxima = max(max(z));
norm_z = (-z.'+maxima)/(maxima-minima);
surf(norm_z-norm_W)
colorbar
title('3D plot')
set(gca,'XTickLabel',[0:nbPixelsPerLens:nbPixels]*7.4);
set(gca,'YTickLabel',[0:nbPixelsPerLens:nbPixels]*7.4);
set(gca,'xtick',[0:1:nbPixelsPerLens]);
set(gca,'ytick',[0:1:nbPixelsPerLens]);
xlabel('x-position in [um]')
ylabel('y-position in [um]')
zlabel('z-height in [um]')

% change grid spacing
% set(gca,'xtick',[0:nbPixelsPerLens:nbPixels]);
% set(gca,'ytick',[0:nbPixelsPerLens:nbPixels]);
% 
% set(gcf, 'PaperPosition', [0 0 20 20]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [20 20]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'test1', 'pdf') %Save figure
