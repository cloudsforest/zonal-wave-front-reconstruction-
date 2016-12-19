function [spot_x spot_y intens] = createSpots_improved(wfs, intensity, x, y, apart_x, apart_y)
% simple function which create spot pattern expected from an SHWFS sensor
%
% wfs: is the z-data of the original wfs data
% intensity: same dimension as wfs; could be gaussian shape
% x,y: the x- and y-position
% apart_x, apart_y: the number of lenses in x-/y-direction (lenslet array)
%
% copyright by:
% Steffen Mauch, (c) 12/2014
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


dim = size(wfs,1);
n=apart_x;
%s1=createMask(n,n,(n-1)/2+1,(n-1)/2+1,(n-1)/2-1);
s=double(createMask(n,n,(n+1)/2,(n+1)/2,(n-1)/2));;
%s = xor(s1,s2);
length_apart_x = floor(dim / apart_x);
length_apart_y = floor(dim / apart_y);

spot_x = zeros(apart_x,apart_x);
spot_y = zeros(apart_y,apart_y);

intens = zeros(apart_x, apart_y);

warning('on','MATLAB:rankDeficientMatrix');

for l=1:apart_y
    for k=1:apart_x

        X = (k-1)*length_apart_x+(1:length_apart_x);
        Y = (l-1)*length_apart_y+(1:length_apart_y);
                
        intens(k,l) = sum( sum(intensity(X,Y),1), 2 ) / (size(X,2)*size(Y,2) );
        temp = wfs(X, Y);

        index = ~isnan(temp);
        temp = temp(index);

        X_value = x(X,Y);  Y_value = y(X,Y);
        X = X_value(index);  Y = Y_value(index);
        X = X(:);  Y = Y(:);

        if( size(~isnan(temp),1) == 0 ) 
            spot_x(k,l) = NaN;
            spot_y(k,l) = NaN;
        else
            %X = [ones(size(x)) x y];
            %a = X\temp * faktor
            
            A = [X,Y,temp,ones(length(X),1)];
            [U,S,V] = svd(A);
            coeff = V(:,4)/V(3,4);
           % if (s(k,l)==1)
                spot_x(k,l) = coeff(1);
                spot_y(k,l) = coeff(2);
           % end
        end
    end
end

spot_x(abs(spot_x) < 1e-8) = 0;
spot_y(abs(spot_y) < 1e-8) = 0;
spot_x=spot_x.*s;
spot_y=spot_y.*s;
% spot_x(abs(spot_x) > 3) = 0;
% spot_y(abs(spot_y) > 3) = 0;
% spot_x(abs(spot_x) < -3) = 0;
% spot_y(abs(spot_y) < -3) = 0;
figure(1)
quiver(1:apart_x,1:apart_y,spot_x,spot_y);

warning('on','MATLAB:rankDeficientMatrix');
