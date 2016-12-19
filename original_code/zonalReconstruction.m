function W = zonalReconstruction(Sx, Sy, ds) 
% zonal reconstruction based on slope data
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


% code from 
% http://books.google.de/books?id=aCC-IciqKkYC&pg=PA122&lpg=PA122&dq=zonalReconstruction+matlab&source=bl&ots=5sSGwNwCdT&sig=EPGccyhgysxa64jdLtIfQ4TXDMg&hl=de&sa=X&ei=l7KtUq2LCMrUtQaMw4CgBw&ved=0CDMQ6AEwAA#v=onepage&q=zonalReconstruction%20matlab&f=false
[n, n]=size(Sx);
S = [reshape(Sx', 1,n*n) reshape(Sy', 1,n*n)]';
A = getA(n);
[U, DD, V] = svd(A, 0);
DD = pinv(DD);
D = getD(n); 
W=V*DD*U'*D*S;
W=reshape(W', n, n)./ds;

%This function obtains the matrix E for the zonal reconstruction function 
function A= getA(n) 
A=zeros(2*n*(n-1), n*n); 
for i=1:n 
	for j=1:(n-1) 
		A((i-1)*(n-1)+j, (i-1)*n+j)=-1; 
		A((i-1)*(n-1)+j, (i-1)*n+j+1)=1; 
		A((n+i-1)*(n-1)+j, i+(j-1)*n)=-1; 
		A((n+i-1)*(n-1)+j, i+j*n)=1; 
	end
end

%This function obtains the matrix C for zonal reconstruction function 
function D=getD(n) 
D=zeros(2*n*(n-1), 2*n*n); 
for i=1:n 
	for j=1:(n-1) 
		D((i-1)*(n-1)+j, (i-1)*n+j)=0.5;
		D((i-1)*(n-1)+j, (i-1)*n+j+1)=0.5; 
		D((n+i-1)*(n-1)+j, n*(n+j-1)+i)=0.5; 
		D((n+i-1)*(n-1)+j, n*(n+j)+i)=0.5; 
	end 
end 
