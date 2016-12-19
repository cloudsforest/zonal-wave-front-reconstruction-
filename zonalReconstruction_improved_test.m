function W = zonalReconstruction_improved_test(Sx, Sy, ds, sig) 
% zonal reconstruction based on slope data
%modified the orignal code according to "Improved zonal wavefront reconstruction algorithm for Hartmann
%type test with arbitrary grid patterns" 



%original copyright by:
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
[m, n]=size(Sx);
zeropos = zeroPos(Sx);
S = [reshape(Sx', 1,n*m) reshape(Sy', 1,n*m)]';
S( ~any(S,2), : ) = [];
[sig1,sig2] = getSig(sig,m,n);
A = getA(m,n);
D = getD(m,n);
[A1,D1] = getWeightedAD(A,D,sig1,sig2,m,n);
[U, DD, V] = svd(A1, 0);
DD = pinv(DD);
 
W=V*DD*U'*D1*S;
W=putZerosBack(W,zeropos);
W=reshape(W', m, n)./ds;

%This function obtains the matrix E for the zonal reconstruction function 


function [sig1,sig2] = getSig(sig,m,n)
sig1 = reshape(sig', 1,n*m);
sig2 = [reshape(sig(:,2:n)',1,n*m-m) reshape(sig(2:m,:)',1,n*m-n)];

%This function obtains the matrix E for the zonal reconstruction function 
function A= getA(m,n) 
A=zeros(2*m*n-n-m, m*n); 
for i=1:m-1 
	for j=1:n 
		A((j-1)*(m-1)+i, (j-1)*m+i)=-1;
		A((j-1)*(m-1)+i, (j-1)*m+i+1)=1; 
	end 
end 
for k=1:m*(n-1) 	 
	A((m-1)*n+k, k)=-1; 
	A((m-1)*n+k, k+m)=1;	
end 
%This function obtains the matrix C for zonal reconstruction function 
function D=getD(m,n) 
D = zeros(2*m*n-n-m, 2*m*n); 
for i=1:m-1 
	for j=1:n 
		D((j-1)*(m-1)+i, (j-1)*m+i)=0.5;
		D((j-1)*(m-1)+i, (j-1)*m+i+1)=0.5; 
	end 
end 
for k=1:m*(n-1) 	 
	D((m-1)*n+k, m*n+k)=0.5; 
	D((m-1)*n+k, m*n+k+m)=0.5;	
end 

function [A1,D1] = getWeightedAD(A,D,sig1,sig2,m,n)
D1=zeros(2*m*n-n-m, 2*m*n); 
A1=zeros(2*m*n-n-m, m*n); 
for i=1:2*m*n-n-m
    for j=1:m*n
        D1(i,j)=D(i,j)*sig2(i)*sig1(j);
        D1(i,j+m*n)=D(i,j+m*n)*sig2(i)*sig1(j);
        A1(i,j)=A(i,j)*sig2(i)*sig1(j);        
    end
end
D1( :, ~any(D1,1) ) = [];  %columns
A1( :, ~any(A1,1) ) = [];  %columns
D1 = delete_multirows(D1);
A1 = delete_multirows(A1);
function zeroPos = zeroPos(D)
[n m] = size(D);
D_temp = reshape(D', 1,n*m);
len = n*m;
zeroPos = [];
for i=1:len
    if D_temp(i) == 0
        zeroPos = [zeroPos i];
    end
end



function D = delete_multirows(D)
[size_row size_col] = size(D);
j = [];
for i=1:size_row
    if nnz(D(i,:))<2
        j = [j i];
    end
end
D(j,:) =[];

function M = putZerosBack(D,zeroPos)
sizeD = size(D');
sizeZeroPos = size(zeroPos);
len = sizeD(2)+sizeZeroPos(2);
t=1;
M = zeros(len,1);
for i =1:len
    if ~ismember(i,zeroPos)
        M(i)=D(t);
        t=t+1;
    end
end


        






        
