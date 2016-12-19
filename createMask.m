
function sigma=createMask(m,n,cx,cy,R)
cc=1:m; 
rr=(1:n).';
cx1=cx; cy1=cy; R1=R;
f=@(xx,yy) (xx-cx1).^2+(yy-cy1).^2 <=R1^2; 
sigma=bsxfun(f,rr,cc); %Logical map of 2 circles
end