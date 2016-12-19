
function zeroPos = zeroPos(D)
[n m] = size(D);
D_temp = reshape(D', 1,n*m);
len = n*m;
zeroPos = [];
for i=1:len
    if D_temp[i] == 0
        zeroPos = [zeroPos i];
    end
end
