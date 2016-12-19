function D = delete_multirows(D)
[size_row size_col] = size(D);
j = [];
for i=1:size_row
    if nnz(D(i,:))<2
        j = [j i];
    end
end
D(j,:) =[];