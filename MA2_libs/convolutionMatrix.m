function con_matrix = convolutionMatrix(vector,Nr)
con_matrix = zeros(size(vector,1),size(vector,1),Nr);
con_matrix(:,1,1:end) = vector(:,1:end);
for i = 2:size(vector,1)
    con_matrix(:,i,1:end) = [zeros(i-1,Nr); vector(1:end-(i-1),1:end)];
end
end

