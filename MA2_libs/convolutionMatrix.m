function con_matrix = convolutionMatrix(vector)
con_matrix = zeros(length(vector),length(vector));
con_matrix(:,1) = vector;
for i = 2:length(vector)
    con_matrix(:,i) = [zeros(i-1,1); vector(1:end-(i-1))];
end
end

