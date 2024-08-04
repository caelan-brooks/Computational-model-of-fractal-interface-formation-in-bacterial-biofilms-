function num_neighbors_edge()
global bio_matrix neighbors row_num col_num

neighbors = zeros(row_num,col_num);
for i=1:row_num
    for j=1:col_num
     if i>1 && i<row_num && j>1 && j<col_num 
        if bio_matrix(i,j) == 2 || bio_matrix(i,j) == 1
            %up
            if bio_matrix(i+1,j) ~= 0
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('up')
            end 
    
            %down
            if bio_matrix(i-1,j) ~= 0
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('down')
            end 
    
            %left 
            if bio_matrix(i,j-1) ~= 0
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('left')
            end 
    
            %right
            if bio_matrix(i,j+1) ~= 0
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('right')
            end 
        end 
     end 
    end 
end 
         return;
end 