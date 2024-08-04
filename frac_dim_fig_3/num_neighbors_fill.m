function num_neighbors_fill()
global mfill neighbors row_num col_num

neighbors = zeros(row_num,col_num);
for i=1:row_num
    for j=1:col_num
     if i>1 && i<row_num && j>1 && j<col_num 
        if mfill(i,j) == 1
            %up
            if mfill(i+1,j) == 1 
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('up')
            end 
    
            %down
            if mfill(i-1,j) == 1
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('down')
            end 
    
            %left 
            if mfill(i,j-1) == 1
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('left')
            end 
    
            %right
            if mfill(i,j+1) == 1
                neighbors(i,j) = neighbors(i,j) + 1;
                %disp('right')
            end 
            
%             %upper right
%             if mfill(i+1,j+1) == 1
%                 neighbors(i,j) = neighbors(i,j) + 1;
%                 %disp('right')
%             end 
%             
%             %upper left
%             if mfill(i+1,j-1) == 1
%                 neighbors(i,j) = neighbors(i,j) + 1;
%                 %disp('right')
%             end 
%             
%             %lower left
%             if mfill(i-1,j-1) == 1
%                 neighbors(i,j) = neighbors(i,j) + 1;
%                 %disp('right')
%             end 
%             
%             %lower right
%             if mfill(i-1,j+1) == 1
%                 neighbors(i,j) = neighbors(i,j) + 1;
%                 %disp('right')
%             end 
            
        end 
     end 
    end 
end 
         return;
end 