function find_perimeter_edge()
global bio_matrix perim row_num col_num neighbors area edge_perimeter

    perim = 0;
    area = 0;
    edge_perimeter = zeros(row_num,col_num);
    
    for i=1:row_num
        for j=1:col_num
            if bio_matrix(i,j) == 2 && neighbors(i,j) ~= 0
                area = area + 1;
                perim = perim + (4 - neighbors(i,j));
            end 
            if neighbors(i,j) < 4 && neighbors(i,j) > 0
                edge_perimeter(i,j) = 1;
            end 
        end 
    end 
    return; 
end 