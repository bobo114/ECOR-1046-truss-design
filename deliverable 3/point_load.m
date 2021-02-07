function load = point_load(letter, list)
    %letter is letter of point load, list is data_list, returns load at
    %point in N
    names = list(1:length(list), 1);
    names = cellstr(names);
    weight_index = contains(names, letter); %generate a logical array of whether row contains needed letter
    load = (sum(cell2mat(list( weight_index, 6))))/2;
    load = load*1.25;
    load = load*1000.0; %convert to kN from N
end
    
    
    