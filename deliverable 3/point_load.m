%% Point load calculator function
% calculates the factored load created at each point by the self weight of
% the truss
function load = point_load(letter, list)
    %letter is letter of point load, list is data_list, returns load at
    %point in N
    names = list(1:length(list), 1);
    names = cellstr(names);
    weight_index = contains(names, letter); %generate a logical array of whether row contains needed letter
    load = (sum(cell2mat(list( weight_index, 6))));
    
     %fix for NO end
    if (eq("O", letter) )
        weight_to_add_index = contains(names, "MO"); %index of NO
        load_to_add = (sum(cell2mat(list( weight_to_add_index, 6)))); % weight to remove
        load = load + load_to_add;
    end
    load = load/2;
    load = load*1.25;
    load = load*1000.0; %convert to kN from N
end
    
    
    