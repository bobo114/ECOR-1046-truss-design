%% Function to design the member based on internal force
% Picks the truss type based on given steel section data, will always pick
% the lightest possible member, also returns the weight and max force that
% the member will support
function [weight, type, max_support, stress, delta_L, strain] = calculate_weight(force, len) 
% force = N, length = m, mass = kN, type = designation, max_support = N
    load('steel_info.mat', 'steel_properties');
    phi = 0.9;
    sigma_y = 370;
    L = len*1000;
    
    % Use steel properties cell array, column info per variables below:
    DESIGNATION = 1;
    DEAD_LOAD = 2; % kN/m
    AREA = 3; % mm^2
    R = 4; % mm
    E = 2E5; % modulus of elcticity in MPa
    
    
    if(force >= 0) % tension
        min_area = force/(phi*sigma_y);
        
        sp_len = size(steel_properties);
        sp_len = sp_len(1);
        areas = cell2mat(steel_properties(1:sp_len, AREA));
        
        option_indexes = ge(areas, min_area); % logical array of whther true
        
        possible_options = steel_properties(option_indexes, DESIGNATION:R);
        
        sp_options_len = size(possible_options);
        sp_options_len = sp_options_len(1);
        
        %find min mass
        weights = cell2mat(possible_options(1:sp_options_len, DEAD_LOAD));
        [min_weight_per_m, min_mass_index] = min(weights);

        weight = min_weight_per_m*len;
        type = possible_options(min_mass_index, DESIGNATION);
        area = cell2mat(possible_options(min_mass_index, AREA));
        max_support = phi*sigma_y*area;
        stress = force/area;
        delta_L = stress*(L/E);
    else % compression
        force = abs(force);
        n = 1.34;
        
        K = 1.0;
        min_r = K*L/200;% lowest possible value
        
        sp_len = size(steel_properties);
        sp_len = sp_len(1);
        r_vals = cell2mat(steel_properties(1:sp_len, R));
        
        %filter for min r
        option_indexes = ge(r_vals, min_r); % logical array of whether true
        option_indexes = find(option_indexes);
        
        possible_options = steel_properties(option_indexes, DESIGNATION:R);
        
        %filter for max compression
        
        areas = cell2mat(possible_options(1:length(option_indexes), AREA));
        r_vals = cell2mat(possible_options(1:length(option_indexes), R));
        
        sigma_e = ((pi^2*E)/((K*L)^2)).*(r_vals.^2);
        lambda = sqrt(sigma_y./sigma_e);
        f = 1./((1+lambda.^(2*n)).^(1/n));
        
        max_loads = (phi*sigma_y).*f.*areas;
        
        %2nd filter operation:
        option_indexes = ge(max_loads, force); % logical array of whether true
        %option_indexes = find(option_indexes)
        possible_options = possible_options(option_indexes, DESIGNATION:R);
        f = f(option_indexes); % filter for later use in calculating max force
        
        %find min mass
        weights = cell2mat(possible_options(1:size(possible_options, 1), DEAD_LOAD));
        [min_weight_per_m, min_mass_index] = min(weights);

        weight = min_weight_per_m*len;
        type = possible_options(min_mass_index, DESIGNATION);
        area = cell2mat(possible_options(min_mass_index, AREA));
        max_support = -phi*sigma_y*area*f(min_mass_index); 
        stress = -force/area;
        delta_L = stress*(L/E);
    end
    strain = delta_L/L;
end