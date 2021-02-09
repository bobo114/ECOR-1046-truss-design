n = 13; % number of joists
Dj = 1.53E3;
Pf = 118.7E3; %pf in newtons

l = 40; % length in m
w = l/(n+1); % joist spacing

height = 4.25;% height in m



loop_length = length(height);
weights = zeros(1, loop_length);
for i = 1:1:loop_length
    h = height(i);
    s = sqrt(h^2 + w^2);
    
    list_different = true; % loop control, check if another iteration is needed
    first_calculation = true; % to set truss self weight to 0 in calculations
    data_lists = []; % Collection of all truss iteraations
    
    index = 1; % index used in determining previous results
    data_list = cell((2+(4*(n-1)/2)+1),6); % data from the last truss iteration, gets overwritten every execution
    
    while(list_different) % iterations loop
        
        % create self weight variables
        if(first_calculation)
            weight_at_A = 0;
            total_weight = 0;
        else
            weight_at_A = point_load("A", data_list); %sum of self weight from AC and AB
            
            %convert to kN from N
            total_weight = total_weight*1000; 
            
            %times 1.25
             total_weight = 1.25*total_weight;
        end
        
        Reaction_At_A_from_roof = (n/2)*Pf; % roof weight /2
        Reaction_At_A_from_OWSJ = (n/2)*Dj; % OWSJ weight /2
        
        R_A = Reaction_At_A_from_roof + Reaction_At_A_from_OWSJ + total_weight/2 - weight_at_A ; % sum of forces at A
        point_loads = Dj + Pf; % This is without the self weight

        member_loads = zeros(2, (n-1)/2); % blank member loads array filled with 0s, with indexes below:
        TOP = 1;
        BOTTOM = 2;
            
        % full truss calculation loop
        for j = 1:1:(n-1)/2

            if(first_calculation)   %make top/bottom weights 0 first time
                top_weight = 0;
                bottom_weight = 0;
            else                    % calculate each point load from truss self weight, point A is already accounted for in R_A
                top_point = string(char('A'+2*j));
                top_weight = point_load(top_point, data_list); 
                
                bottom_point = string(char('B'+2*(j-1)));
                bottom_weight = point_load(bottom_point, data_list);
            end
            
            % add loads to list
            member_loads(TOP, j) = top_weight;
            member_loads(BOTTOM, j) = bottom_weight;
            
            % Self weight of section
            self_weight = sum(member_loads, 'all');
            
            distances = ((j-1):-1:(j-6)).*w; % creates list of distances, all negative numbers will be multiplied by 0
            moments_from_self_weight = sum( distances.*member_loads, 'all' ); % sum moments created by member lads
            
            %calculate ends of truss:
            if(j==1)
                AB = double((s/h)*(R_A));
                AC = double(-(w/s)*AB);
                
                data_list(1:2,1:3) = [ % name, force, length
                    {"AB", AB, s};
                    {"AC", AC, w};
                    ];
            end
            
            % method of sections solve for rest of truss:
            syms A B C;
            if(mod(j, 2) == 1) % If odd section
%                 -j*p + R_A  - self_weight;
%                 A=vpa(A);
                Fy = -j*point_loads + R_A + (h/s)*A - self_weight == 0.0;
                m = -(j*w)*R_A + point_loads*w*sum(0:(j-1)) - h*B + moments_from_self_weight == 0.0;

                D = -point_loads - member_loads(TOP, j);

                A_name = string(char(64+2*j))+string(char(67+2*j));
            else % If even section
                Fy = -j*point_loads + R_A - (h/s)*A - self_weight == 0;
                m = -(j*w)*R_A + point_loads*w*sum(0:(j-1)) + h*C + moments_from_self_weight== 0.0;

                D = member_loads(BOTTOM, j);

                A_name = string(char(65+2*j))+string(char(66+2*j));

            end
            %sum forces in x equation and solve forces
            Fx = (w/s)*A + B + C == 0;
            [A, B, C] = solve([Fy, Fx, m]);
            A = double(A);
            B = double(B);
            C = double(C);

            %create forces labels
            B_name = string(char(65+2*j))+string(char(67+2*j));
            C_name = string(char(64+2*j))+string(char(66+2*j));
            D_name = string(char(64+2*j))+string(char(65+2*j));

            %add forces to data list
            data_list((3+(j-1)*4):(6+(j-1)*4),1:3) = [
                {A_name, A, s};
                {B_name, B, w};
                {C_name, C, w};
                {D_name, D, h}
                ];
        end
        
        % Calculate load at NO, accounting for added O
        if(first_calculation)
            o_load = 0;
        else
            o_load = point_load("O",data_list);
        end
        data_list((7+(j-1)*4),1:3) = {"NO", -point_loads-o_load, h};
        
        
        total_weight = 0; % reset total weigjht
        
        % sum weight and add chosen type, its max support and its weight to
        % the data_list
        for j = 1:length(data_list)
            force = cell2mat(data_list(j, 2));
            len = cell2mat(data_list(j, 3));
            [member_mass, type, max_support] = calculate_weight(force, len);
            data_list(j, 4) = {type};
            data_list(j, 5) = {max_support};
            data_list(j, 6) = {double(member_mass)};

            if(j < length(data_list)) % multiply all except last one by 2 due to symmetry
                member_mass = member_mass*2;
            end

            total_weight = total_weight + member_mass; % add to total mass
        end

        
        data_lists = [data_lists, data_list]; % add this last data list to all data lists
        
        %check if any change of members from previous iteration:
        if(index > 1)
            member_designation_column = 4 + (index-2)*6; % column to use in datalists that represents the last iteration
            last_members = data_lists(1:27, member_designation_column); % Last members used list
            
            
            if(isequal(last_members, data_list(1:27, 4))) %compare last list to current list
                list_different = false; % stop iterating
                weights(i) = total_weight; % add to weights list
            end
        end
        
        % Add 1 to index and set first calculation to false
        index = index+1;
        first_calculation = false;
    end
end

% display analysis result
[min_mass, ind] = min(weights);
disp("best height:"+string(height(ind))+"m"+", mass at height:"+string(min_mass)+"kN")

% Create graph
figure
plot(height, weights)
xlabel('Height(m)')
ylabel('Weight(kN)')
title('Steel Used in Warren Truss')

