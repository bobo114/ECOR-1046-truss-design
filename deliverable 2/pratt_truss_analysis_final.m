n = 13; % number of joists
p = 118.7E3; %pf in newtons

l = 40; % length in m
w = l/(n+1); % joist spacing

height = 3.75; % height in m


R_A = (n/2)*p; % reaction at joint A

loop_length = length(height);
masses = zeros(loop_length, 1);
for i = 1:1:loop_length
    h = height(i);
    s = sqrt(h^2 + w^2);
    
    
    AB = (s/h)*R_A;
    AC = -(w/s)*AB;
    
    data_list = cell((2+(4*(n-1)/2)+1),5);
    data_list(1:2,1:3) = [ % name, force, length
        {"AB", AB, s};
        {"AC", AC, w};
        ];
    
    for j = 1:1:(n-1)/2
        syms A B C;

        Fy = -j*p + R_A + (h/s)*A == 0;
        Fx = (w/s)*A + B + C == 0;
        m = -(j*w)*R_A + p*w*sum(0:(j-1)) - h*B == 0;
        
        sol = solve([Fy, Fx, m], [A, B, C]);
        
        A = double(sol.A);
        B = double(sol.B);
        C = double(sol.C);
        
        if( j == 1 )
            D = -p;
        else
            D = -(h/s)*A;
        end
        A_name = string(char(64+2*j))+string(char(67+2*j));
        B_name = string(char(65+2*j))+string(char(67+2*j));
        C_name = string(char(64+2*j))+string(char(66+2*j));
        D_name = string(char(64+2*j))+string(char(65+2*j));
        
        data_list((3+(j-1)*4):(6+(j-1)*4),1:3) = [
            {A_name, A, s};
            {B_name, B, w};
            {C_name, C, w};
            {D_name, D, h}
            ];
    end
    
    data_list((7+(j-1)*4),1:3) = {"NO", 0, h};
    
    mass = 0;
    for j = 1:length(data_list)
        force = cell2mat(data_list(j, 2));
        len = cell2mat(data_list(j, 3));
        [member_mass, type, max_support] = calculate_mass(force, len);
        data_list(j, 4) = {type};
        data_list(j, 5) = {max_support};
        
        if(j < length(data_list))
            member_mass = member_mass*2;
        end
        
        mass = mass + member_mass;
    end
   
    masses(i) = mass;
end

%display analysis result
[min_mass, ind] = min(masses);
disp("best height:"+string(height(ind))+"m"+", mass at height:"+string(min_mass)+"kN")

% plot graph
figure
plot(height, masses)
xlabel('Height(m)')
ylabel('Weight(kN)')
title('Steel Used in Pratt Truss')

