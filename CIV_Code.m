% Initialize constants variables
% distances in mm, forces in N, stress in MPa
step = 1; % change in train loading position in each iteration
% length of bridge; left support: x=0, right support: x=1200
bridge_length = 1200; 
% positions of loading relative to first load
% assume train drives from left to right
loading_pos = [0,176,340,516,680,856]; 

% first location, increase by step each loop
% should start with load 0 at x=0
% positions where not all loads on bridge are ignored
start_pos = 0; 
end_pos = bridge_length-856; 
% end position such that load @856 is at x=1200

list_pos = start_pos:step:end_pos; % all possible positions with entire train on bridge
list_x = 0:1:1200; % x-axis as discrete values
% to store max SF and BM along each point on x=axis
SFE_empty = zeros(1,length(list_x));
BME_empty = zeros(1,length(list_x));

% Load Case 1
% 400N distributed evenly across 6 axels; ~66.7N per axel
load_case1 = ones(6,1)'*400/6;
[SFE1, V_max1, BME1, M_max1]=make_SFE_BME(list_pos, list_x, load_case1,loading_pos, SFE_empty, BME_empty);

% Load Case 2 - Base Case
% 2 scenarios: locomotive first or locomotive last
% base case: locomotive = 182N, freights = 135 x 2
% locomotive first: [fr/2, fr/2, fr/2, fr/2, loco/2, loco/2]
load_case2_loco_first = [135/2, 135/2, 135/2, 135/2, 182/2, 182/2];
load_case2_loco_last = flip(load_case2_loco_first);
[SFE2a, V_max2a, BME2a, M_max2a]=make_SFE_BME(list_pos, list_x, load_case2_loco_first, loading_pos, SFE_empty, BME_empty);
[SFE2b, V_max2b, BME2b, M_max2b]=make_SFE_BME(list_pos, list_x, load_case2_loco_last, loading_pos, SFE_empty, BME_empty);
for pos = list_x+1
    if abs(SFE2a(pos)) > abs(SFE2b(pos))
        SFE2b(pos) = SFE2a(pos);
    end
    if abs(BME2a(pos)) > abs(BME2b(pos))
        BME2b(pos) = BME2a(pos);
    end
end
SFE2=SFE2b;
BME2=BME2b;
V_max2=max(abs(V_max2a),abs(V_max2b));
M_max2=max(abs(M_max2a),abs(M_max2b));

% actually getting the graphs and outputs
% Case 1
subplot(2,1,1)
plot_data(list_x,SFE1,"Shear Force Envelope","Distance Along Bridge (mm)","Max Shear Force (N)","Shear Force Envelope Case 1")
disp("Case 1 - Max Shear:")
disp(V_max1)

subplot(2,1,2)
plot_data(list_x,BME1,"Bending Moment Envelope","Distance Along Bridge (mm)","Max Bending Moment (Nmm)","Bending Moment Envelope Case 1")
disp("Case 1 - Max Bending Moment:")
disp(M_max1)

% Case 2 base case
figure
subplot(2,1,1)
plot_data(list_x,SFE2,"Shear Force Envelope","Distance Along Bridge (mm)","Max Shear Force (N)","Shear Force Envelope Case 2 Base Case")
disp("Case 2 Base Case - Max Shear:")
disp(V_max2)

subplot(2,1,2)
plot_data(list_x,BME2,"Bending Moment Envelope","Distance Along Bridge (mm)","Max Bending Moment (Nmm)","Bending Moment Envelope Case 1")
disp("Case 2 Base Case - Max Bending Moment:")
disp(M_max2)

% exporting data as txt file
data_diagrams = [SFE1;BME1;SFE2;BME2];
data_max = [V_max1;M_max1;V_max2;M_max2];
writematrix(data_diagrams,"SFE_BME.txt");
writematrix(data_max,"Vmax_Mmax.txt");

% Functions

% calulate the reaction forces for any loading
function [left_sup, right_sup] = calc_rxn(cur_pos, loading)
    % cur_pos contains position of loads as x-coords
    % sum moments at x=0 (left support) --> obtains right support
    load_M = 0;
    % multiply position of load (wrt 0) by load to get moment
    for M = cur_pos.*loading
        load_M = load_M + M;
    end
    right_sup = load_M/1200; % since right support is 1200 away
    left_sup = sum(loading) - right_sup;
end

% function to find SFD (at x values from 1-1200, with step 1)
function SFD = calc_SFD(left_sup, right_sup, loading, pos_list, list_x)
    % set entire SFD to left_sup force
    SFD = ones(1,length(list_x)).*left_sup;
    % keeping track of current values
    cur_V = left_sup;
    load_num = 1;
    % iterate through the positions of loading
    for cur_pos = pos_list
        % changes current shear from corresponding value in load list
        cur_V = cur_V-loading(load_num);
        % updates so next load is applied
        load_num = load_num + 1;
        % changes all SFD values after current position to new shear
        SFD(cur_pos+1:length(SFD)) = cur_V; 
        % +1 to account for the zero; x=1200 has position of 1201 in vector
    end
    % finishes SFD with the right support
    SFD(length(SFD)) = SFD(length(SFD))+right_sup;
end

% function to find BMD by integrating SFD
function BMD = calc_BMD(SFD)
    % SFD as parameter
    % integrates SFD using trapezoid approx.
    BMD=cumtrapz(SFD);
end

% function to make SFE and BME
function [SFE, V_max,BME, M_max] = make_SFE_BME(list_pos, list_x, loading, loading_pos, empty_SFE, empty_BME)
    SFE=empty_SFE;
    BME=empty_BME;
    for pos_adjust = list_pos
        cur_pos_list = loading_pos+pos_adjust;
        [L_sup, R_sup] = calc_rxn(cur_pos_list, loading);
        cur_SFD = calc_SFD(L_sup,R_sup, loading, cur_pos_list, list_x);
        cur_BMD = calc_BMD(cur_SFD);
        
        for pos = list_x+1
            if abs(cur_SFD(pos)) > abs(SFE(pos))
                SFE(pos) = cur_SFD(pos);
            end
            if abs(cur_BMD(pos)) > abs(BME(pos))
                BME(pos) = cur_BMD(pos);
            end
        end
    end
    V_max=max(abs(SFE));
    M_max=max(abs(BME));
end

% function to plot the data
function plot_data(x,y,legend,x_lab,y_lab,the_title)
    plot(x,y,"b.-")
    
    axis on;
    grid on;
    xlabel(x_lab);
    ylabel(y_lab);
    title(the_title);

    current_ylim = ylim;
    ylim([current_ylim(1) - 20, current_ylim(2) + 20]);
    current_xlim = xlim;
    xlim([current_xlim(1) - 20, current_xlim(2) + 20]);
    
end