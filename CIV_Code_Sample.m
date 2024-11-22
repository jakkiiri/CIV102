% ----- ----- ----- ----- ----- ----- CONSTANTS ----- ----- ----- ----- ----- -----
resolution = 0.01;
bridge_length = 1200;
train_support_tolerance = 0.01;
relative_load_positions = [0,176,340,516,680,856];
% ----- ----- ----- ----- ----- ----- INITIALIZE VARS ----- ----- ----- ----- ----- -----
startpt = train_support_tolerance; % first x-val of envelope
endpt = bridge_length-train_support_tolerance; % last x-val of envelope
diagramvec = startpt:resolution:endpt; % all x-vals of envelope
diagramlen = length(diagramvec); % number of x-vals of envelope
SFE = zeros(1,diagramlen); % initial SFE y-vals
BME = zeros(1,diagramlen); % initial BME y-vals
% ----- ----- ----- ----- ----- ----- LOAD CASE 2: BASE CASE ----- ----- ----- ----- ----- -----
train_load_case_2_base = 452;
factor_heavier = 1.3481481481481481481481481481;
freight_load_case_2_base = train_load_case_2_base/(2+factor_heavier);
train_loads = [0,-freight_load_case_2_base/2,-freight_load_case_2_base/2,-freight_load_case_2_base/2,-freight_load_case_2_base/2,-factor_heavier*freight_load_case_2_base/2,-factor_heavier*freight_load_case_2_base/2,0];
% ----- ----- ----- ----- ----- ----- LOAD CASE 1 ----- ----- ----- ----- ----- -----
train_load_case_1 = 400;
freight_load_case_1 = train_load_case_1/3;
train_loads_case_1 = [0,-freight_load_case_1/2,-freight_load_case_1/2,-freight_load_case_1/2,-freight_load_case_1/2,-freight_load_case_1/2,-freight_load_case_1/2,0];
disp("----- ----- ----- ----- ----- -----")
% ----- ----- ----- ----- ----- ----- SHEAR FORCE ENVELOPE ----- ----- ----- ----- ----- -----
figure
[SFE,max_SF_at]= makeSFE(SFE,startpt,resolution,bridge_length,relative_load_positions,train_support_tolerance,train_loads_case_1,train_load_case_1,diagramvec,diagramlen);
pretty_plot_COLOR(diagramvec,SFE,"Shear Force Envelope",'Distance along span (mm)','Maximum Shear Force V (N)','Shear Force Envelope',0,0,0,0,[0,0.5,0],2);
fprintf("MAX SF sustained by bridge is: %f",max(SFE));
disp("----- ----- ----- ----- ----- -----")
% ----- ----- ----- ----- ----- ----- BENDING MOMENT ENVELOPE ----- ----- ----- ----- ----- -----
figure
[BME,max_BM_at] = makeBME(BME,startpt,resolution,bridge_length,relative_load_positions,train_support_tolerance,train_loads_case_1,train_load_case_1,diagramvec,diagramlen);
pretty_plot_COLOR(diagramvec,BME,"Bending Moment Envelope",'Distance along span (mm)','Maximum Bending Moment M (Nmm)','Bending Moment Envelope',0,10000,0,0,[0,0,1],2);
fprintf("MAX BM sustained by bridge is: %f",max(BME));
disp("----- ----- ----- ----- ----- -----")
% ----- ----- ----- ----- ----- ----- SFD/BMD FOR TRAIN @ MAX SF ----- ----- ----- ----- ----- -----
figure
fprintf("MAXIMUM SHEAR FORCE!!! @ %f mm",max_SF_at)
mySFD = SFD(relative_load_positions+max_SF_at,train_loads_case_1,diagramvec,diagramlen,bridge_length,train_load_case_1);
pretty_plot_COLOR(diagramvec,mySFD,"Shear Force Diagram",'Distance along span (mm)','Shear Force V (N)','Shear Force Diagram for MAX SF',10,10,0,0,[0,0.5,0],2);
fprintf("MAX SF is: %f",max(mySFD));
myBMD = BMD(diagramvec,mySFD);
pretty_plot_COLOR(diagramvec,myBMD,"Bending Moment Diagram",'Distance along span (mm)','Bending Moment M (Nmm)','Bending Moment Diagram for MAX SF',0,10000,0,0,[0,0,1],2);
fprintf("MAX BM is: %f",max(myBMD));
disp("----- ----- ----- ----- ----- -----")
% ----- ----- ----- ----- ----- ----- SFD/BMD FOR TRAIN @ MAX BM ----- ----- ----- ----- ----- -----
figure
fprintf("MAXIMUM BENDING MOMENT!!! @ %f mm",max_BM_at);
mySFD = SFD(relative_load_positions+max_BM_at,train_loads_case_1,diagramvec,diagramlen,bridge_length,train_load_case_1);
pretty_plot_COLOR(diagramvec,mySFD,"Shear Force Diagram",'Distance along span (mm)','Shear Force V (N)','Shear Force Diagram for MAX BM',10,10,0,0,[0,0.5,0],2);
fprintf("MAX SF is: %f",max(mySFD));
myBMD = BMD(diagramvec,mySFD);
pretty_plot_COLOR(diagramvec,myBMD,"Bending Moment Diagram",'Distance along span (mm)','Bending Moment M (Nmm)','Bending Moment Diagram for MAX BM',0,10000,0,0,[0,0,1],2);
fprintf("MAX BM is: %f",max(myBMD));
disp("----- ----- ----- ----- ----- -----")
% ----- ----- ----- ----- ----- ----- MAKING THE ENVELOPES ----- ----- ----- ----- ----- -----
function [SFE,dis_at_max_SF] = makeSFE(SFE,startpt,resolution,bridge_length,relative_load_positions,train_support_tolerance,train_loads,train_load,diagramvec,diagramlen)
    dis_at_max_SF = -1;
    max_SF = -1;
    for train_pos = startpt:resolution:bridge_length-relative_load_positions(end)-train_support_tolerance
        loads_at = relative_load_positions+train_pos;
        tempSFD = SFD(loads_at,train_loads,diagramvec,diagramlen,bridge_length,train_load);
        
        for datapt = 1:diagramlen
            if (abs(tempSFD(datapt))>abs(max_SF))
                dis_at_max_SF = train_pos;
                max_SF = abs(tempSFD(datapt));
            end
            if (abs(tempSFD(datapt))>abs(SFE(datapt)))
                SFE(datapt) = tempSFD(datapt);
            end
        end   
    end
end
function [BME,dis_at_max_BM] = makeBME(BME,startpt,resolution,bridge_length,relative_load_positions,train_support_tolerance,train_loads,train_load,diagramvec,diagramlen)
    dis_at_max_BM = -1;
    max_BM = -1;
    for train_pos = startpt:resolution:bridge_length-relative_load_positions(end)-train_support_tolerance
        loads_at = relative_load_positions+train_pos;
        
        tempSFD = SFD(loads_at,train_loads,diagramvec,diagramlen,bridge_length,train_load);
        tempBMD = BMD(diagramvec,tempSFD);
        
        for datapt = 1:diagramlen
            if (abs(tempBMD(datapt))>abs(max_BM))
                dis_at_max_BM = train_pos;
                max_BM = abs(tempBMD(datapt));
            end
            if (abs(tempBMD(datapt))>abs(BME(datapt)))
                BME(datapt) = tempBMD(datapt);
            end
        end   
    end
end
% ----- ----- ----- ----- ----- ----- RXN FORCES AT SUPPORTS ----- ----- ----- ----- ----- -----
function result = right_support_reaction(train_loads,bridge_length,loads_at)
    result = -(sum(loads_at.*train_loads(2:7)))/bridge_length;
end
function result = left_support_reaction(right_support,train_load)
    result = train_load - right_support;
end
% ----- ----- ----- ----- ----- ----- MAKE THE DIAGRAMS ----- ----- ----- ----- ----- -----
function SFD = SFD(loads_at,train_loads,diagramvec,diagramlen,bridge_length,train_load)
    train_loads(end) = right_support_reaction(train_loads,bridge_length,loads_at);
    train_loads(1) = left_support_reaction(train_loads(end),train_load);
    
    loads_at = [0,loads_at];
    loads_at = [loads_at,bridge_length];
    
    curr_shear = train_loads(1);
    cur_load = 2;
    SFD_index = 1;

    SFD = ones(1,diagramlen);
    for x = diagramvec
        SFD(SFD_index) = curr_shear;
        SFD_index=SFD_index+1;
        if(x>loads_at(cur_load))
            curr_shear = curr_shear+train_loads(cur_load);
            cur_load=cur_load+1;
        end
    end
end
function BMD = BMD(diagramvec,SFD)
    BMD = cumtrapz(diagramvec,SFD);
end
% ----- ----- ----- ----- ----- ----- PRETTY PLOTTING ----- ----- ----- ----- ----- -----
function pretty_plot_COLOR(x,y,legend_data,x_axis,y_axis,graph_title,y_down,y_up,x_left,x_right, color,width)
    plot(x,y,'Color',color,'LineStyle', '-','LineWidth', width);

    axis on;
    grid on;

    yline(0);
    xline(0);

    legend(legend_data, 'FontSize', 14)
    xlabel(x_axis, 'FontSize', 15);
    ylabel(y_axis, 'FontSize', 15);
    title(graph_title, 'FontSize', 18);
    
    
    current_ylim = ylim;
    ylim([current_ylim(1) - y_down, current_ylim(2) + y_up]);
    current_xlim = xlim;
    xlim([current_xlim(1) - x_left, current_xlim(2) + x_right]);
end
