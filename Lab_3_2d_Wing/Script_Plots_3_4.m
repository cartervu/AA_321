%% Plotting Scripts 3, 4, 5, 6
clear all; close all; clc;
 
% - Create function that gets a unit vector along the average slope of line between three points, and then a unit vector normal to this (normal to the surface).
% - get "average" pressure vectors along each given surface
% - Then decompose the vector into x and y components (this gives A and N components).
% - Then add all the A and N components together in a for loop
% - Then translate to L and D
% - Then get coefficients


section = 'E'
absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
% Abs changes for each user that pulls from the github
cd(absp)
load(strcat("2D_Wing_Data_Aggregate",section,".mat"),'Aggregate')


load("portlocations.mat")
port_loc = port_loc * c;
port_loc = [1:44; port_loc(:,:)];
port_loc(:,44) = []

qfolder = ["q20"; "q25"; "q30"; "q35"; "q40"; "q45"];
alpha = ["-4", "-2", "+0", "+2", "+4", "+6", "+8", "+10", "+12", "+14", "+16", "+18", "+20"];

% Initialize the arrays that will contain the normal and axial forces
CLArray = zeros(length(qfolder),length(alpha));
CDArray = zeros(length(qfolder),length(alpha));
xcpArray = zeros(length(qfolder),length(alpha));
Mc4Array = zeros(length(qfolder),length(alpha));


% list of unit normal vectors
vec = [];
for port = 1:43
    vec(port,:) = [port, getnormal(port)]; % [port, xcomp, ycomp]
%     disp('----')
%     disp(port)
%     disp(vec)
end
disp(vec)


% list of distances
distances = [];
distances(:,1) = vec(:,1); % 1st column is port number
for sidenum = 1:43
    distances(sidenum,2) = dist(sidenum);
end


% list of appropriately scaled normal vectors
for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    for aindex = 1:length(alpha)
        % every q/alpha combination has a set of pressure readings
        % we neglect viscous effects
        a = alpha(aindex);
        disp(a)
        Dataset = Aggregate{qindex,aindex};
        if size(Dataset) ~= [0 0]

            % iterate over each port to get the pressure and pressure
            % uncertainty
            % create an array with the port number, uncertainty, and scaled normal vector
            pressures = vec;
            pressures(:,2) = pressures(:,2) .* Dataset(:,2); % x component
            pressures(:,3) = pressures(:,3) .* Dataset(:,2); % y component
            pressures(:,4) = pressures(:,2) .* Dataset(:,3); % uncert x
            pressures(:,5) = pressures(:,3) .* Dataset(:,3); % uncert y
            pressures(end+1,:) = pressures(1); % make last port (44) a repeat
            % of port 1, so that the list of ports comes full circle
            
            % check on uncertainty
%             if max(pressures(:,4)) > 500
%                    disp(pressures)
%             end


            % average the normal vectors to get a vector for each line
            % segment instead of for each point
            lines = [];
            lines(:,1) = pressures(:,1); % port number / line number
            lines(44,:) = [];
            for segment = 1:43
                % average the vectors to get a vector "in between" that's
                % normal to the surface of interest
                lines(segment,2) = (pressures(segment,2)+pressures(segment+1,2))/2;
                lines(segment,3) = (pressures(segment,3)+pressures(segment+1,3))/2;
                lines(segment,4) = sqrt(pressures(segment,4)^2+pressures(segment+1,4)^2);
                lines(segment,5) = sqrt(pressures(segment,5)^2+pressures(segment+1,5)^2);
            end
            % This gives a list of averaged normal vectors for each line 
            % segment instead of each point

            % Multiply each vector by the length of the appropriate
            % line segment
            forceperspan = [];
            forceperspan(:,1) = lines(:,1);
            forceperspan(:,2) = lines(:,2) .* distances(:,2); %axial force
            forceperspan(:,3) = lines(:,3) .* distances(:,2); % normal force
            forceperspan(:,4) = lines(:,4) .* distances(:,2); %axial force uncert
            forceperspan(:,5) = lines(:,5) .* distances(:,2); % normal force uncert
            
            % Take sum of the normal force from each component
            % Take sum of the axial force from each component
            
            A = sum(forceperspan(:,2));
            N = sum(forceperspan(:,3));
%             disp(A)
%             disp(N)
            Auncert = sum(forceperspan(:,4));
            Nuncert = sum(forceperspan(:,5));
            

            
            M1 = - forceperspan(:,2) .* transpose(port_loc(3,:)); % axial force * y distance
            M2 = - forceperspan(:,3) .* transpose(port_loc(2,:)); % normal force * x distance
            M_LE = sum(M1 + M2); %moment about leading edge per unit span
            
            % Note that a positive pitching moment tends to decrease the
            % angle of attack
            
            % unstable at low angles of attack, can get some wild values
            % for x_CP (why?) but it tends to be constant for higher values
            
            % Dataset D looks the best: high xcp for low angles of attack,
            % low xcp for high angles of attack.
            % since the airfoil is symmetrical, we should see a symmetrical
            % x_cp distribution about alpha = 0, we do see this a little
            % bit. Perhaps this is feasible if statistical uncertainties
            % are taken into account.
            
            x_CP = - M_LE / N;
            
            M_c4 = M_LE - c/4*N; 
            cM_c4 = M_c4/(avg(Dataset(:,4)) * c^2);  % moment / q*c^2
            
            xcpArray(qindex,aindex) = x_CP;
            Mc4Array(qindex,aindex) = cM_c4;
            
            
            a = str2double(a)*2*pi/360;
            L = N*cos(a) - A*sin(a);
            D = N*sin(a) + A*cos(a);
            Luncert = sqrt((Nuncert*cos(a))^2 + (Auncert*sin(a))^2);
            Duncert = sqrt((Nuncert*sin(a))^2 + (Auncert*cos(a))^2);
            % Finish the uncertainties later, right now, focus on getting
            % the cl and cd values into an array
            
            CL = L / (avg(Dataset(:,4)) * c); % per unit span
            CD = D / (avg(Dataset(:,4)) * c); % per unit span
            
            CLArray(qindex,aindex) = CL;
            CDArray(qindex,aindex) = CD;
        end
    end
end




if section == "E"
%     CLArray( ~any(CLArray,2), :) = []; % rows
    CLArray( :, ~any(CLArray,1)) = []; % columns
    
%     CDArray( ~any(CDArray,2), :) = []; % rows
    CDArray( :, ~any(CDArray,1)) = []; % columns
    
    Mc4Array( :, ~any(Mc4Array,1)) = []; % columns
    xcpArray( :, ~any(xcpArray,1)) = []; % columns
    
    alpha(2) = [];
    alpha(3) = [];
    alpha(4) = [];
end

cd("Figures/Fig_3")
figure()
for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    plot(str2double(alpha),CLArray(qindex,:),'DisplayName',q)
    hold on
end
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %y-axis
line(xL, [0 0],'Color','black');  %x-axis
title('CL vs alpha')
legend('Location','southeast')
legend
set(gca,'FontSize',20)
hold off
% saveas(gcf,strcat('CL_vs_alpha_',section,'.jpg'))

cd("../Fig_4")
figure()
for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    plot(str2double(alpha),CDArray(qindex,:),'DisplayName',q)
    hold on
end
title('CD vs alpha')
legend('Location','southeast')
legend
set(gca,'FontSize',20)
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %y-axis
line(xL, [0 0],'Color','black');  %x-axis
hold off
% saveas(gcf,strcat('CD_vs_alpha_',section,'.jpg'))



cd("../Fig_5")
figure()
for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    plot(str2double(alpha),xcpArray(qindex,:),'DisplayName',q)
    hold on
end
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %y-axis
line(xL, [0 0],'Color','black');  %x-axis
title('xcp vs alpha')
legend('Location','southeast')
legend
set(gca,'FontSize',20)
hold off
% saveas(gcf,strcat('xcp_vs_alpha_',section,'.jpg'))


cd("../Fig_6")
figure()
for qindex = 1:length(qfolder)
    q = qfolder(qindex);
    plot(str2double(alpha),Mc4Array(qindex,:),'DisplayName',q)
    hold on
end
xL = xlim;
yL = ylim;
line([0 0], yL,'Color','black');  %y-axis
line(xL, [0 0],'Color','black');  %x-axis
title('Mc4 vs alpha')
legend('Location','southeast')
legend
set(gca,'FontSize',20)
hold off
% saveas(gcf,strcat('CL_vs_alpha_',section,'.jpg'))

%% Get Distances

function d = dist(sidenum)
    % the first side clockwise from the leading edge is side 1, goes up to
    % side 43.
    % portnum = sidenum for the port counterclockwise from the side.
    % portnum = sidenum+1 for the port clockwise from the side.
    absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
    % Abs changes for each user that pulls from the github
    cd(absp)
    load("portlocations.mat")
    port_loc = port_loc * c;
    port_loc = [1:44; port_loc(:,:)];
    
    x1 = port_loc(2,sidenum);
    x2 = port_loc(2,sidenum+1);
    y1 = port_loc(3,sidenum);
    y2 = port_loc(3,sidenum+1);
    
    d = sqrt((x2-x1)^2+(y2-y1)^2);
end

%% GetVectors
function vec = getnormal(port)
    absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
    % Abs changes for each user that pulls from the github
    cd(absp)
    load("portlocations.mat")
    port_loc = port_loc * c;
    port_loc = [1:44; port_loc(:,:)];
    
    if 2 <= port & port <= 21
        vec = normvec3pt_upper(port);
    elseif 24 <= port & port <= 43
        vec = normvec3pt_lower(port);
    elseif port == 22
        x1 = port_loc(2,port-1);
        x2 = port_loc(2,port);
        y1 = port_loc(3,port-1);
        y2 = port_loc(3,port);
        slope = (y2-y1)/(x2-x1);
        perp = -1/slope;
        mag = sqrt(1+perp^2);
        vec = -[1/mag, perp/mag];
    elseif port == 23
        x1 = port_loc(2,port);
        x2 = port_loc(2,port+1);
        y1 = port_loc(3,port);
        y2 = port_loc(3,port+1);
        slope = (y2-y1)/(x2-x1);
        perp = -1/slope;
        mag = sqrt(1+perp^2);
        vec = -[1/mag, perp/mag];
    elseif port == 1
        vec = [1, 0];
    end
end


function vec = normvec3pt_upper(port)
% port are all numbers between 2 and 21, inclusive. 21 - 1 = 20 ports total.
% we load the port matrix to find the x and y positions
    absp = "/Users/carter/Software/MATLAB/321_Matlab_Scripts/Lab_3_2d_Wing/";
    % Abs changes for each user that pulls from the github
    
    cd(absp)
    load("portlocations.mat")
    port_loc = port_loc * c;
    port_loc = [1:44; port_loc(:,:)];
    
    
% returns a unit vector normal to the surface
% slope of each line
% y2-y1/x2-x1
    x1 = port_loc(2,port-1);
    y1 = port_loc(3,port-1);
    x2 = port_loc(2,port);
    y2 = port_loc(3,port);
    x3 = port_loc(2,port+1);
    y3 = port_loc(3,port+1);
    sl1 = (y2-y1)/(x2-x1);
    sl2 = (y3-y2)/(x3-x2);
    slavg = avg([sl1,sl2]);
    slperp = -1/slavg;
    
    mag = sqrt(1+slperp^2);
    if slavg > 0
        vec = [1/mag slperp/mag];
    elseif slavg < 0
        vec = -[1/mag slperp/mag];
    end
    % we now have the slope of a line, and need to determine a normal

    % we now have the slope of a line, and need to determine a vector
    % we default to having the vector point downward, and so for the lower
    % surface, we will have to add a negative sign (see function
    % normvec3pt_lower)
    
    % Vector goes 1 over in x, slperp in y, [1,slperp]

    
    % If slavg is pos, slperp is neg
    % if slavg is neg, slperp is pos
end

function vec = normvec3pt_lower(portnum)
% ports 24 thru 43, inclusive 43-23 = 20
    OGvec = normvec3pt_upper(portnum);
    vec = -[OGvec(1), OGvec(2)];
end

%% Statistics Functions
function average = avg(dataset)
    numvals = length(dataset);
    total = 0;
    for valindex = 1:numvals
       total = total + dataset(valindex);
    end
    average = total/numvals;
end

function uncertainty = uncert_digital(dataset)
    uncertainty = (max(dataset)-min(dataset))/(2*sqrt(3));
end

function statistics = stat_digital(dataset)
    statistics = [avg(dataset) uncert_digital(dataset)];
end

function statistics = stat_analog(dataset)
    statistics = [avg(datset) uncert_digital(dataset)/sqrt(2)];
end
