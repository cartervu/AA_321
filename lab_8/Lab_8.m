%  AA 321 Lab 8 Dynamic Analysis

clear all; close all; clc;

cd '/Users/carter/Software/MATLAB/321_Matlab_Scripts/lab_8/data'


% directory = dir;
% s = size(directory);

% 
% for dindex = 5:s(1)
%     DataMatrix = readmatrix(directory(dindex).name);
%     DataMatrix = DataMatrix(7:end,1:10); %format it nicely
%     DataArray{dindex-4} = DataMatrix; % store it in the cell array
% end


% Read in data
DataArray = cell(2,2);
% each column is a material (column 1 is Al, column 2 is Steel)
% each row is a test (row 1 is no mass, row 2 is with mass)
DataArray{1,1} = readmatrix("aluminum1.csv");
DataArray{2,1} = readmatrix("aluminum_mass.csv");
DataArray{1,2} = readmatrix("steel1.csv");
DataArray{2,2} = readmatrix("steel_wmass.csv");

NameArray = ["aluminum no mass" "steel no mass"; "aluminum mass" "steel mass"];


for colnum = 1:2
    for rownum = 1:2
        figure()
        plot(DataArray{rownum,colnum}(:,2),DataArray{rownum,colnum}(:,1))
        xlabel("Time (sec)")
        ylabel("Voltage (V)")
        grid on
        set(gca,'FontSize',20)
        title(strcat("Oscilloscope Voltage Trace ",NameArray(rownum,colnum)))
        saveas(gcf,strcat("VoltageTrace_",NameArray(rownum,colnum),".jpg"))
        % A very unuseful check that the sampling frequency is constant
%         for val = [2.44531, 1.83789, 2.00488, 1.92188]
%             check = ismembertol(DataArray{1,1}(:,2),val,0.00001);
%             linearIndices = find(check==1)
%         end
    end
end





% Attempt at getting peaks out without manually looking 
% This is very difficult, actually. Perhaps I should try looking at
% something like peak width.

% Note: peak width actually worked surprisingly well. Huh.

Q1Stats = cell(2,2);

Fs = 1/0.0009765625;


for colnum = 1:2
    for rownum = 1:2
        [pks,locs,w] = findpeaks(DataArray{rownum,colnum}(:,1),Fs);
        % this takes the voltage data, assumes a given sampling frequency,
        % and outputs the times and 
        % see https://www.mathworks.com/matlabcentral/answers/519463-error-using-findpeaks-expected-y-to-be-a-vector
        
        
        % I want to remove all the peaks that are really small
        s = size(pks);
        s = s(1);
        index = 1;
        wavg = avg(w);
        
        [M,I] = max(pks);
        pks(1:I) = [];
        locs(1:I) = [];
        w(1:I) = [];
        
        while index < s
            s = size(pks);
            s = s(1);
            if w(index) < wavg
                pks(index) = [];
                locs(index) = [];
                w(index) = [];
            else
                index = index + 1;
            end
        end
        
        s = size(pks);
        s = s(1);
        index = 1;
        
        while index < s
            s = size(pks);
            s = s(1);
            if pks(index) < 2 || pks(index) > 7
                pks(index) = [];
                locs(index) = [];
                w(index) = [];
            else
                index = index +1;
            end
        end
        
%         check = issorted(flip(pks)) % okay that the last one doesn't pass -- fluctuations
        Q1Stats{rownum,colnum} = [pks,locs,w];
    end
end



% For these values, we get all the peaks, determined primarily by peak
% width and with a voltage range (4 to 7 volts) to avoid low-amplitude noise or 
% high-amplitude spikes from affecting the result too much. 


% No voltage cut:
% --------------
% aluminum no mass
% average frequency
%    21.6410    2.5085
% 
% --------------
% aluminum mass
% average frequency
%    17.3581    0.6795
% 
% --------------
% steel no mass
% average frequency
%    22.2903   11.2055
% 
% --------------
% steel mass
% average frequency
%    20.7855   24.2373


% with Voltage cut 4 < V < 7

% --------------
% aluminum no mass
% average frequency
%    21.5754    0.7874
% 
% --------------
% aluminum mass
% average frequency
%    17.3339    0.1694
% 
% --------------
% steel no mass
% average frequency
%    22.1965    0.4167
% 
% --------------
% steel mass
% average frequency
%    20.1884    0.9192

% Note that the reduction in uncertainty is really only due to eliminating
% outlier measurements, since a triangular PDF is used to model the data,
% which could definitely be improved -- choose a poisson, ideally.
% However, this amplitude range only corresponds to about a 1-2 second range, which for a 3-10 second
% decay, isn't bad, but could be improved.
   
% Higher frequencies can be observed, (see the very high frequency
% oscillations), but they are very difficult to model based solely on
% peak-finding. The high-amplitude, low-frequency stuff is really all that
% can be picked out -- really only the first mode can be picked out
% visually


% Q2: logarithmic decrement
% c = ln(a1/a2) / (t2-t1)
reportedcvals = cell(2,2);
for colnum = 1:2
    for rownum = 1:2
        pks = Q1Stats{rownum,colnum}(:,1);
        times = Q1Stats{rownum,colnum}(:,2);
%         widths = Q1Stats{rownum,colnum}(:,3);
        s = size(pks);
        s = s(1);
        cvals = zeros(s-1);
        cvals = cvals(:,1); % column vector
        for peakindex = 1:s-1
            cvals(peakindex) = log(pks(peakindex)/pks(peakindex+1))/(times(peakindex+1)-times(peakindex));
            % Note that the high uncertainties are okay -- these are really
            % only a measurement of how far the most extreme values are
            % apart, and do not use the proper PDF
        end
        % Plot the 2-point modeled damping coefficient as a function of time
        figure()
        plot(times(1:s-1),cvals)
        grid on
        reportedcvals{rownum,colnum} = stat_digital(cvals);
        title("Variation in Damping Coefficient vs. Time")
        subtitle(NameArray(rownum,colnum))
        xlabel("Time (sec)")
        ylabel("Damping Coefficient c")
        set(gca,'FontSize',20)
        saveas(gcf,strcat("dampingcoeff_",NameArray(rownum,colnum),".jpg"))
    end
end



for colnum = 1:2
    for rownum = 1:2
        disp("--------------")
        disp(NameArray(rownum,colnum))
        d = delta(Q1Stats{rownum,colnum}(:,2));
%         disp("average Delta t:")
%         disp(d)
        disp("average frequency & uncertainty")
        disp([1/d(1), d(2)/d(1)^2])
        disp("damping coefficient c-vals")
        disp(reportedcvals{rownum,colnum})
    end
end


% Q3 fast fourier transform (fft)
for colnum = 1:2
    for rownum = 1:2
        data_ft = fft(DataArray{rownum,colnum}(:,1));
        s = size(DataArray{rownum,colnum}(:,1),1);
        f = (Fs)*(0:s/2)/s;
        half = abs(data_ft(1:ceil(s/2)+1));
        
        figure()
        semilogy(f,half)
        axis([0,500,0,max(half)])
        xlabel('Frequency [Hz]')
        ylabel('Amplitude |V|')
        grid on
        title(strcat("Frequency Analysis of ",NameArray(rownum,colnum)))
        set(gca,'FontSize',20)
        saveas(gcf,strcat("fourier_",NameArray(rownum,colnum),".jpg"))
    end
end


C = 39.37; % inches per meter
C = 1/C; % meters per inch
rad = 1/2/pi;

% Al, use E = 70 GPa
% 0.507 tall by 1.005 wide by 25 and 1/8 long (I believe this is inches)
EAl = 70 * 10^9;
hAl = 0.507 * C;
wAl = 1.005 * C;
lAl = 25.125 * C;
IAl = wAl * hAl^3 / 12;
rAl = 2710; % kg/m^3
AAl = wAl * hAl;
constAl = sqrt(EAl * IAl / rAl / AAl)/lAl^2;

w1Al = 3.516 * constAl*rad % 25.96 Hz
w2Al = 22.03 * constAl*rad % 162.66 Hz
w3Al = 61.70 * constAl*rad % 455.55 Hz



% Steel, use E = 200 GPa
ES = 200*10^9;
% 0.5 tall by 0.996 wide by 25 and 3/32 long (I belive this is inches)
hS = 0.5 * C;
wS = 0.996 * C;
lS = 25.09375 * C;
IS = wS * hS^3 /12;
rS = 8050; % kg / m^3
AS = wS * hS;
constS = sqrt(ES * IS / rS / AS)/lS^2;

w1S = 3.516*constS*rad % 25.1710 Hz
w2S = 22.03*constS*rad % 157.7122 Hz
w3S = 61.70*constS*rad % 441.7088 Hz

omegas(70*10^9,2710,25.125*C,1.005*C,0.507*C)
omegas(200*10^9,8050,25.09375*C,0.996*C,0.5*C)

%% Check spread
for colnum = 1:2
    for rownum = 1:2
        disp(min(Q1Stats{rownum,colnum}(:,2)))
        disp(max(Q1Stats{rownum,colnum}(:,2)))
    end
end

%% functions

function deltaav = delta(dataset)
    s = size(dataset);
    s = s(1); % assuming column-based dataset
    deltas = zeros(s-1);
    deltas = deltas(:,1);
    for i = 1:s-1
        deltas(i) = abs(dataset(i)-dataset(i+1));
    end
    deltaav = stat_digital(deltas);
end


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

function vec = omegas(E,r,L,w,h)
    I = w*h^3/12;
    A = w*h;
    const = sqrt(E*I/r/A)/L^2;
    omega1 = 3.516*const/2/pi;
    omega2 = 22.03*const/2/pi;
    omega3 = 61.70*const/2/pi;
    vec = [omega1 omega2 omega3];
end