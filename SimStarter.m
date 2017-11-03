clear all; clc; close all;
%Simulation Starter for Maximum Drag Orbital Propagator
%M. McDonald 10/29/17

%Vehicle params
a_max = 0.185806; %m^2
msc = 50.0; %kg.

%ECEF SV:
r_ECEF_i = [-4292653.41; 955168.47; 5139356.57]; %inertial position represented in ECEF frame, m.
v_ECEF_i = [109.649663; -7527.726490; 1484.521489]; %inertial velocity represented in ECEF frame, m/s.
time_vec = [2017 10 29 0 0 0];
Recef_i = dcmeci2ecef('IAU-2000/2006',time_vec);

r_ECI_i = transpose(Recef_i)*r_ECEF_i;
v_ECI_i = transpose(Recef_i)*v_ECEF_i;

initial_conditions = [r_ECI_i; v_ECI_i]';

%sim setup:
dt = 1; %sim discrete timestep, seconds
deci = 86400;
years_duration = 10; 
seconds_year = 31557600;
t_end = years_duration.*seconds_year; %whole number only.

%setup timevec needed for sim:
year = time_vec(1,1);
month = time_vec(1,2);
day = time_vec(1,3);
hour = time_vec(1,4);
minute = time_vec(1,5);
second = time_vec(1,6);

%1. determine if year is a leap year:
if (mod((year-2000),4) == 0)
    if month == 1
        month_doy = 0;
    elseif month == 2
        month_doy = 31;
    elseif month == 3
        month_doy = 31 + 29;
    elseif month == 4
        month_doy = 31 + 29 + 31;
    elseif month == 5
        month_doy = 31 + 29 + 31 + 30;
    elseif month == 6
        month_doy = 31 + 29 + 31 + 30 + 31;
    elseif month == 7
        month_doy = 31 + 29 + 31 + 30 + 31 + 30;
    elseif month == 8
        month_doy = 31 + 29 + 31 + 30 + 31 + 30 +31;
    elseif month == 9
        month_doy = 31 + 29 + 31 + 30 + 31 + 30 +31 + 31;
    elseif month == 10
        month_doy = 31 + 29 + 31 + 30 + 31 + 30 +31 + 31 + 30;
    elseif month == 11
        month_doy = 31 + 29 + 31 + 30 + 31 + 30 +31 + 31 + 30 +31;
    else
        month_doy = 31 + 29 + 31 + 30 + 31 + 30 +31 + 31 + 30 +31 + 30;
    end
else % not a leap year
     if month == 1
        month_doy = 0;
    elseif month == 2
        month_doy = 31;
    elseif month == 3
        month_doy = 31 + 28;
    elseif month == 4
        month_doy = 31 + 28 + 31;
    elseif month == 5
        month_doy = 31 + 28 + 31 + 30;
    elseif month == 6
        month_doy = 31 + 28 + 31 + 30 + 31;
    elseif month == 7
        month_doy = 31 + 28 + 31 + 30 + 31 + 30;
    elseif month == 8
        month_doy = 31 + 28 + 31 + 30 + 31 + 30 +31;
    elseif month == 9
        month_doy = 31 + 28 + 31 + 30 + 31 + 30 +31 + 31;
    elseif month == 10
        month_doy = 31 + 28 + 31 + 30 + 31 + 30 +31 + 31 + 30;
    elseif month == 11
        month_doy = 31 + 28 + 31 + 30 + 31 + 30 +31 + 31 + 30 +31;
    else
        month_doy = 31 + 28 + 31 + 30 + 31 + 30 +31 + 31 + 30 +31 + 30;
     end
end

timevec_o = [year; (day + month_doy); ((hour.*3600) + (minute.*60) + second)];
%output IC data to file for sim to read:
fileID = fopen('IC.txt','w');
A = [a_max msc r_ECI_i' v_ECI_i' timevec_o' dt t_end deci];
fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f',A);
fclose(fileID);

save('time_vec.mat','time_vec');
%Ready to start sim.