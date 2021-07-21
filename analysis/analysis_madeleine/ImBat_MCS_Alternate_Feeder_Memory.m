%ImBat_Alternate_Day_Flights


% fd creatde in Load_Flight_Sequence_Data has the indexes of the flights
% for each day

for i=1:size(flightPaths.Dates,2)
    if fd(i) == i
        temp_day = c_s_34(i)
end