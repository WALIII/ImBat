function [] = savevars(KC,pruned_dataset,takeoff_locations_cpu)

save('KC','KC');
save('takeoff_locations_cpu','takeoff_locations_cpu');
save('pruned_dataset','pruned_dataset');
save('outliers','outliers');

end