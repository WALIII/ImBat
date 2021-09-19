% Stationarity and Flight Types w highest P1 structure?

% KLD_mean KLD_mean_all_windows{1}
% sig_flight_types 
for xx =1:size(KLD_mean_all_windows,2)

    norm_proportion_of_P1_structure = sig_flight_types(1:length(KLD_mean_all_windows{xx}));
    KLD_mean_window = KLD_mean_all_windows{xx};
    R = corrcoef(norm_proportion_of_P1_structure,KLD_mean_window);

    % Plot the training set of data (our noisy y values).
    format long g;
    format compact;
    fontSize = 15;
    figure();
    subplot(2, 1, 1);
    plot(norm_proportion_of_P1_structure,KLD_mean_window, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    grid on;
    xlabel('X', 'FontSize', fontSize);
    ylabel('Y', 'FontSize', fontSize);
    title('Linear Fit', 'FontSize', fontSize);

    % Do the regression with polyfit.  Fit a straight line through the noisy y values.
    linearCoefficients = polyfit(norm_proportion_of_P1_structure,KLD_mean_window, 1)
    xFit = linspace(0, 0.5, 39);
    yFit = polyval(linearCoefficients, xFit);

    hold on;
    plot(xFit, yFit, 'b.-', 'MarkerSize', 15, 'LineWidth', 1);
    legend('Training Set', 'Fit', 'Location', 'Northwest');
end

