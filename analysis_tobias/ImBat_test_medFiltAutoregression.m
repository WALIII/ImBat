%function ImBat_test_medFiltAutoregression
Afull = full(results.A); %convert matrix from sparse to full
Amat = reshape(Afull,60,80,[]); %reshape so 4800 is now 60x80 pixels
trace = zeros(2,length(Y(1,1,:)));
trace_M3 = zeros(2,length(Y(1,1,:)));
trace_M5 = zeros(2,length(Y(1,1,:)));
trace_M9 = zeros(2,length(Y(1,1,:)));

options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

for i = 1:2
trace_temp = mean((Amat(:,:,i).*Y),[1,2]); %take mean of each frame within the first cell ID from cnmfe
trace(i,:) = reshape(trace_temp,1,[]);
trace_M3(i,:) = medfilt1(trace(i,:),3);
trace_M5(i,:) = medfilt1(trace(i,:),5);
trace_M9(i,:) = medfilt1(trace(i,:),9);

[c1(i), s1(i), options.b, options.g] = foopsi_oasisAR1(y-options.b, options.pars, options.lambda, ...
                options.optimize_b, options.optimize_pars, [], options.maxIter);

end





