%% 3 C
function Z = calc_z(X, Z0, theta)
    T = length(X);
    Z = zeros(T, 1);
    Z_prev = Z0;
    for t = 1:T
        Z(t) = X(t) - theta * Z_prev;
        Z_prev = Z(t);
    end
end


function l = loglik(X, theta)
    Z = calc_z(X, 0, theta);      
    SSE = sum(Z .^ 2);           
    T = length(X);
    l = -log(SSE / T);             
end

%% 4 A

% Load the data
load('data.mat');  % Assumes variable is called 'data' or a table with named columns

% Subsample: use data up to August 2022
dataS = data(data.Date <= datetime(2022, 8, 31), :);

% Extract the Spread series
Spread = dataS.Spread;
T = length(Spread);

AIC = NaN(5,5);
BIC = NaN(5,5);

for p = 0:5
    for q = 0:5
        try
            model = arima(p, 0, q);
            [fit,~,logL] = estimate(model, Spread);
            
            numParams = p + q + 2;
            [aic, bic] = aicbic(logL, numParams, T);
            
            AIC(p,q+1) = aic;
            BIC(p+1,q+1) = bic;
        catch
            continue;
        end
    end
end

disp('AIC Matrix:');
disp(AIC);
disp('BIC Matrix:');
disp(BIC);


%%
model = arima(0, 0, 4);
[fit, ~, ~, info] = estimate(model, Spread, 'Display', 'off');

theta = cell2mat(fit.MA);  % FIXED
phi = 0;

h = 20; 
theoretical_acf = acf(phi, theta, h);
theoretical_pacf = pacf(phi, theta, h);

figure;
subplot(2,1,1);
stem(1:h, theoretical_acf, 'filled');
xlabel('Lag');
ylabel('ACF');
title('Theoretical ACF of Fitted ARMA(0,4)');

subplot(2,1,2);
stem(1:h, theoretical_pacf, 'filled');
xlabel('Lag');
ylabel('PACF');
title('Theoretical PACF of Fitted ARMA(0,4)');

disp(fit);


disp('Estimated coefficients:');
disp(fit);

residuals = infer(fit, Spread);
figure;
subplot(2,1,1);
plot(residuals);
title('Residuals of ARMA(0,4)');

subplot(2,1,2);
autocorr(residuals);
title('ACF of Residuals');
[h_lb, pValue_lb] = lbqtest(residuals, 'Lags', [5, 10, 15]);

fprintf('Ljung-Box test results:\n');
fprintf('Lag 5: h = %d, p-value = %.4f\n', h_lb(1), pValue_lb(1));
fprintf('Lag 10: h = %d, p-value = %.4f\n', h_lb(2), pValue_lb(2));
fprintf('Lag 15: h = %d, p-value = %.4f\n', h_lb(3), pValue_lb(3));

h = 7;
[YF, YMSE] = forecast(fit, h, 'Y0', Spread);

futureData = data(data.Date > datetime(2022, 8, 31), :);
actuals = futureData.Spread(1:h); 

if length(actuals) >= h
    MSPE = mean((YF - actuals).^2);
    fprintf('MSPE over 7-step forecast: %.4f\n', MSPE);
else
    warning('Not enough actual observations to compute MSPE.');
end

upper = YF + 1.96 * sqrt(YMSE);
lower = YF - 1.96 * sqrt(YMSE);

forecastDates = futureData.Date(1:h);

figure;
plot(forecastDates, YF, 'b-', 'LineWidth', 1.5); hold on;
plot(forecastDates, upper, 'r--');
plot(forecastDates, lower, 'r--');
if length(actuals) >= h
    plot(forecastDates, actuals, 'k-o'); 
end
legend('Forecast', '95% CI', '', 'Actual');
xlabel('Date');
ylabel('Spread');
title('Forecasts with 95% Confidence Interval');
grid on;
