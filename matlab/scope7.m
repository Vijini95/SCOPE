% ============================================================
%  SCOPE rule for Exponential-Power Distributions
%  β = 0.5, 1.5, 2.5
%  δβ(x) = x |Fβ*(λ x)|^k
%  Author: Brani Vidakovic (2025)
% ============================================================

clear; clc;

% Grid for x
x = linspace(-3,3,1000);

lambda = 1.0;    % scale inside F*( )
k = 3;           % SCOPE exponent

% Exponential-power CDF (standardized)
F_exp_pow = @(x,beta) arrayfun(@(t) ...
    (t<0) * 0.5 * gammainc(abs(t)^beta, 1/beta, 'upper') + ...
    (t>=0)*(1 - 0.5 * gammainc(abs(t)^beta, 1/beta, 'upper')), x);

% Centered CDF F*(x) = 2F(x) - 1
Fstar = @(x,beta) 2*F_exp_pow(x,beta) - 1;

% Values of beta to plot
betas = [1.5, 2.5, 10];

figure; hold on; %grid on;
%title('SCOPE Rule for Exponential-Power Family', 'FontSize', 14);

colors = lines(length(betas));

for j = 1:length(betas)
    beta = betas(j);
    Fmix = arrayfun(@(t) Fstar(lambda*t, beta), x);
    delta = x .* abs(Fmix).^k;
    plot(x, delta, 'LineWidth', 2.5, 'Color', colors(j,:), ...
        'DisplayName', sprintf('\\beta = %.1f', beta));
end

% Identity line
plot(x, x, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Identity');

xlabel('x', 'FontSize', 12);
ylabel('\delta_\beta(x)', 'FontSize', 12);
legend('Location','SouthEast');

% Save figure
print('SCOPE_EP','-dpng','-r300');
print('SCOPE_EP','-dpdf','-r300');

disp('Saved SCOPE_ExpPower_betas.png and .pdf');
