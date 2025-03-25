function [metric_min, metric_max, lambda_min, lambda_max, fig_handle] = Metric_Graph(lambda, metric, title_string, param_string, metric_string, fig_index)
% Metric_Graph realizes the graph of a metric with lambda (parameters) as x-axis.
% It computes then the points (lambda, metric(lambda))
% associated with the minimum and the maximum of the metric.
% N.B.: This function has to be used when performing the analysis of the metric 
% concerning only one parameter. In the case of two parameters, use Metric_HeatMap

% lambda: vector of lambdas (parameters)
% metric: accounted metric (e.g. RMSE)
% title string: title of the plot
% fig_index: index of the current figure (if no figure, 0)

% metric_min: minimum value of the metric
% metric_max: maximum value of the metric
% lambda_min: lambda value associated with the minimum value of the metric
% lambda_max: lambda value associated with the maximum value of the metric
% fig_handle: handle of the figure

% Minimum 
metric_min = min(metric);
% Maximum 
metric_max = max(metric);
% Identifiying lambda corresponding to the minimum and the maximum of the metric
% Minimum:
lambda_min = lambda(1, metric == metric_min);
% Maximum:
lambda_max = lambda(1, metric == metric_max);

% scatter size
size = 200;

if  fig_index ~=0
    fig_handle = figure(fig_index);
    sct = scatter(lambda, metric, size, 'b', 'filled', 'DisplayName', 'Metric');
    hold on, grid on,
    xlabel(param_string,'Interpreter','latex');
    ylabel(metric_string,'Interpreter','latex');
    title(title_string, 'Interpreter','latex', 'FontWeight', 'bold');
    set(gca, 'XScale', 'log')
% Plotting the minimum and the maximum, treating the multiple cases for the legend
for i = 1 : length(lambda_min)
    if i > 1 
        scatter(lambda_min(i), metric_min, size, 'g', 'filled', 'DisplayName', 'Metric minima');
    else
        sct_min = scatter(lambda_min(i), metric_min, size, 'g', 'filled', 'DisplayName', 'Metric minima');
    end
end
%
for i = 1 : length(lambda_max)
    if i > 1 
        scatter(lambda_max(i), metric_max, size, 'r', 'filled', 'DisplayName', 'Metric maxima'),
    else
        sct_max = scatter(lambda_max(i), metric_max, size, 'r', 'filled', 'DisplayName', 'Metric maxima');
    end
end
legend([sct, sct_min, sct_max]),
hold off
else 
    fig_handle = NaN;
end

end

