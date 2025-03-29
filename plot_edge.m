function plot_edge(a, b, theta_star, theta_end, T, color_label, inner_outer_label, legend_label)
% plot_edge 根据弧起点、弧长以及厚度绘制对应的中线及内外边界
%
% 输入参数：
%   a, b         - 圆弧初始半径及单位角增量（径向变化）
%   theta_star   - 起始角（弧度）
%   theta_end    - 弧长（角度差，单位：弧度）
%   T            - 厚度
%   color_label  - 绘图颜色
%   inner_outer_label - 布尔值，若为 true 则绘制内外边界
%   legend_label - 图例名称

plot_theta = linspace(theta_star, theta_star + theta_end, 10000);
plot_r = a + b .* (plot_theta - theta_star);
plot_x = plot_r .* cos(plot_theta);
plot_y = plot_r .* sin(plot_theta);

% 计算中线关系的内边界和外边界
r_inner = plot_r - T/2;
r_outer = plot_r + T/2;
x_inner = r_inner .* cos(plot_theta);
y_inner = r_inner .* sin(plot_theta);
x_outer = r_outer .* cos(plot_theta);
y_outer = r_outer .* sin(plot_theta);

% 绘制中线
plot(plot_x, plot_y, 'LineWidth', 2, 'Color', color_label, 'DisplayName', sprintf('%s-mid', legend_label));
if inner_outer_label
    plot(x_inner, y_inner, '--', 'LineWidth', 1, 'Color', color_label, 'DisplayName', sprintf('%s-inner', legend_label));
    plot(x_outer, y_outer, '--', 'LineWidth', 1, 'Color', color_label, 'DisplayName', sprintf('%s-outer', legend_label));
end
end