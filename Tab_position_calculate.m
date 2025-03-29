function Tab_position_calculate(options)
% Tab_position_calculate 计算双极耳位置并绘制隔膜、负极、正极及其层间接触的卷绕几何图形
% 用法：
%     options.neg_thickness = 47.5*2/1000;  % 负极双面料厚[mm]
%     options.sep_thickness = 12/1e3;       % 隔膜厚度[mm]
%     options.pos_thickness = 38.5*2/1000;  % 正极双面料厚[mm]
%     options.foil_pos_thickness = 12/1000; % 正极箔材厚度[mm]
%     options.foil_neg_thickness = 8/1000;  % 负极箔材厚度[mm]
%     options.sep_arc_length=16;            % 隔膜入位空卷长度
%     options.L_neg_star=6;                 % 负极入位长度[mm]
%     options.nail_a = 1.75;                % 卷针半径[mm]
%     options.pos_length = 1411;            % 正极滚压后长度[mm]
%     options.neg_length = 1461;            % 负极滚压后长度,含留白宽度[mm]
%     options.w=11;                         % 极耳滚压前留白宽度[mm]
%     options.L_first=340;                  % 负极第一段涂布滚压后长度,约总长1/4[mm]
%     options.L_second=700;                 % 负极第二段涂布滚压后长度,约总长1/2[mm]
%     options.el_rate=1.02;                 % 极片延展率
%     options.calculate_flag='pos'          % 计算正极或者负极标签，正极pos,负极'neg'  
%     Tab_position_calculate(options);
%  极片尺寸示意图：
%  |L_neg_star----------|---|------------|---|-----------|
%  |------L_neg_first---|-w-|L_neg_second|-w-|-----------|
%  |--------------------neg_length-----------------------|
arguments
    options.neg_thickness = 47.5*2/1000;  % 负极双面料厚[mm]
    options.sep_thickness = 12/1e3;       % 隔膜厚度[mm]
    options.pos_thickness = 38.5*2/1000;  % 正极双面料厚[mm]
    options.foil_pos_thickness = 12/1000; % 正极箔材厚度[mm]
    options.foil_neg_thickness = 8/1000;  % 负极箔材厚度[mm]
    options.sep_arc_length=16;            % 隔膜入位空卷长度
    options.L_neg_star=6;                 % 负极入位长度[mm]
    options.nail_a = 1.75;                % 卷针半径[mm]
    options.pos_length = 1411;            % 正极滚压后长度[mm]
    options.neg_length = 1461;            % 负极滚压后长度,含留白宽度[mm]
    options.w=11;                         % 极耳滚压前留白宽度[mm]
    options.L_first=340;                  % 负极第一段涂布滚压后长度,约总长1/4[mm]
    options.L_second=700;                 % 负极第二段涂布滚压后长度,约总长1/2[mm]
    options.el_rate=1.02;                 % 极片延展率
    options.calculate_flag='pos'          % 计算正极或者负极标签，正极pos,负极'neg'  
end
%% 1. 参数提取
neg_thickness      = options.neg_thickness;
sep_thickness      = options.sep_thickness;
pos_thickness      = options.pos_thickness;
foil_pos_thickness = options.foil_pos_thickness;
foil_neg_thickness = options.foil_neg_thickness;
sep_arc_length     =options.sep_arc_length;
pos_length         = options.pos_length;
neg_length         = options.neg_length;
nail_a             = options.nail_a;
w                  =options.w* options.el_rate;
L_neg_star         =options.L_neg_star;
L_first            =options.L_first;
L_second           =options.L_second;
el_rate            =options.el_rate;
calculate_flag     =options.calculate_flag;

%% 2. 隔膜中线及厚度边界计算
% 单位角增量对应的径向厚度变化（隔膜空卷）
sep_b = (sep_thickness*2) / (2*pi);

% 调用 solve_theta 得到隔膜弧对应的极角及半径
[sep_theta, sep_r] = solve_theta(nail_a + sep_thickness, sep_b, 0, sep_arc_length);

% 隔膜中线坐标（极坐标）
plot_theta_list_sep = linspace(0, sep_theta, 2000);
plot_r_sep   = nail_a + sep_b .* plot_theta_list_sep;
plot_x_sep   = plot_r_sep .* cos(plot_theta_list_sep);
plot_y_sep   = plot_r_sep .* sin(plot_theta_list_sep);

% 隔膜厚度（内、外边界）——简单进行径向内缩与外延扩展
T_sep = sep_thickness * 2;
plot_r_sep_inner = plot_r_sep - T_sep/2;
plot_r_sep_outer = plot_r_sep + T_sep/2;

plot_x_sep_inner = plot_r_sep_inner .* cos(plot_theta_list_sep);
plot_y_sep_inner = plot_r_sep_inner .* sin(plot_theta_list_sep);
plot_x_sep_outer = plot_r_sep_outer .* cos(plot_theta_list_sep);
plot_y_sep_outer = plot_r_sep_outer .* sin(plot_theta_list_sep);

%% 3. 绘制隔膜空卷
figure;
plot(plot_x_sep, plot_y_sep, 'LineWidth', 1.5, 'Color', 'b', 'DisplayName', 'sep-star-mid');
hold on;
plot(plot_x_sep_inner, plot_y_sep_inner, '--', 'LineWidth', 1, 'Color', 'b', 'DisplayName', 'sep-star-inner');
plot(plot_x_sep_outer, plot_y_sep_outer, '--', 'LineWidth', 1, 'Color', 'b', 'DisplayName', 'sep-star-outer');

%% 4. 负极几何计算
% 负极入位参数：其起始半径依赖于隔膜末端半径以及负极、箔材的厚度中心位置
neg_a = sep_r - sep_thickness + (neg_thickness + foil_neg_thickness) / 2;
neg_b = (neg_thickness + foil_neg_thickness) / (2*pi);
neg_arc_length = L_neg_star; % 负极入位长度[mm]
[neg_theta, neg_r] = solve_theta(neg_a, neg_b, 0, neg_arc_length);

neg_total_thickness = neg_thickness + foil_neg_thickness;
% 负极入位中线角度从 sep_theta 开始
plot_theta_list_neg = linspace(sep_theta, sep_theta + neg_theta, 2000);
plot_r_neg = neg_a + neg_b .* (plot_theta_list_neg - sep_theta);
plot_x_neg = plot_r_neg .* cos(plot_theta_list_neg);
plot_y_neg = plot_r_neg .* sin(plot_theta_list_neg);

% 负极厚度处理：分别求得内层与外层边界（相对于中线延伸 T_neg/2）
T_neg = neg_total_thickness;
plot_r_neg_inner = plot_r_neg - T_neg/2;
plot_r_neg_outer = plot_r_neg + T_neg/2;
plot_x_neg_inner = plot_r_neg_inner .* cos(plot_theta_list_neg);
plot_y_neg_inner = plot_r_neg_inner .* sin(plot_theta_list_neg);
plot_x_neg_outer = plot_r_neg_outer .* cos(plot_theta_list_neg);
plot_y_neg_outer = plot_r_neg_outer .* sin(plot_theta_list_neg);

%% 5. 绘制负极入位
plot(plot_x_neg, plot_y_neg, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'neg-mid');
plot(plot_x_neg_inner, plot_y_neg_inner, '--', 'LineWidth', 1, 'Color', 'r', 'DisplayName', 'neg-inner');
plot(plot_x_neg_outer, plot_y_neg_outer, '--', 'LineWidth', 1, 'Color', 'r', 'DisplayName', 'neg-outer');

%% 6. 绘制负极入位处与隔膜的衔接
% 内层隔膜（图例隐藏）
plot((plot_r_neg_inner - sep_thickness/2) .* cos(plot_theta_list_neg), ...
    (plot_r_neg_inner - sep_thickness/2) .* sin(plot_theta_list_neg), ...
    'LineWidth',2, 'Color', 'b', 'HandleVisibility','off');
plot((plot_r_neg_inner - sep_thickness) .* cos(plot_theta_list_neg), ...
    (plot_r_neg_inner - sep_thickness) .* sin(plot_theta_list_neg), ...
    '--', 'LineWidth',1, 'Color', 'b', 'HandleVisibility','off');
plot(plot_x_neg_inner, plot_y_neg_inner, '--', 'LineWidth',1, 'Color', 'b', 'HandleVisibility','off');

% 外层隔膜（图例隐藏）
plot((plot_r_neg_outer + sep_thickness/2) .* cos(plot_theta_list_neg), ...
    (plot_r_neg_outer - sep_thickness/2) .* sin(plot_theta_list_neg), ...
    'LineWidth',2, 'Color', 'b', 'HandleVisibility','off');
plot((plot_r_neg_outer + sep_thickness) .* cos(plot_theta_list_neg), ...
    (plot_r_neg_outer - sep_thickness) .* sin(plot_theta_list_neg), ...
    '--', 'LineWidth',1, 'Color', 'b', 'HandleVisibility','off');
plot(plot_x_neg_outer, plot_y_neg_outer, '--', 'LineWidth',1, 'Color', 'b', 'HandleVisibility','off');

%% 7. 后续水平——同卷时隔膜、负极和正极
% 同卷时隔膜1与隔膜2（依次取负极内外层末端作为起点）
a_sep1 = plot_r_neg_inner(end) - sep_thickness/2;
b_sep1 = (sep_thickness*2 + neg_thickness + foil_neg_thickness + pos_thickness + foil_pos_thickness) / (2*pi);
a_sep2 = plot_r_neg_outer(end) + sep_thickness/2;
b_sep2 = b_sep1;

% 同卷时负极与正极
a_neg1 = plot_r_neg(end);
b_neg1 = b_sep1;
a_pos = a_sep2 + sep_thickness/2 + (pos_thickness + foil_pos_thickness) / 2;
b_pos = b_sep1;

arc_pos = pos_length;
arc_neg = neg_length;
[theta_pos_end, r_pos_end] = solve_theta(a_pos, b_pos, 0, arc_pos);
[theta_neg_end, r_neg_end] = solve_theta(a_neg1, b_neg1, 0, arc_neg - neg_arc_length);

% 调用子函数绘制各段边界（DisplayName 用于 legend）
plot_edge(a_pos, b_pos, sep_theta + neg_theta, theta_pos_end, (pos_thickness + foil_pos_thickness), 'k', true, 'pos');
plot_edge(a_neg1, b_neg1, sep_theta + neg_theta, theta_neg_end, (neg_thickness + foil_neg_thickness), 'r', true, 'neg');
plot_edge(a_sep1, b_sep1, sep_theta + neg_theta, theta_neg_end + pi/3, sep_thickness, 'g', false, 'sep1');
plot_edge(a_sep2, b_sep2, sep_theta + neg_theta, theta_neg_end + pi/3, sep_thickness, 'm', false, 'sep2');

%% 8. 计算双极耳位置
if calculate_flag=="neg"
    L_total=neg_length;
    [theta_neg_first, r_neg_first] = solve_theta(a_neg1, b_neg1, 0, L_first - neg_arc_length+w/2);
    num=ceil(theta_neg_end/(2*pi)-theta_neg_first/(2*pi)); %剩余可卷绕圈数
    arc_length_second=inf.*ones(num,1);
    for k=1:num
        arc_length_second(k,1)=solve_arc(a_neg1, b_neg1, 0,theta_neg_first+k*2*pi);
    end
    diff=arc_length_second-(L_second+L_first - neg_arc_length + w+w/2);
    [min_diff,index]=min(abs(diff));

    second_true=arc_length_second(index)-w-(L_first - neg_arc_length)-w/2;

    %% 可视化
    %第一留白区中心位置
    [point1_theta, ~] = solve_theta(a_neg1, b_neg1, 0, L_first - neg_arc_length+w/2);
    plot_point1_theta=point1_theta+(sep_theta + neg_theta);
    point1_r=a_neg1+b_neg1*(plot_point1_theta-(sep_theta + neg_theta));
    scatter(point1_r*cos(plot_point1_theta),point1_r*sin(plot_point1_theta),'r','filled',DisplayName='tab1');

    [point2_theta, ~] = solve_theta(a_neg1, b_neg1, 0, L_first-neg_arc_length+w+second_true+w/2);
    plot_point2_theta=point2_theta+(sep_theta + neg_theta);
    point2_r=a_neg1+b_neg1*(plot_point2_theta-(sep_theta + neg_theta));
    scatter(point2_r*cos(plot_point2_theta),point2_r*sin(plot_point2_theta),'r','filled',DisplayName='tab2');

    legend('Location','eastoutside');
    axis square;grid on
    hold off;

elseif calculate_flag=="pos"
    L_total=pos_length;
    [theta_pos_first, r_pos_first] = solve_theta(a_pos, b_pos, 0, L_first +w/2);
    num=ceil(theta_pos_end/(2*pi)-theta_pos_first/(2*pi)); %剩余可卷绕圈数
    arc_length_second=inf.*ones(num,1);
    for k=1:num
        arc_length_second(k,1)=solve_arc(a_pos, b_pos, 0,theta_pos_first+k*2*pi);
    end
    diff=arc_length_second-(L_second+L_first  + w+w/2);
    [min_diff,index]=min(abs(diff));
    second_true=arc_length_second(index)-w-(L_first )-w/2;

    %% 可视化
    %第一留白区中心位置
    [point1_theta, ~] = solve_theta(a_pos, b_pos, 0, L_first +w/2);
    plot_point1_theta=point1_theta+(sep_theta + neg_theta);
    point1_r=a_pos+b_pos*(plot_point1_theta-(sep_theta + neg_theta));
    scatter(point1_r*cos(plot_point1_theta),point1_r*sin(plot_point1_theta),'r','filled',DisplayName='tab1');

    [point2_theta, ~] = solve_theta(a_pos, b_pos, 0, L_first+w+second_true+w/2);
    plot_point2_theta=point2_theta+(sep_theta + neg_theta);
    point2_r=a_neg1+b_neg1*(plot_point2_theta-(sep_theta + neg_theta));
    scatter(point2_r*cos(plot_point2_theta),point2_r*sin(plot_point2_theta),'r','filled',DisplayName='tab2');

    legend('Location','eastoutside');
    axis square;grid on
    hold off;
else
    error('请输入正确的计算标签，neg或者pos')
end


%% 9. 涂布后滚压前数据
%压后尺寸：
L_B_1=L_first;
L_B_2=w;
L_B_3=second_true;
L_B_4=w;
L_B_5=L_total-L_B_1-L_B_2-L_B_3-L_B_4;
%压前尺寸
L_A_1=L_B_1/el_rate;
L_A_2=L_B_2/el_rate;
L_A_3=L_B_3/el_rate;
L_A_4=L_B_4/el_rate;
L_A_5=L_B_5/el_rate;

width=200;
figure
% 初始化 x 坐标位置
x_start = 0;
% 绘制压后尺寸矩形
% subplot(2, 1, 1); % 第一行绘制压后尺寸
rectangle('Position', [x_start, 0, L_B_1, width], 'FaceColor', 'r', 'EdgeColor', 'k');
x_start = x_start + L_B_1;
rectangle('Position', [x_start, 0, L_B_2, width], 'FaceColor', 'g', 'EdgeColor', 'k');
x_start = x_start + L_B_2;
rectangle('Position', [x_start, 0, L_B_3, width], 'FaceColor', 'b', 'EdgeColor', 'k');
x_start = x_start + L_B_3;
rectangle('Position', [x_start, 0, L_B_4, width], 'FaceColor', 'y', 'EdgeColor', 'k');
x_start = x_start + L_B_4;
rectangle('Position', [x_start, 0, L_B_5, width], 'FaceColor', 'm', 'EdgeColor', 'k');
title('压后尺寸矩形排列');
axis equal;
xlabel('Length (mm)');
ylabel('Height (mm)');
% 准备打印的字符内容
line1 = '示意图：   |---L1---|-w-|---L3---|-w-|---L5---|';
line2 = sprintf('压后尺寸：|---%.2f---|-%.2f-|---%.2f---|-%.2f-|---%.2f---|',L_B_1,L_B_2,L_B_3,L_B_4,L_B_5);
line3 = sprintf('压前尺寸：|---%.2f---|-%.2f-|---%.2f---|-%.2f-|---%.2f---|',L_A_1,L_A_2,L_A_3,L_A_4,L_A_5);

% 在绘图窗口中显示文本框
annotation('textbox', [0.15, 0.8, 0.7, 0.1], 'String', {line1, line2,line3}, ...
    'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'left', 'BackgroundColor', 'w');

fprintf('示意图：  |---L1---|-w-|---L3---|-w-|---L5---|\n')
fprintf('压后尺寸：|---%.2f---|-%.2f-|---%.2f---|-%.2f-|---%.2f---|\n',L_B_1,L_B_2,L_B_3,L_B_4,L_B_5)
fprintf('压前尺寸：|---%.2f---|-%.2f-|---%.2f---|-%.2f-|---%.2f---|\n',L_A_1,L_A_2,L_A_3,L_A_4,L_A_5)
end
