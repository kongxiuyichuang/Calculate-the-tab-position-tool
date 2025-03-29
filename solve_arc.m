function [arc_length]=solve_arc(a,b,theta_start,theta_end)

L = @(a,b,theta) (b.^2*theta+a.*b)./(2*b.^2).*sqrt(b.^2*theta.^2+2*a.*b.*theta+a.^2+b.^2)...
    +b./2.*log(abs(b.*theta+a+sqrt(b.^2*theta.^2+2*a.*b.*theta+a.^2+b.^2)));

arc_length = L(a,b,theta_end)-L(a,b,theta_start);
end

% %% 绘图
% theta_plot=linspace(theta_start,theta_end,2000);
% r_plot=a+b.*theta_plot;
% 
% x_plot=r_plot.*cos(theta_plot);
% y_plot=r_plot.*sin(theta_plot);
% 
% figure
% plot(x_plot,y_plot,'-',LineWidth=2)
% xlabel('x');ylabel('y');
% maxR = max(r_plot);
% xlim([-maxR*1.1, maxR*1.1]);
% ylim([-maxR*1.1, maxR*1.1]);
% axis square;
% grid on;
% title(sprintf('arc_{length}=%.3f',arc_length));
