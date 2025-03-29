function [theta_solution,r_solution]=solve_theta(a,b,theta_start,s_know)

L = @(a,b,theta) (b.^2*theta+a.*b)./(2*b.^2).*sqrt(b.^2*theta.^2+2*a.*b.*theta+a.^2+b.^2)...
    +b./2.*log(abs(b.*theta+a+sqrt(b.^2*theta.^2+2*a.*b.*theta+a.^2+b.^2)));

arc_length =@(a,b,theta_start,theta_end) L(a,b,theta_end)-L(a,b,theta_start);

theta_guess = s_know / sqrt((a + b*0)^2 + b^2);

options = optimset('TolX',1e-6);
theta_solution = fzero(@(theta) arc_length(a,b,theta_start,theta) - s_know, theta_guess, options);
r_solution=a+b*theta_solution;
end
% %% 绘图
% theta_plot=linspace(theta_start,theta_solution,2000);
% r_plot=a+b.*theta_plot;
% x_plot=r_plot.*cos(theta_plot);
% y_plot=r_plot.*sin(theta_plot);
% figure
% plot(x_plot,y_plot,'-',LineWidth=2)
% xlabel('x');ylabel('y');
% 
% maxR = max(r_plot);
% xlim([-maxR*1.1, maxR*1.1]);
% ylim([-maxR*1.1, maxR*1.1]);
% axis square;
% grid on;
% title(sprintf('theta_{solution}=%.3f(%.2f pi) r_{solution}=%.3f',theta_solution,theta_solution/pi,r_solution))
