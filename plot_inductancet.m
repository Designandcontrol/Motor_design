num_iter=21;
figure('Name','i_d');
for x = 1:num_iter
    plot(L.i_d(x,:))
    hold on
end
figure('Name','i_q');
title('i_q');
for x = 1:num_iter
    plot(L.i_q(x,:))
    hold on
end

figure('Name','flux_q');
for x = 1:num_iter
    plot(L.i_q(:,x),L.flux_q(:,x))
    hold on
    xlabel('L.i_q')
    ylabel('flux_q')
    grid on
end
figure('Name','flux_d');
for x = 1:num_iter
    plot(L.i_d(x,:),L.flux_d(x,:))
    hold on
    xlabel('L.i_d')
    ylabel('flux_d')
    grid on
end
figure('Name','L_dd');
for x = 1:num_iter
      
    plot(L.i_d(x,:), L.L_dd(x, :));
    hold on
    xlabel('L.i_d')
    ylabel('L_dd,L_qq')
    grid on
     ylim([0 1])
end
figure('Name','L_qq');
for x = 1:num_iter
       plot(L.i_q(:,x), L.L_qq(:, x));
    hold on
    xlabel('L.i_q')
    ylabel('L_qq')
    grid on
end
figure('Name','L_qd');
for x = 1:num_iter

    plot(L.i_d(x,:), L.L_qd(:, x));
    hold on
    xlabel('L.i_d')
    ylabel('L_qd')
    grid on
end
figure('Name','L_dq');
for x = 1:num_iter
   
        % 좌변과 우변의 크기를 맞추기 위해 L.i_d 배열의 크기를 21x21로 만듭니다.
    plot(L.i_q(:,x), L.L_dq(x, :));
    hold on
    xlabel('L.i_d')
    ylabel('L_dq')
    grid on

end

%%%%%% 3차원 plot %%%%%%%%%%%%%%%%%
figure('Name','L_dd');

  
    surf(L.i_d,L.i_q, L.L_dd);
    hold on
    xlabel('i_d')
    ylabel('i_q')
    zlabel('L.dd')
    grid on
    
figure('Name','L_qq');

    surf(L.i_q,L.i_d, L.L_qq);
    hold on
    xlabel('i_q')
    ylabel('i_d')
    zlabel('L.qq')
    grid on
    
figure('Name','L_dq');

    
    surf(L.i_q,L.i_d,L.L_dq);
    hold on
    xlabel('i_q')
    ylabel('i_d')
    zlabel('L.dq')
    grid on
    
figure('Name','L_qd');

    surf(L.i_d,L.i_q, L.L_qd);
    hold on
    xlabel('i_d')
    ylabel('i_q')
    zlabel('L.qd')
    grid on
        
