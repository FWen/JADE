function [azm, ele] = delay_sum_init(R, azm_range, ele_range, K_pos, num)
    
%     Rx2 = inv(sqrtm(R));
     Rx2 = sqrtm(R);

    for ke=1:length(ele_range)
        %spec(ke,:) = 1./sum(abs(arst(K_pos,azm_range,ele_range(ke))'*Rx2).^2, 2);
        spec(ke,:) = sum(abs(arst(K_pos,azm_range,ele_range(ke))'*Rx2).^2, 2);
    end
%     [AZ,EL] = meshgrid(azm_range, ele_range);
%     figure(1); mesh(AZ,EL,spec); xlabel('Azimuth (degree)'); ylabel('Elevation (degree)');  

    [x_pos, y_pos, ~] = find2dpeaks(spec);
    azm = sort(azm_range(y_pos(1:num)),'ascend').';
    ele = sort(ele_range(x_pos(1:num)),'ascend').';
end