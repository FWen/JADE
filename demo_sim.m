
% This is a demo for CSI based joint azimuth, elevation and timde delay estimation

clear all; %close all; clc;

fc = 5.32e9;
delta_f = 312.5e3;
sub_fc = [-64:1:63];

lambda = 3e8/fc;
M = 16;
K = length(sub_fc);

d2pi = pi/180;
theta = [30, 40];
phi = [50, 60];
tof = [50, 100]*1e-9;
L = length(theta);

% attenuation factors
beta = exp(1i*2*pi*rand(1,L));

% position matrix of the array for UCA
R = 1.5*lambda;
array_pos = R*[cos(2*pi*(0:M-1)'/M), sin(2*pi*(0:M-1)'/M)];

K_tof = -1i*2*pi*(fc+sub_fc*delta_f).';
K_pos = -1i*2*pi/lambda*array_pos;

% Array manifold
A = arst(K_pos, theta, phi);
F = exp(K_tof*tof) * diag(beta);

X0 = A*F.';

SNR = -15:2.5:20;

for snr=1:length(SNR)
    
    % ---CRB-----------------------------------
    % CRB for joint DOA and TD estimation
    rho_az = [-sin(theta*d2pi).*cos(phi*d2pi); cos(theta*d2pi).*cos(phi*d2pi)];
    rho_el = [-cos(theta*d2pi).*sin(phi*d2pi); -sin(theta*d2pi).*sin(phi*d2pi)];
    
    A_az = K_pos*rho_az.*A;
    A_el = K_pos*rho_el.*A;
    for kf=1:K
        D((kf-1)*M+1:kf*M,:) = A*diag(exp(K_tof(kf)*tof));
        E((kf-1)*M+1:kf*M,:) = [A_az*diag(exp(K_tof(kf)*tof)), A_el*diag(exp(K_tof(kf)*tof))];
        G((kf-1)*M+1:kf*M,:) = A*diag(K_tof(kf)*exp(K_tof(kf)*tof));
    end
    
    P_D = eye(M*K) - D*inv(D'*D)*D';
    F1 = real( (E'*P_D*E).*([beta, beta]'*[beta, beta]) );
    F2 = real( (E'*P_D*G).*([beta, beta]'*beta) );
    F3 = real( (G'*P_D*G).*(beta'*beta) );

    thou_n = var(X0(:))*10^(-SNR(snr)/10);
    CRB_ang = sqrt(thou_n/2*inv(F1-F2*inv(F3)*F2.'));
    CRB_tof = sqrt(thou_n/2*inv(F3-F2.'*inv(F1)*F2));
    CRB_J(:,snr) = [diag(CRB_ang)/d2pi;diag(CRB_tof)*1e9]; 
    
    % CRB for DOA-only estimation
    PHI = [A_az, A_el];
    Pc = [F,F]'*[F,F];
    F4 = real((PHI'*(eye(M)-A*inv(A'*A)*A')*PHI) .* Pc);
    CRB_O(:,snr) = diag( sqrt(thou_n/2*inv(F4)) )/d2pi;
    % ---CRBs End--------------------------------------
    
    
    for mc=1:100
        disp(sprintf('SNR: %.1f,  turn: %d',SNR(snr),mc));
        
        X = awgn(X0, SNR(snr), 'measured');
        
       
        %-- ML algorithm for DOA-only estimation--------------
        theta_e = theta + 2*rand(size(theta));
        phi_e   = phi + 2*rand(size(theta));
        Rs = X*X'/K;
        if SNR(snr)<=-10
            [theta_e, phi_e] = delay_sum_init(Rs, 0:2:180, 0:2:90, K_pos, L);
        end
        if SNR(snr)>10
            theta_e = theta + 1*rand(size(theta));
            phi_e   = phi + 1*rand(size(theta));
        end

        detaDOA = [2, 0.2, 0.02];
        for iter=1:3 % iteration, alternating between the two paths
            theta_e = flip(theta_e); 
            phi_e   = flip(phi_e); 
            for step=1:length(detaDOA)
                azm_range = theta_e(1)-detaDOA(step)*5:detaDOA(step):theta_e(1)+detaDOA(step)*5;
                ele_range = phi_e(1)-detaDOA(step)*5:detaDOA(step):phi_e(1)+detaDOA(step)*5;
                for ka=1:length(azm_range)
                    for ke=1:length(ele_range)
                        As = [arst(K_pos,azm_range(ka),ele_range(ke)), arst(K_pos,theta_e(2),phi_e(2))];
                        spec1(ke,ka) = real(trace((As*inv(As'*As)*As')*Rs));
                    end
                end
                [theta_e(1), phi_e(1)] = find_max_peak(spec1, azm_range, ele_range);
            end
        end
        ml_azm(mc,:) = sort(theta_e,'ascend');
        ml_ele(mc,:) = sort(phi_e,'ascend');



        %-- AML algorithm for jiont DOA and TD estimation--------------
        iters = [2,4];
        for n_iter=1:length(iters)
            theta_l = [];  phi_l   = [];
            tof_l   = [];  beta_l  = [];
            for l=1:L

                % initialization
                if isempty(theta_l)
                    X_res = X;
                else
                    for k=1:K
                        X_res(:,k) = X(:,k) - arst(K_pos,theta_l,phi_l)*(exp(K_tof(k)*tof_l).*beta_l);
                    end
                end
                Rx  = X_res*X_res'/K;
                Rx2 = sqrtm(Rx);
                azm_range = 0:10:180;
                ele_range = 0:10:89;
                for ke=1:length(ele_range)
                    spec(ke,:) = sum(abs(arst(K_pos,azm_range,ele_range(ke))'*Rx2).^2, 2);
                end
                [azm, ele] = find_max_peak(spec, azm_range, ele_range);
    
                theta_l = [theta_l, azm];
                phi_l   = [phi_l, ele];
                for rep=1:iters(n_iter)
                    A_l = arst(K_pos, theta_l, phi_l);
                    U = inv(A_l'*A_l)*A_l'*X;

                    %--- time delay estimation-------------
                    detaTOF = [10, 1, 0.1]*1e-9;
                    rangTOF = [0, 1000; 5, 5; 0.5, 0.5]*1e-9;
                    for l2=1:size(U,1)
                        tof_l(l2) = 0;
                        for step=1:length(detaTOF)
                            tof_range = tof_l(l2)-rangTOF(step,1):detaTOF(step):tof_l(l2)+rangTOF(step,2);
                            [mv,mi] = max(abs(exp(K_tof*tof_range)'*U(l2,:).'));
                            tof_l(l2) = tof_range(mi);
                        end
                        beta_l(l2) = exp(K_tof*tof_l(l2))' * U(l2,:).'/K;
                    end

                    %--- angle estimaton--------------------
                    rr = exp(-1i*2*pi*(fc+sub_fc*delta_f).'*tof_l).' ;
                    B = X*rr'*inv(rr*rr');
                    detaDOA = [2, 0.2, 0.02];
                    for l3=1:size(B,2)
                        for step=1:length(detaDOA)
                            azm_range = theta_l(l3)-detaDOA(step)*5:detaDOA(step):theta_l(l3)+detaDOA(step)*5;
                            ele_range = phi_l(l3)-detaDOA(step)*5:detaDOA(step):phi_l(l3)+detaDOA(step)*5;
                            [theta_l(l3), phi_l(l3)] = az_el_estimation(B(:,l3), azm_range, ele_range, K_pos);
                        end                    
                    end
                end
            end
            aml_azm(mc,(1:L)+(n_iter-1)*L) = sort(theta_l,'ascend');
            aml_ele(mc,(1:L)+(n_iter-1)*L) = sort(phi_l,'ascend');
            aml_tof(mc,(1:L)+(n_iter-1)*L) = sort(tof_l,'ascend');
        end
    end
    RMSE_aml_azm(:,snr) = sqrt(mean((aml_azm(:,:) - repmat(theta,[mc,length(iters)])).^2))';
    RMSE_aml_ele(:,snr) = sqrt(mean((aml_ele(:,:) - repmat(phi,[mc,length(iters)])).^2))';
    RMSE_aml_tof(:,snr) = sqrt(mean((aml_tof(:,:) - repmat(tof,[mc,length(iters)])).^2))'*1e9;
    RMSE_ml_azm(:,snr)  = sqrt(mean((ml_azm(:,:) - repmat(theta,[mc,1])).^2))';    % DOA only est
    RMSE_ml_ele(:,snr)  = sqrt(mean((ml_ele(:,:) - repmat(phi,[mc,1])).^2))';      % DOA only est
end

figure(1); 
subplot(2,3,1);
semilogy(SNR,RMSE_ml_azm(1,:),'g--+',SNR,RMSE_aml_azm(1,:),'b--x',SNR,RMSE_aml_azm(3,:),'b--o',SNR,CRB_O(1,:),'g--',SNR,CRB_J(1,:),'r-','linewidth',1);grid on;
legend('ML (DOA-only)','AML (2 iter.)','AML (converged)','CRB (DOA-only)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (deg)');title('\theta_1');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(1,:))/2,max(RMSE_aml_azm(1,:))*2]);
subplot(2,3,2);
semilogy(SNR,RMSE_ml_ele(1,:),'g--+',SNR,RMSE_aml_ele(1,:),'b--x',SNR,RMSE_aml_ele(3,:),'b--o',SNR,CRB_O(3,:),'g--',SNR,CRB_J(3,:),'r-','linewidth',1);grid on;
legend('ML (DOA-only)','AML (2 iter.)','AML (converged)','CRB (DOA-only)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (deg)');title('\phi_1');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(3,:))/2,max(RMSE_aml_ele(1,:))*2]);
subplot(2,3,3);
semilogy(SNR,RMSE_aml_tof(1,:),'b--x',SNR,RMSE_aml_tof(3,:),'b--o',SNR,CRB_J(5,:),'r-','linewidth',1);grid on;
legend('AML (2 iter.)','AML (converged)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (ns)');title('\tau_1');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(5,:))/2,max(RMSE_aml_tof(1,:))*2]);
subplot(2,3,4);
semilogy(SNR,RMSE_ml_azm(2,:),'g--+',SNR,RMSE_aml_azm(2,:),'b--x',SNR,RMSE_aml_azm(4,:),'b--o',SNR,CRB_O(2,:),'g--',SNR,CRB_J(2,:),'r-','linewidth',1);grid on;
legend('ML (DOA-only)','AML (2 iter.)','AML (converged)','CRB (DOA-only)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (deg)');title('\theta_2');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(2,:))/2,max(RMSE_aml_azm(2,:))*2]);
subplot(2,3,5);
semilogy(SNR,RMSE_ml_ele(2,:),'g--+',SNR,RMSE_aml_ele(2,:),'b--x',SNR,RMSE_aml_ele(4,:),'b--o',SNR,CRB_O(4,:),'g--',SNR,CRB_J(4,:),'r-','linewidth',1);grid on;
legend('ML (DOA-only)','AML (2 iter.)','AML (converged)','CRB (DOA-only)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (deg)');title('\phi_2');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(4,:))/2,max(RMSE_aml_ele(2,:))*2]);
subplot(2,3,6);
semilogy(SNR,RMSE_aml_tof(2,:),'b--x',SNR,RMSE_aml_tof(4,:),'b--o',SNR,CRB_J(5,:),'r-','linewidth',1);grid on;
legend('AML (2 iter.)','AML (converged)','CRB (joint est.)','Location','NorthEast');
xlabel('SNR (dB)');ylabel('RMSE (ns)');title('\tau_2');
xlim([SNR(1),SNR(end)]);ylim([min(CRB_J(6,:))/2,max(RMSE_aml_tof(2,:))*2]);
