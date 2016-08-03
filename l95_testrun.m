clear all, close all, format shortg, format compact, clc;
nens = 2;   % number of ensemble elements
time = (0:0.1:20)';
sigma_x = 0.1;
sigma_y = 0.01;
sigma_init = 0.1;
sigma_F = 0.1;
sigma_Fy = 0.1;

K = 40;
J = 8;

x0 = 3*ones(K,1)+sigma_x*randn(K,1);
y0 = 3*ones(K*J,1)+sigma_y*randn(K*J,1);
x0ens = repmat(x0,1,nens).*(ones(K,nens)+repmat(sigma_init*randn(1,nens),K,1));
y0ens = repmat(y0,1,nens).*(ones(K*J,nens)+sigma_init*randn(K*J,nens));

F_npar = 40;
F = 8+sigma_F*randn(F_npar,nens);
F = repmat(F,floor(K/F_npar),1);

Fy_npar = 320;
Fy = 8 + sigma_Fy*randn(Fy_npar,nens);
Fy = repmat(Fy,floor(K*J/Fy_npar),1);

b = 5.0;
c = 1.0;
h = 0.1;

dt = 0.0025;

init = [x0ens;y0ens];
params = [F;Fy;b*ones(1,nens);c*ones(1,nens);h*ones(1,nens)];
opt = dt;

delete_files = false;

l95_out = l95_fortran(time,init,params,opt,delete_files);

figure; hold on;
plot(time,l95_out(:,1,1),'r'); 
plot(time,l95_out(:,41,1),'k');
legend('slow var','fast var');
grid on;
axis tight;





% function create_test_data(time,nstates,nparams,nens,nopt)
%     init_ens = eye(nstates,nens);
%     params_ens = eye(nparams,nens);
%     opt = ones(1,nopt);
%     
%     save time.mat time -v4;
%     save init_ens.mat init_ens -v4;
%     save params_ens.mat params_ens -v4;
%     save opt.mat opt -v4;
% end
    