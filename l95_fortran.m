function s = l95_fortran(time,init,params,opt,delete_files)
    if nargin < 5
        delete_files = true;
    end
    ntime = length(time);
    [K,nens] = size(init);
    save time.mat time -v4;
    save init.mat init -v4;
    save params.mat params -v4;
    save opt.mat opt -v4
    !ode_solver.exe time.mat init.mat params.mat opt.mat
    % load results and reshape
    load solution.mat;
    s = reshape(solution,ntime,K,nens);
    if delete_files 
        delete('time.mat','init.mat','params.mat','opt.mat','solution.mat');
    end
end