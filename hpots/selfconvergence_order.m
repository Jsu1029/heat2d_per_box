
% self convergence for heat solver
format long

pot = load('self_convergence_test.data');

nt = 4; %number of calling heat solver
ntarg=1000000; %number of targets

pot_num = zeros(nt,ntarg);
for i = 1:nt
    pot_num(i,1:ntarg)=pot((i-1)*ntarg+1:ntarg*i);
end

k=zeros(nt-1,1);

for i = 1:nt-1
    err=0.0;
    err_scale=0.0;
    for ii =1:ntarg
        err=err+abs(pot_num(i,ii)-pot_num(i+1,ii))^2;
        err_scale=err_scale+abs(pot_num(i,ii))^2;
    end
    k(i)=sqrt(err)/sqrt(err_scale);
end


    order = k(1:end-1) ./ k(2:end);
    log2(order)

