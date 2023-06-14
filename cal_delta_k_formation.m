    while pos_y_line<lb
        pos_i = 2*[lb, lb, lb, lb, lb, lb, lb];
        tic
        parfor i=1:n
            y2=sym('y2'); 
            f1=int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,y2)==Ts*v;
            pos_i(i) = vpasolve(f1,y2, pos_y_line+Ts*v/(((rou_up+rou_low)/2)*la/n));              
        end
        toc
        k_line = k_line+1;
        pos_y_line = min(pos_i);
    end
    fprintf('for Ts = %f, Delta_K_line is', Ts);
    Delta_K_line = k_line-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))