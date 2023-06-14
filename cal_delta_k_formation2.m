    for i=1:n
        x0=(i-1)*la/n;
        x1=i*la/n;
        f(i)=@(y) 1/2*(rou_up+rou_low)*(x1-x0)-1/2*(rou_up-rou_low)*(cos(x1+y)-cos(x0+y));
    end
    y2=sym('y2'); 

    while pos_y_line<lb
        pos_i = 2*[lb, lb, lb, lb, lb, lb, lb];
        tic
        for i=1:n
            if i==1
%                 f1=int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,y2)==Ts*v;
                pos_i(i) = vpasolve(int(f(i),y,pos_y_line,y2)==Ts*v,y2);
            elseif integral2(rou,(i-1)*la/n,i*la/n,pos_y_line,min(pos_i)) > Ts*v
%                 f1=int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,y2)==Ts*v;
                pos_i(i) = vpasolve(int(f(i),y,pos_y_line,y2)==Ts*v,y2);
            else
                pos_i(i) = min(pos_i);
            end                
        end
        toc

        k_line = k_line+1;
        pos_y_line = min(pos_i);
    end
    fprintf('for Ts = %f, Delta_K_line is', Ts);
    Delta_K_line = k_line-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))