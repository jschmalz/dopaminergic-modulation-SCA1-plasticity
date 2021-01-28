function [slope_EPSP, F_t] = make_slopeEPSP(v,time,freq)
dt = time(2) - time(1);
len = floor(dt*size(v,1)/(freq))
slope_EPSP = zeros(len,1);
F_t = zeros(len,1);

% find first EPSP ignoring AP
check = 1;
for j = (60*1000-100)/dt:(60*1000+100)/dt
    if v(j+2) > v(j)*0.999 && check == 1 
        slope_EPSP(1) = (v(j+2) - v(j))/(dt*1);
        F_t(1) = j*dt;
        check = 0;
    end
end


for i = 2:len
    a = ((i-1)*freq/dt-100/dt);
    check = 1;
    for j = a:a+100/dt
        if v(j+2) > v(j)*0.999 && check == 1
            slope_EPSP(i) = (v(j+2) - v(j))/(dt*1);
            F_t(i) = (i-1)*freq;         
            check = 0;
        end
    end
end

end