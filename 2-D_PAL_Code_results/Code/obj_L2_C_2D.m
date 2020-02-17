function [obj_L2_C] = fun(initial,y,numsamp,comp)

for pa=1:comp
    a(pa)=initial(pa,1);
    om(pa)=initial(pa,2);
    be(pa)=initial(pa,3);
end
obj_L2_C=0.0;
for ls=1:numsamp
    for lt=1:numsamp
        yp((lt-1)*numsamp+ls)=0.0;
        for lp=1:comp  
            yp((lt-1)*numsamp+ls)=yp((lt-1)*numsamp+ls)+(a(lp)*exp(i*2*pi*(om(lp)*ls+be(lp)*lt)));
        end
        diff_sq((lt-1)*numsamp+ls)=abs(y((lt-1)*numsamp+ls)-yp((lt-1)*numsamp+ls))^2;
        obj_L2_C=obj_L2_C+diff_sq((lt-1)*numsamp+ls);
    end
end

% CALCULATION OF VANDERMONDE MATRIX AND CONCENTRATED LIKELIHOOD
% for lr=1:numsamp
%  for lc=1:numpar
%   a(lr,lc)=cos(freqd(lc*exp(i*trfreq(lc)*lr);
%  end;
% end;
% CALCULATION OF OBJECTIVE FUNCTION
%obj=-(yd'*a*pinv(a'*a)*a'*yd);