function [obj_L2_C] = fun(initial,y,numsamp,comp)

for pa=1:comp
    a(pa)=initial(pa,1);
    om(pa)=initial(pa,2);
end
obj_L2_C=0.0;
for ln=1:numsamp
    yp(ln)=0.0;
    for lp=1:comp  
        yp(ln)=yp(ln)+(a(lp)*exp(i*2*pi*om(lp)*ln));
    end
    diff_sq(ln)=abs(y(ln)-yp(ln))^2;
    obj_L2_C=obj_L2_C+diff_sq(ln);
end

% CALCULATION OF VANDERMONDE MATRIX AND CONCENTRATED LIKELIHOOD
% for lr=1:numsamp
%  for lc=1:numpar
%   a(lr,lc)=cos(freqd(lc*exp(i*trfreq(lc)*lr);
%  end;
% end;
% CALCULATION OF OBJECTIVE FUNCTION
%obj=-(yd'*a*pinv(a'*a)*a'*yd);