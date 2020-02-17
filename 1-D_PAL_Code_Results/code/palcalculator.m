function [I]=palcalculator(y,maxcomp,numsamp);

yadj=y;
comp=maxcomp;
gridnum=1000; % number of grid points for periodogram maximizer 

for ks=1:maxcomp
   inits=per_max_c(yadj,numsamp,gridnum); 
   seq_a(ks)=inits(1);
   seq_om(ks)=inits(2);
for ns=1:numsamp
    yadj(ns)=(yadj(ns)-inits(1)*exp(i*2*pi*inits(2)*ns));
end

end

initial=[seq_a.',seq_om'];
L2_full=fminsearch('obj_L2_C',initial,[],y,numsamp,comp);

sigsq_null=0;
sigsq_full=0;
for k=1:numsamp
    yp(k)=0.0;
    for kp=1:maxcomp
        yp(k)=yp(k)+(L2_full(kp,1)*exp(i*2*pi*L2_full(kp,2)*k));
    end
sigsq_null=sigsq_null+abs(y(k))^2;    
diff_full_sq(k)=abs(y(k)-yp(k))^2;
sigsq_full=sigsq_full+diff_full_sq(k);
end
sigsq_null=sigsq_null/(numsamp);
sigsq_full=sigsq_full/(numsamp);


for comp=1:maxcomp
    nummax=comp;
yadj=y;

gridnum=1000; % number of grid points for periodogram maximizer 

for ks=1:nummax
   inits=per_max_c(yadj,numsamp,gridnum);
   sequ_a(ks)=inits(1);
   sequ_om(ks)=inits(2);
for ns=1:numsamp
    yadj(ns)=(yadj(ns)-inits(1)*exp(i*2*pi*inits(2)*ns));
end
end

initial_m=[sequ_a.',sequ_om'];
L2_model=fminsearch('obj_L2_C',initial_m,[],y,numsamp,comp);

sum_sq=0.0;
for k=1:numsamp
    yp(k)=0.0;
    for kp=1:comp
        yp(k)=yp(k)+(L2_model(kp,1)*exp(i*2*pi*L2_model(kp,2)*k));
    end
        diff(k)=(y(k)-yp(k));
        sum_sq=sum_sq+abs(diff(k))^2;

end
est_var=sum_sq/(numsamp);
pal_est_sigsq(comp)=est_var;


if comp==1
    pal_rn=2*numsamp*log(sigsq_null/sigsq_null);
    pal_rhon=2*numsamp*log(sigsq_null/sigsq_full);
else
    pal_rn=2*numsamp*log(sigsq_null/pal_est_sigsq(comp-1));
    pal_rhon=2*numsamp*log(pal_est_sigsq(comp-1)/sigsq_full);
end

%bic_ic_obj(comp)=2*numsamp*log(est_var)+(log(numsamp))*(3*comp+1);
%bic_cor_obj(comp)=2*numsamp*log(est_var)+(log(numsamp))*(5*comp+1);
pal_obj(comp)=2*numsamp*log(est_var)+(3*comp+1)*log(3*maxcomp+1)*(log((pal_rn)+1)/log((pal_rhon)+1));
%aic_obj(comp)=2*numsamp*log(est_var)+2*(3*comp+1);
%aicc_obj(comp)=numsamp*log(est_var)+2*(3*comp*numsamp)/(numsamp-3*comp-1);
%ca=2*(3*comp*numsamp)/(numsamp-3*comp-1); cb=(log(numsamp))*(5*comp);  % check the expression
%wic_obj(comp)=((ca)/(ca+cb))*aicc_obj(comp)+((cb)/(ca+cb))*bic_cor_obj(comp);
%hq_obj(comp)=2*numsamp*log(est_var)+2*(3*comp+1)*log(log(numsamp));

end % loop for all possible number of components
% PAL
[CR,I]=min(pal_obj);
estcomp=I;