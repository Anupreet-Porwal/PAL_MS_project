
% MODEL SELECTION OF CISOIDS WITH PAL 
% January 10,2017
% CISOIDS MODEL SIMULATIONS
clear;

%a(1)=7+i*3;
%a(2)=5+4*i;
%a(3)=3+2*i;

a(1)=1+i*sqrt(2);
a(2)=2+2*i;
om(1)=0.13; om(2)=0.31; 
be(1)=0.13; be(2)=0.31;  % 5 4 
freqd=[om; be];%*pi;
sz=size(a);
numpar=sz(1,2);
count=0;
maxcomp=5; % maximum number of components
% for numsamp =100;%:25:75;
for numsamp =50:50:100;    
% 'Sample Size',numsamp
sigsq=15;
df=2; % df 2 for t error; 

stdev1=sqrt(1); stdev2=sqrt(4); 
% Number of simulations
nsim=100;

%Initialization of counters

for comp=1:maxcomp
   bic(comp)=0;
   bicc(comp)=0;

   pal(comp)=0;
   eef(comp)=0;
   aic(comp)=0;
   hq(comp)=0;
end


for isim=1:nsim
   'simulation #'; isim

%Data generation

%independent normal error
y=data_n_c_2D(numsamp,numpar,a,freqd,sigsq);      

%independent t error
%y=data_t(numsamp,a,b,om,df,numpar);      

% normal mixture error
%y=data_n_mix(numsamp,numpar,a,b,om,stdev1,stdev2); 

%Estimation of sigma from full model
yadj=y;

gridnum=100; % number of grid points for periodogram maximizer 

for ks=1:maxcomp
   inits=per_max_c_2D(yadj,numsamp,gridnum); 
   seq_a(ks)=inits(1);
   seq_om(ks)=inits(2);
   seq_be(ks)=inits(3);

for ns=1:numsamp
    for nt=1:numsamp
    yadj((nt-1)*numsamp+ns)=(yadj((nt-1)*numsamp+ns)-inits(1)*exp(i*2*pi*(inits(2)*ns+inits(3)*nt)));
    end
end

end

initial=[seq_a.',seq_om',seq_be']; 
L2_full=fminsearch('obj_L2_C_2D',initial,[],y,numsamp,comp);
sigsq_null=0;
sigsq_full=0;
for k=1:numsamp
    for l=1:numsamp
    yp((l-1)*numsamp+k)=0.0;
    for kp=1:maxcomp
        yp((l-1)*numsamp+k)=yp((l-1)*numsamp+k)+(L2_full(kp,1)*exp(i*2*pi*(L2_full(kp,2)*k+L2_full(kp,3)*l)));
    end
    sigsq_null=sigsq_null+abs(y((l-1)*numsamp+k))^2;    
    diff_full_sq((l-1)*numsamp+k)=abs(y((l-1)*numsamp+k)-yp((l-1)*numsamp+k))^2;
    sigsq_full=sigsq_full+diff_full_sq((l-1)*numsamp+k);
    end
end
sigsq_null=sigsq_null/(numsamp*numsamp);
sigsq_full=sigsq_full/(numsamp*numsamp);


% MODEL SELECTION USING NONROBUST METHODS
for comp=1:maxcomp
    nummax=comp;
yadj=y;

gridnum=100; % number of grid points for periodogram maximizer 

for ks=1:nummax
   inits=per_max_c_2D(yadj,numsamp,gridnum);
   sequ_a(ks)=inits(1);
   sequ_om(ks)=inits(2);
   sequ_be(ks)=inits(3);
for ns=1:numsamp
    for nt=1:numsamp
        yadj((nt-1)*numsamp+ns)=yadj((nt-1)*numsamp+ns)-inits(1)*exp(i*2*pi*(inits(2)*ns+inits(3)*nt));
    end
end
end

initial_m=[sequ_a.',sequ_om',sequ_be'];
L2_model=fminsearch('obj_L2_C_2D',initial_m,[],y,numsamp,comp);

sum_sq=0.0;
for k=1:numsamp
    for l=1:numsamp
        yp((l-1)*numsamp+k)=0.0;
        for kp=1:comp
            yp((l-1)*numsamp+k)=yp((l-1)*numsamp+k)+(L2_model(kp,1)*exp(i*2*pi*(L2_model(kp,2)*k+L2_model(kp,3)*l)));
        end
            diff((l-1)*numsamp+k)=(y((l-1)*numsamp+k)-yp((l-1)*numsamp+k));
            sum_sq=sum_sq+abs(diff((l-1)*numsamp+k))^2;
    end
end
est_var=sum_sq/(numsamp*numsamp);
pal_est_sigsq(comp)=est_var;

if comp==1
    pal_rn=2*numsamp*numsamp*log(sigsq_null/sigsq_null);
    pal_rhon=2*numsamp*numsamp*log(sigsq_null/sigsq_full);
else
    pal_rn=2*numsamp*numsamp*log(sigsq_null/pal_est_sigsq(comp-1));
    pal_rhon=2*numsamp*numsamp*log(pal_est_sigsq(comp-1)/sigsq_full);
end

bic_ic_obj(comp)=2*numsamp*numsamp*log(est_var)+(log(numsamp*numsamp))*(4*comp+1);
%bic_cor_obj(comp)=2*numsamp*log(est_var)+(log(numsamp))*(5*comp+1);
pal_obj(comp)=2*numsamp*numsamp*log(est_var)+(4*comp+1)*log(4*maxcomp+1)*(log((pal_rn)+1)/log((pal_rhon)+1));
aic_obj(comp)=2*numsamp*numsamp*log(est_var)+2*(4*comp+1);
%aicc_obj(comp)=numsamp*log(est_var)+2*(3*comp*numsamp)/(numsamp-3*comp-1);
%ca=2*(3*comp*numsamp)/(numsamp-3*comp-1); cb=(log(numsamp))*(5*comp);  % check the expression
%wic_obj(comp)=((ca)/(ca+cb))*aicc_obj(comp)+((cb)/(ca+cb))*bic_cor_obj(comp);
%hq_obj(comp)=2*numsamp*log(est_var)+2*(3*comp+1)*log(log(numsamp));

end % loop for all possible number of components

% USUAL BIC
[CR,I]=min(bic_ic_obj);
estcomp=I;
%'USUAL BIC', estcomp

    if estcomp == 1
        bic(1)=bic(1)+1;
    elseif estcomp == 2
        bic(2)=bic(2)+1;
    elseif estcomp == 3
        bic(3)=bic(3)+1;
    elseif estcomp == 4
        bic(4)=bic(4)+1;
    elseif estcomp == 5
        bic(5)=bic(5)+1;
    elseif estcomp == 6
        bic(6)=bic(6)+1;
    elseif estcomp == 7
        bic(7)=bic(7)+1;
    elseif estcomp == 8
        bic(8)=bic(8)+1;
    elseif estcomp == 9
        bic(9)=bic(9)+1;
    else estcomp == 10
        bic(10)=bic(10)+1;
    end

% CORRECTED BIC/MAP
% [CR,I]=min(bic_cor_obj);
% estcomp=I;
% %'CORRECTED BIC', estcomp
% 
% 
%     if estcomp == 1
%         bicc(1)=bicc(1)+1;
%     elseif estcomp == 2
%         bicc(2)=bicc(2)+1;
%     elseif estcomp == 3
%         bicc(3)=bicc(3)+1;
%     elseif estcomp == 4
%         bicc(4)=bicc(4)+1;
%     elseif estcomp == 5
%         bicc(5)=bicc(5)+1;
%     elseif estcomp == 6
%         bicc(6)=bicc(6)+1;
%     elseif estcomp == 7
%         bicc(7)=bicc(7)+1;
%     elseif estcomp == 8
%         bicc(8)=bicc(8)+1;
%     elseif estcomp == 9
%         bicc(9)=bicc(9)+1;
%     else estcomp == 10
%         bicc(10)=bicc(10)+1;
%     end
% 

% PAL
[CR,I]=min(pal_obj);
estcomp=I;
%'PAL', estcomp


    if estcomp == 1
        pal(1)=pal(1)+1;
    elseif estcomp == 2
        pal(2)=pal(2)+1;
    elseif estcomp == 3
        pal(3)=pal(3)+1;
    elseif estcomp == 4
        pal(4)=pal(4)+1;
    elseif estcomp == 5
        pal(5)=pal(5)+1;
    elseif estcomp == 6
        pal(6)=pal(6)+1;
    elseif estcomp == 7
        pal(7)=pal(7)+1;
    elseif estcomp == 8
        pal(8)=pal(8)+1;
    elseif estcomp == 9
        pal(9)=pal(9)+1;
    else estcomp == 10
        pal(10)=pal(10)+1;
    end

    
% USUAL AIC
[CR,I]=min(aic_obj);
estcomp=I;
%'USUAL AIC', estcomp


    if estcomp == 1
        aic(1)=aic(1)+1;
    elseif estcomp == 2
        aic(2)=aic(2)+1;
    elseif estcomp == 3
        aic(3)=aic(3)+1;
    elseif estcomp == 4
        aic(4)=aic(4)+1;
    elseif estcomp == 5
        aic(5)=aic(5)+1;
    elseif estcomp == 6
        aic(6)=aic(6)+1;
    elseif estcomp == 7
        aic(7)=aic(7)+1;
    elseif estcomp == 8
        aic(8)=aic(8)+1;
    elseif estcomp == 9
        aic(9)=aic(9)+1;
    else estcomp == 10
        aic(10)=aic(10)+1;
    end


% USUAL HANNAN QUINN
% [CR,I]=min(hq_obj);
% estcomp=I;
% 'USUAL HANNAN QUINN', estcomp
% 
% 
%     if estcomp == 1
%         hq(1)=hq(1)+1;
%     elseif estcomp == 2
%         hq(2)=hq(2)+1;
%     elseif estcomp == 3
%         hq(3)=hq(3)+1;
%     elseif estcomp == 4
%         hq(4)=hq(4)+1;
%     elseif estcomp == 5
%         hq(5)=hq(5)+1;
%     elseif estcomp == 6
%         hq(6)=hq(6)+1;
%     elseif estcomp == 7
%         hq(7)=hq(7)+1;
%     elseif estcomp == 8
%         hq(8)=hq(8)+1;
%     elseif estcomp == 9
%         hq(9)=hq(9)+1;
%     else estcomp == 10
%         hq(10)=hq(10)+1;
%     end

result_all=[bic;pal;aic]
end % Loop for number of simulations
% '-----------------------------------------------------------------------',
% 
% '-----------------------------------------------------------------------',
% 
% 
% 'USUAL BIC'
% bic,
% 
% '-----------------------------------------------------------------------',
% 'CORRECTED BIC'
% bicc,
% 
% 
% '-----------------------------------------------------------------------',
% 
% 'PAL'
% pal,
% 
% '-----------------------------------------------------------------------',
% 
% 'USUAL AIC'
% aic,
% 
% '-----------------------------------------------------------------------',
%  'USUAL HANNAN QUINN'
%  hq,
% 
% '-----------------------------------------------------------------------',
result_all=[bic;pal;aic];
result_all=result_all/nsim
%save 'freqL2_AR_no_out_sig_1txt' freq_L2 -ASCII;
result_diff_n((count+1):(count+3),1)=numsamp;
result_diff_n((count+1):(count+3),2)=['B'; 'P'; 'A'];
result_diff_n((count+1):(count+3),3:7)=result_all;


count=count+3;
end% loop for numsamp
xlswrite('result_numsamp_50_100_sigmasq_15.xlsx', result_diff_n,'Sheet1');