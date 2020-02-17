
% MODEL SELECTION OF CISOIDS WITH PAL AND BOOTSTRAPING
% April 03, 2017
% CISOIDS MODEL  SIMULATIONS 
clear;
clc;
%a(1)=7+i*3;
%a(2)=5+4*i;
%a(3)=3+2*i;

a(1)=3+i*2;
a(2)=2+1.66*i;
a(3)=1.75+i;
om(1)=0.4; om(2)=0.6; 
om(3)=0.7;
% a(1)=3+i*2;
% a(2)=0.2+0.16*i;
% a(3)=1.75+i;
% om(1)=0.42; om(2)=0.45; om(3)=0.7;% 5 4 
freqd=om;%*pi;
sz=size(a);
numpar=sz(1,2);
count=0;
maxcomp=10; % maximum number of components
% for numsamp =100;%:25:75;
for numsamp =50;
% 'Sample Size',numsamp
sigsq=30;
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
y=data_n_c(numsamp,numpar,a,om,sigsq);      

%independent t error
%y=data_t(numsamp,a,b,om,df,numpar);      

% normal mixture error
%y=data_n_mix(numsamp,numpar,a,b,om,stdev1,stdev2); 
% 
% %Estimation of sigma from full model
% yadj=y;
% 
% gridnum=1000; % number of grid points for periodogram maximizer 
% 
% for ks=1:maxcomp
%    inits=per_max_c(yadj,numsamp,gridnum); 
%    seq_a(ks)=inits(1);
%    seq_om(ks)=inits(2);
% for ns=1:numsamp
%     yadj(ns)=(yadj(ns)-inits(1)*exp(i*2*pi*inits(2)*ns));
% end
% 
% end
% 
% initial=[seq_a.',seq_om'];
% L2_full=fminsearch('obj_L2_C',initial,[],y,numsamp,comp);
% 
% sigsq_null=0;
% sigsq_full=0;
% for k=1:numsamp
%     yp(k)=0.0;
%     for kp=1:maxcomp
%         yp(k)=yp(k)+(L2_full(kp,1)*exp(i*2*pi*L2_full(kp,2)*k));
%     end
% sigsq_null=sigsq_null+abs(y(k))^2;    
% diff_full_sq(k)=abs(y(k)-yp(k))^2;
% sigsq_full=sigsq_full+diff_full_sq(k);
% end
% sigsq_null=sigsq_null/(numsamp);
% sigsq_full=sigsq_full/(numsamp);
% 
% 
% for comp=1:maxcomp
%     nummax=comp;
% yadj=y;
% 
% gridnum=1000; % number of grid points for periodogram maximizer 
% 
% for ks=1:nummax
%    inits=per_max_c(yadj,numsamp,gridnum);
%    sequ_a(ks)=inits(1);
%    sequ_om(ks)=inits(2);
% for ns=1:numsamp
%     yadj(ns)=(yadj(ns)-inits(1)*exp(i*2*pi*inits(2)*ns));
% end
% end
% 
% initial_m=[sequ_a.',sequ_om'];
% L2_model=fminsearch('obj_L2_C',initial_m,[],y,numsamp,comp);
% 
% sum_sq=0.0;
% for k=1:numsamp
%     yp(k)=0.0;
%     for kp=1:comp
%         yp(k)=yp(k)+(L2_model(kp,1)*exp(i*2*pi*L2_model(kp,2)*k));
%     end
%         diff(k)=(y(k)-yp(k));
%         sum_sq=sum_sq+abs(diff(k))^2;
% 
% end
% est_var=sum_sq/(numsamp);
% pal_est_sigsq(comp)=est_var;
% 
% 
% if comp==1
%     pal_rn=2*numsamp*log(sigsq_null/sigsq_null);
%     pal_rhon=2*numsamp*log(sigsq_null/sigsq_full);
% else
%     pal_rn=2*numsamp*log(sigsq_null/pal_est_sigsq(comp-1));
%     pal_rhon=2*numsamp*log(pal_est_sigsq(comp-1)/sigsq_full);
% end
% 
% %bic_ic_obj(comp)=2*numsamp*log(est_var)+(log(numsamp))*(3*comp+1);
% %bic_cor_obj(comp)=2*numsamp*log(est_var)+(log(numsamp))*(5*comp+1);
% pal_obj(comp)=2*numsamp*log(est_var)+(3*comp+1)*log(3*maxcomp+1)*(log((pal_rn)+1)/log((pal_rhon)+1));
% %aic_obj(comp)=2*numsamp*log(est_var)+2*(3*comp+1);
% %aicc_obj(comp)=numsamp*log(est_var)+2*(3*comp*numsamp)/(numsamp-3*comp-1);
% %ca=2*(3*comp*numsamp)/(numsamp-3*comp-1); cb=(log(numsamp))*(5*comp);  % check the expression
% %wic_obj(comp)=((ca)/(ca+cb))*aicc_obj(comp)+((cb)/(ca+cb))*bic_cor_obj(comp);
% %hq_obj(comp)=2*numsamp*log(est_var)+2*(3*comp+1)*log(log(numsamp));
% 
% end % loop for all possible number of components
% % PAL
% [CR,I]=min(pal_obj);
% estcomp=I;
%'PAL', estcomp
estcomp=palcalculator(y,maxcomp,numsamp);
nbootsim=100;

    
    nummax=estcomp;
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
    L2_model=fminsearch('obj_L2_C',initial_m,[],y,numsamp,nummax);
%sum_sq=0.0;
    for k=1:numsamp
        yp(k)=0.0;
        for kp=1:nummax
            yp(k)=yp(k)+(L2_model(kp,1)*exp(i*2*pi*L2_model(kp,2)*k));
        end
        diff(k)=(y(k)-yp(k));
%        sum_sq=sum_sq+abs(diff(k))^2;
    end
for bootsim=1:nbootsim    
    diffboot(bootsim,:)=datasample(diff,numsamp);
    bootsamp(bootsim,:)=yp+diffboot(bootsim,:);
    pal(isim,bootsim)=palcalculator(bootsamp(bootsim,:),maxcomp,numsamp);
end
%     if estcomp == 1
%         pal(1)=pal(1)+1;
%     elseif estcomp == 2
%         pal(2)=pal(2)+1;
%     elseif estcomp == 3
%         pal(3)=pal(3)+1;
%     elseif estcomp == 4
%         pal(4)=pal(4)+1;
%     elseif estcomp == 5
%         pal(5)=pal(5)+1;
%     elseif estcomp == 6
%         pal(6)=pal(6)+1;
%     elseif estcomp == 7
%         pal(7)=pal(7)+1;
%     elseif estcomp == 8
%         pal(8)=pal(8)+1;
%     elseif estcomp == 9
%         pal(9)=pal(9)+1;
%     else estcomp == 10
%         pal(10)=pal(10)+1;
%     end
% 
palquantiles(isim,:) = quantile(pal(isim,:),[0.025 0.05 0.50 0.95 0.975])
% USUAL BIC
% [CR,I]=min(bic_ic_obj);
% estcomp=I;
% %'USUAL BIC', estcomp
% 
%     if estcomp == 1
%         bic(1)=bic(1)+1;
%     elseif estcomp == 2
%         bic(2)=bic(2)+1;
%     elseif estcomp == 3
%         bic(3)=bic(3)+1;
%     elseif estcomp == 4
%         bic(4)=bic(4)+1;
%     elseif estcomp == 5
%         bic(5)=bic(5)+1;
%     elseif estcomp == 6
%         bic(6)=bic(6)+1;
%     elseif estcomp == 7
%         bic(7)=bic(7)+1;
%     elseif estcomp == 8
%         bic(8)=bic(8)+1;
%     elseif estcomp == 9
%         bic(9)=bic(9)+1;
%     else estcomp == 10
%         bic(10)=bic(10)+1;
%     end
% 
% % CORRECTED BIC/MAP
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



% USUAL AIC
% [CR,I]=min(aic_obj);
% estcomp=I;
% %'USUAL AIC', estcomp
% 
% 
%     if estcomp == 1
%         aic(1)=aic(1)+1;
%     elseif estcomp == 2
%         aic(2)=aic(2)+1;
%     elseif estcomp == 3
%         aic(3)=aic(3)+1;
%     elseif estcomp == 4
%         aic(4)=aic(4)+1;
%     elseif estcomp == 5
%         aic(5)=aic(5)+1;
%     elseif estcomp == 6
%         aic(6)=aic(6)+1;
%     elseif estcomp == 7
%         aic(7)=aic(7)+1;
%     elseif estcomp == 8
%         aic(8)=aic(8)+1;
%     elseif estcomp == 9
%         aic(9)=aic(9)+1;
%     else estcomp == 10
%         aic(10)=aic(10)+1;
%     end
% 

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

%result_all=[pal]
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
% result_all=[bic;bicc;pal;aic];
% result_all=result_all/nsim
% %save 'freqL2_AR_no_out_sig_1txt' freq_L2 -ASCII;
% result_diff_n((count+1):(count+4),1)=numsamp;
% result_diff_n((count+1):(count+4),2)=['B'; 'c'; 'P'; 'A'];
% result_diff_n((count+1):(count+4),3:12)=result_all;


%count=count+4;
end% loop for numsamp
xlswrite('result_numsamp_50_sigmasq_30_bootstrap_quantiles.xlsx', palquantiles,'Sheet1');
xlswrite('result_numsamp_50_sigmasq_30_bootstrap_distributions.xlsx', pal,'Sheet1');
%xlswrite('result_numsamp_5_15_25_sigmasq_0_5.xlsx', result_diff_n,'Sheet1');