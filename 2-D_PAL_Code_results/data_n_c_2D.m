function y=data_n_c_2D(numsamp,numpar,a,freqd,sigsq);

for k=1:numsamp
    for l =1:numsamp
        x(k,l)=0.0;
        for j=1:numpar
            x(k,l)=x(k,l)+a(j)*exp(i*2*pi*freqd(1,j)*k+i*2*pi*freqd(2,j)*l);
        end 
    y_mat(k,l)=x(k,l)+(normrnd(0,sqrt(sigsq/2))+i*normrnd(0,sqrt(sigsq/2)));
    end
end
y=reshape(y_mat,[numsamp*numsamp,1]);

    


    