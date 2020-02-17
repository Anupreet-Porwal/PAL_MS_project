function y=data_n_c(numsamp,numpar,a,om,sigsq);

for k=1:numsamp
    x(k)=0.0;
    for j=1:numpar
      x(k)=x(k)+a(j)*exp(i*2*pi*om(j)*k);
    end 
    y(k,1)=x(k)+(normrnd(0,sqrt(sigsq/2))+i*normrnd(0,sqrt(sigsq/2)));
end

    