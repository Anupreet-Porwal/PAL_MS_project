function inits=per_max_c(yadj,numsamp,gridnum);
start_om=0;
end_om=1;
for k=1:gridnum+1
    om_grid(k)=start_om+((end_om-start_om)/gridnum)*(k-1);
    per(k)=0.0;
    for kk=1:numsamp
        per(k)=per(k)+yadj(kk)*exp(-i*kk*2*pi*om_grid(k));
    end
    per(k)=per(k)/numsamp;
    per(k)=abs(per(k))^2;
end
[Y,I]=max(per);
om=om_grid(I);

for kk=1:numsamp
    a_om(kk,1)=exp(i*2*pi*om*kk);
    yvec(kk,1)=yadj(kk);
end
lin=(a_om'*yvec)/(a_om'*a_om);
inits=[lin;om];