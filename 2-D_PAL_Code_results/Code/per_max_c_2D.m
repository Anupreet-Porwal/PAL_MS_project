function inits=per_max_c_2D(yadj,numsamp,gridnum);
start_om=0;
end_om=1;
start_be=0;
end_be=1;
for k=1:gridnum+1
    for l=1:gridnum+1
        om_grid(k)=start_om+((end_om-start_om)/gridnum)*(k-1);
        be_grid(l)=start_be+((end_be-start_be)/gridnum)*(l-1);
    per(k,l)=0.0;
    for kk=1:numsamp
        for ll=1:numsamp
            per(k,l)=per(k,l)+yadj((ll-1)*numsamp+kk)*exp(-i*2*pi*(kk*om_grid(k)+ll*be_grid(l)));
        end
    end
    per(k,l)=per(k,l)/(numsamp*numsamp);
    per(k,l)=abs(per(k,l))^2;
    end
end
[maxper,ind] = max(per(:));
[m,n] = ind2sub(size(per),ind);

om=om_grid(m);
be=be_grid(n);

for kk=1:numsamp
    for ll=1:numsamp
        a_om_mat(kk,ll)=exp(i*2*pi*(om*kk+be*ll));
        %yvec((ll-1)*numsamp+kk,1)=yadj((ll-1)*numsamp+kk);
    end
end
yvec=yadj;
a_om=reshape(a_om_mat,[numsamp*numsamp, 1]);

lin=(a_om'*yvec)/(a_om'*a_om);
inits=[lin;om;be];