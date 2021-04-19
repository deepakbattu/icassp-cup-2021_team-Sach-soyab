clc
clear all
load('h_eff_final.mat')
user_no=2;
h_eff_final_1=h_eff_final(:,:,user_no);
%phi_1=sign(real(exp(1i*randn(4096,1)*pi))); %initialization of phases
phi_1=ones(4096,1);
SNR_db=60; %SNR per subcarrier
SNR=10^(SNR_db/10);  
B=10*10^(6);
K=500;
M=20;
sub_carriers=500;
sum_1=0;
for iter=1:3
    iter
    sum_1=0;
for i=1:sub_carriers
    i
H=SNR*h_eff_final_1(:,i)*h_eff_final_1(:,i)';
H_tilde=SNR*h_eff_final_1(:,i)*h_eff_final_1(:,i)'*phi_1;
g=2*real(H_tilde);
c=1-phi_1'*H*phi_1+phi_1'*2*real(H*phi_1);
metric_1=g/c;
c_k=1-phi_1'*H*phi_1-2*sum(abs(g));
if c_k>0;
metric_2=2*phi_1*norm(g)^2/(1-phi_1'*H*phi_1-2*sum(abs(g)))^3*c;
else 
    metric_2=0;
end
sum_1=sum_1+metric_1+metric_2;%2*H_tilde(:,i)/(c(:,i));

end
phi=sign(real(sum_1));

phi_1=phi;

R_achieved(iter,1)=0;
for i=1:sub_carriers;
    metric_3(:,i)=phi'*SNR*h_eff_final_1(:,i)*h_eff_final_1(:,i)'*phi;
    R_achieved(iter,1)=R_achieved(iter,1)+B/(K+M-1)*log2(1+phi'*SNR*h_eff_final_1(:,i)*h_eff_final_1(:,i)'*phi);
end
R_achieved(iter,1)
end
R_upper_bound=0;
for i=1:sub_carriers
   
    
    
    R_upper_bound=R_upper_bound+B/(K+M-1)*log2(1+(sqrt(4096)*norm(sqrt(SNR)*h_eff_final_1(:,i)))^2);
end
R_upper_bound
