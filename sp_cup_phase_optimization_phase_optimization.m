clc
clear all
load('h_eff_final.mat')
SNR_db=60; %SNR per subcarrier
SNR=10^(SNR_db/10);  
B=10*10^(6);
K=500;
M=20;
sub_carriers=500;
for user_no=1:50;
    user_no
%phi_1=sign(real(exp(1i*randn(4096,1)*pi))); %initialization of phases
phi_1=ones(4096,1); %initialization of phases

R_achieved_k_1(user_no)=0;
R_achieved_k(user_no)=0;
sum_1=0;
while abs(R_achieved_k_1(user_no)-R_achieved_k(user_no))>=0.001*R_achieved_k_1(user_no);
    R_achieved_k_1(user_no)=R_achieved_k(user_no);
    sum_1=0;
for i=1:sub_carriers

H_tilde=SNR*h_eff_final(:,i,user_no)*h_eff_final(:,i,user_no)'*phi_1;
g=2*real(H_tilde);
c=1-phi_1'*SNR*h_eff_final(:,i,user_no)*h_eff_final(:,i,user_no)'*phi_1+phi_1'*2*real(SNR*h_eff_final(:,i,user_no)*h_eff_final(:,i,user_no)'*phi_1);
term_1=g/c;
c_k=1-phi_1'*SNR*h_eff_final(:,i,user_no)*h_eff_final(:,i,user_no)'*phi_1-2*sum(abs(g));
if c_k>0;
term_2=2*phi_1*norm(g)^2/(c_k)^3*c;
else 
    term_2=0;
end
sum_1=sum_1+term_1;%+term_2;

end
phi=sign(real(sum_1));
phi_1=phi;

R_achieved_k(user_no)=0;
for i=1:sub_carriers;
   
    R_achieved_k(user_no)=R_achieved_k(user_no)+B/(K+M-1)*log2(1+phi'*SNR*h_eff_final(:,i,user_no)*h_eff_final(:,i,user_no)'*phi);
end

end
phi_optimal(:,user_no)=phi;



R_upper_bound_1(user_no)=0;
for i=1:sub_carriers
   
    
    
    R_upper_bound_1(user_no)=R_upper_bound_1(user_no)+B/(K+M-1)*log2(1+(sqrt(4096)*norm(sqrt(SNR)*h_eff_final(:,i,user_no)))^2);
end

R_upper_bound_2(user_no)=0;
for i=1:sub_carriers
   
    
    
    R_upper_bound_2(user_no)=R_upper_bound_2(user_no)+B/(K+M-1)*log2(1+(norm(sqrt(SNR)*h_eff_final(:,i,user_no),1))^2);
end
end

save('phase_configuration','phi_optimal')
save('Rate_achieved','R_achieved_k')
save('Rate_upperbound','R_upper_bound_1','R_upper_bound_2')

