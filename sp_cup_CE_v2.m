clc;
clear all;

load dataset2.mat

P=pilotMatrix.';
for user_no=1:50;
    user_no
h_eff_est=zeros(4096,500);
for k=1:500;
    
y=receivedSignal(k,:,user_no);
y_tilde=receivedSignal(k,:,user_no)/transmitSignal(k);
h_eff_est(:,k)=(P.')*y_tilde.';
h_eff_final(:,k,user_no)=h_eff_est(:,k);
end

end

save('h_eff_final','h_eff_final')