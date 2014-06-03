% find the mode

load chaintest.mat
disp([' '])
disp(['alpha_mod = '         num2str(xplot(find(yplot(:,1) == max(yplot(:,1))),1))])
disp(['tau_mod = '           num2str(xplot(find(yplot(:,2) == max(yplot(:,2))),2))])
disp(['gamma_mod = '         num2str(xplot(find(yplot(:,3) == max(yplot(:,3))),3))])
disp(['m_mod = '             num2str(xplot(find(yplot(:,4) == max(yplot(:,4))),4))])
disp(['S0_mod = '            num2str(xplot(find(yplot(:,5) == max(yplot(:,5))),5))])
disp(['E0 = '                num2str(xplot(find(yplot(:,6) == max(yplot(:,6))),6))])