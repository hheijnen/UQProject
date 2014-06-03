% find the mode

load chaintest.mat
disp([' '])
disp(['alpha_mod = '         num2str(xplot(find(yplot(:,1) == max(yplot(:,1))),1))])
disp(['tau_mod = '           num2str(xplot(find(yplot(:,2) == max(yplot(:,2))),2))])
disp(['m_mod = '             num2str(xplot(find(yplot(:,3) == max(yplot(:,3))),3))])
disp(['S0_mod = '            num2str(xplot(find(yplot(:,4) == max(yplot(:,4))),4))])
