% find the mode

load chaintest.mat
disp(['modes:'])
disp(['alpha = '         num2str(xplot(find(yplot(:,1) == max(yplot(:,1))),1))])
disp(['tau = '           num2str(xplot(find(yplot(:,2) == max(yplot(:,2))),2))])
disp(['m = '             num2str(xplot(find(yplot(:,3) == max(yplot(:,3))),3))])
disp(['S0 = '            num2str(xplot(find(yplot(:,4) == max(yplot(:,4))),4))])
disp([' '])

dx2 = diff(xplot,1).^2;
d2y = diff(diff(yplot,1),1);

uncertainties = -1./(d2y./dx2(1:(end-1),:));

disp(['uncertanties:'])
disp(['alpha = ' num2str(uncertainties(find(yplot(:,1) == max(yplot(:,1)))-1,1))])
disp(['tau = ' num2str(uncertainties(find(yplot(:,2) == max(yplot(:,2)))-1,2))])
disp(['m = ' num2str(uncertainties(find(yplot(:,3) == max(yplot(:,3)))-1,3))])
disp(['S0 = ' num2str(uncertainties(find(yplot(:,4) == max(yplot(:,4)))-1,4))])
disp([' '])
