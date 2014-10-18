
numSamps = 10;  %values used in appendix 
mL = 380;
ML = 720;

%numSamps = 20;  %values used to generate curves in main body of the paper
%mL = 360;
%ML = 800;



%A is the matrix from spectra to RGB
%curves holds the basis functions used by the algorithm
%(white cyan magenta yellow red green blue.
[curves,A] = colorGen(mL, ML, numSamps);  
curves



%%%%%% Create the figures used in the paper

curves(numSamps+1,1) = 0;  

figure(1);

clf;
white = [1 1 1]';
An = null(A)
curves(1:numSamps,8) = .5 * pinv(A) *white;
curves(1:numSamps,9) = .5 * (nnls(A,white)); % - An(:,1) - 2*An(:,2) +An(:,4) - 2*An(:,5) + .4*An(:,6) + .4*An(:,7) - 2.5*An(:,8) + An(:,10) + An(:,11)+ An(:,12)+ An(:,13)+ An(:,14)+ An(:,15)+ An(:,16)+ An(:,17))
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[8 9]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
%xlabel('Wavelength');
set(h,'LineWidth',1,{'LineStyle'}, {'-';'-.'});
set(h,{'Color'}, {'k';'b'});
legend(h,'Grey1','Grey2',0);
axis([mL, ML, 0, 2]);
a = input('Continue?');
%%% %%% set(1,'PaperPosition',[.91,4.3,6.666,2.222]);
%print '/home/bes/local/docs/paper/color/greys.eps' -depsc;

curves(:,10) = curves(:,8) .* curves(:,8);
curves(:,11) = curves(:,8) .* curves(:,9);
curves(:,12) = curves(:,9) .* curves(:,9);
(A * curves(1:numSamps,9:12))'
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[10 11 12]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
%xlabel('Wavelength');
set(h,'LineWidth',1,{'LineStyle'}, {'-';'-.';'--'});
set(h,{'Color'}, {'k';'b';'r'});
legend(h,'Grey11','Grey12','Grey22',0);
axis([mL, ML, 0, 4]);
a = input('Continue?');
%%% %%% set(1,'PaperPosition',[.91,4.3,6.666,2.222]);
%print '/home/bes/local/docs/paper/color/product.eps' -depsc;


clf;
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[5 6 7]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
%xlabel('Wavelength');
set(h,'LineWidth',1,{'LineStyle'}, {'-.';'-';'--'});
set(h,{'Color'}, {'r';'g';'b'});
legend(h,'Red','Green','Blue',0);
axis([mL, ML, 0, 1.1]);
%a = input('Continue?');
%%% %%% set(1,'PaperPosition',[.91,4.3,6.666,2.222]);
%print '/home/bes/local/docs/paper/color/primaries.eps' -depsc;


clf;
sum = curves(:,1);
sum(:,2) = curves(:,5) + curves(:,6) + curves(:,7);
[a,b] = stairs(linspace(mL,ML,numSamps+1), sum);
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
set(h,'LineWidth',1,{'LineStyle'}, {'-';':'});
set(h,{'Color'}, {'k';'b'});
legend(h,'White','R+G+B',0);
axis([mL, ML, .8, 1.2]);
%xlabel('Wavelength');
%a = input('Continue?');
%%% %%% set(1,'PaperPosition',[.91,4.3,6.666,2.222]);
%print '/home/bes/local/docs/paper/color/sum.eps' -depsc;

clf;
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[2 3 4]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
set(h,'LineWidth',1,{'LineStyle'}, {'-.';':';'-'});
set(h,{'Color'}, {'c';'m';'y'});
legend(h,'Cyan', 'Magenta', 'Yellow',0);
axis([mL, ML, 0, 1.1]);
%xlabel('Wavelength');
%title('Spectra for cyan, magenta, and yellow.');
%a = input('Continue?');
%%% %%% set(1,'PaperPosition',[.91,4.3,6.666,2.222]);
%print '/home/bes/local/docs/paper/color/secondary.eps' -depsc;

clf;
subplot(1,2,1);
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[1 5 6 7]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
set(h,'LineWidth',1,{'LineStyle'}, {'-';'-.';'--';':'});
set(h,{'Color'}, {'k';'r';'g';'b'});
legend(h,'White','Red','Green','Blue',0);
axis([mL, ML, 0, 1.1]);
%xlabel('Wavelength');
%title('Spectra for white, red, green, and blue.');

subplot(1,2,2);
[a,b] = stairs(linspace(mL,ML,numSamps+1), curves(:,[2 3 4]));
h = plot(a(1:2*numSamps,:),b(1:2*numSamps,:));
set(h,'LineWidth',1,{'LineStyle'}, {'-.';':';'--'});
set(h,{'Color'}, {'c';'m';'y'});
legend(h,'Cyan', 'Magenta', 'Yellow',0);
axis([mL, ML, 0, 1.1]);
xlabel('Wavelength');
%title('Spectra for cyan, magenta, and yellow.');

%subplot(2,2,3);
%plot(A');


print 'tst.eps' -depsc;
