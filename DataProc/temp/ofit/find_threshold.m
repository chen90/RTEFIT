function threshold=find_threshold(img,showfig)
img=medfilt2(img);
xmax=max(max(img));
[counts,temp]=imhist(img/xmax,255);
x=temp*xmax;
x=x(2:end);
counts=counts(2:end);

ft1 = fittype( 'gauss2' );
opts_g = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_g.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts_g.StartPoint = [234670 76.0629921259843 118.145967089219 4722.74299652904 380.314960629921 239.486101175188];
[fitresult_g, ~] = fit(x, counts, ft1, opts_g );
b1=fitresult_g.b1;
b2=fitresult_g.b2;
center=max(b1,b2);
g_fit=feval(fitresult_g,x);
if b1<0 || b2<0
   center=2*max(b1,b2); 
end

x2=x(x>center);
counts2=counts(x>center);
%---------polyfit of the light erea---------------
ft2 = fittype( 'poly1' );
opts_p = fitoptions( 'Method', 'LinearLeastSquares' );
opts_p.Normalize = 'on';
opts_p.Robust = 'Bisquare';
[fitresult_p, ~] = fit(x2, counts2, ft2, opts_p);
p_fit = feval(fitresult_p,x2);
%--------Exponential fit of the background--------------
ft3 = fittype( 'exp1' );
opts_e = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_e.Normalize = 'on';
opts_e.StartPoint = [680.917268177929 -1.4450527287455];
[fitresult_e, ~] = fit( x2, counts2, ft3, opts_e );
e_fit = feval(fitresult_e,x2);
%%
new=abs(e_fit-p_fit);
new=new(1:100);
threshold=x2(new==min(new));
%%
if showfig
figure('unit','normalized','DefaultAxesFontSize',14,'DefaultAxesFontWeight','normal','DefaultAxesLineWidth',1.5,'position',[0.1,0.1,0.6,0.6]);
plot(x,counts,':o','MarkerSize',5,'Linewidth',2,'MarkerFace',[0 0.45 0.74],'MarkerFaceColor',[0 0.45 0.74]);
hold on;plot(x,g_fit,'k',x2,e_fit,'r',x2,p_fit,'g','Linewidth',2.5);
legend('Hist','Double Gaussian Fitting','Linear Fitting','Exponential fitting')
hold on;plot([b1,b1],[0,1e5],'--k','linewidth',1)
hold on;plot([b2,b2],[0,1e5],'--k','linewidth',1)
xlim([0,2000])
ylim([0,0.6e5])
legend('Hist','Double Gaussian Fitting','Linear Fitting','Exponential fitting')
end

end