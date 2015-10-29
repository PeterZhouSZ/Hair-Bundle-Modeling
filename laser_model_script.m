clf
t = (0:5e-5:0.4);
for p=linspace(0,5,5)
   [pr, nr, tr] = laser_model(t,0.1,0.14,p);
   subplot_tight(3,1,1,[0.03 0.1])
   plot(t,pr,'color',[0 0.8 0.2],'linewidth',2)
   set(gca,'XTickLabel',[])
   hold on
   axis([0 0.4 -150 210])
   subplot_tight(3,1,2,[0.03 0.1])
   plot(t,nr,'r','linewidth',2)
   set(gca,'XTickLabel',[])
   hold on
   axis([0 0.4 -150 210])
   subplot_tight(3,1,3,[0.03 0.1])
   plot(t,tr,'color',[0 0.3 0.9],'linewidth',2)
   hold on
   axis([0 0.4 -150 210])
end

