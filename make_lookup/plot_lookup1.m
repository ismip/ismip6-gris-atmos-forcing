figure
hold on; box on;
time = t;
look1 = look;
for b=1:25
    subplot(5,5,b)        
    look1(2,:) = table(b,:,time);
    plot(look1(1,:),look1(2,:)*31556926/1000)
    title(['B' num2str(b) ])
    axis([0 3500 -5 1])
end
xlabel('Elevation (m)')
ylabel('DSMB (m yr-1)')

