function [y1,ysel,x1,xsel] = find_local_average(x,y,x0,dx)
% calculate local average y for a given x range

eid=find(abs(x-x0)<dx);
if(~isempty(eid))
    ysel=y(eid);
    xsel=x(eid);
    y1=nanmedian(ysel);
    x1=nanmedian(xsel);
else
    ysel=NaN;
    xsel=NaN;
    y1=NaN;
    x1=NaN;
end
