function [xsorted,ysorted,indexsorted,xtopsorted,ytopsorted,indextop,xbottomsorted,ybottomsorted,indexbottom] = findbijectiveregion(x,y)
% interp1 fails when we have a multivalued function. Here, we assume that
% X has a maximal value XM corresponding to unique Y value YM.
% Below XM, Y follows two paths so that for any X, we get Y1 and Y2 with
% Y1(X) < YM and Y2(X) > YM.
% Given a vector of data points x and y, this function first sorts the
% points into vectors xsorted and ysorted and gives their indices from the original xy data
% indexsorted. It also outputs the same data for only the upper or lower
% portions of the curve.
% Interp1 can then be applied on either the top or bottom sorted data.
% (c) Copyright Dan Byrom and Alexander Darlington 2024.

[~,xmindex] = max(x);
ym = y(xmindex);

xtop = x(y>=ym);
ytop = y(y>=ym);
xbottom = x(y<ym);
ybottom = y(y<ym);

[xtopsorted,Itop] = sort(xtop);
ytopsorted = ytop(Itop);
[~,indextop] = ismember(xtopsorted,x);

[xbottomsorted,Ibottom] = sort(xbottom,'descend');
ybottomsorted = ybottom(Ibottom);
[~,indexbottom] = ismember(xbottomsorted,x);

xsorted = [xtopsorted;xbottomsorted];
ysorted = [ytopsorted;ybottomsorted];
[~,indexsorted] = ismember(xsorted,x);
end