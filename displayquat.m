
pts = [x,y];
QT = QuaTree(pts, 'binCapacity', 64, 'minSize', 0.05, 'style', 'equal');
figure(1);
boxH = QT.plot;
cols = lines(QT.BinCount);
doplot3 = @(p,varargin)plot(p(:,1),p(:,2),varargin{:});
for i = 1:QT.BinCount
%     set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
%     doplot3(pts(QT.PointBins==i,:),'.','Color',cols(i,:))
end
axis image, view(2)