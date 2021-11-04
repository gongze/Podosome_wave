function movieplot2(Na, d0, zpodt)
h = figure();
[Podoconnect, xpod, ypod, Npod] = PodoConnectivity_hexagon(Na, d0);
scatter(xpod, ypod);
grid on;
axis equal
%axis tight manual
ax = gca;
set(ax,'visible','off')
ax.NextPlot = 'replaceChildren';


loops = size(zpodt,1);
M(loops) = struct('cdata',[],'colormap',[]);
% length = 5;
Dblue = [0, 0, 1];
Lblue = [200, 200, 255]/255;
len_color = size(jet, 1);
% colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)',...
%linspace(red(3),pink(3),length)'];
% colormap(colors_p)
%zpod = xpod*0.1+0.1;
zpodmax = 850; %max(max(zpodt,[],1));
zpodmin = 350; %min(min(zpodt,[],1));
h.Visible = 'off';
video_podo = VideoWriter('Podo_collect_center.avi');
video_podo.FrameRate = 20; video_podo.Quality = 80;
open(video_podo);
for j = 1:loops
    zpod = zpodt(j,:)';
    sizepod = (zpod-zpodmin)*170/(zpodmax-zpodmin)+ones(size(xpod,1),1)*5;
    sizepod(sizepod<2)=2; sizepod(sizepod>180)=180;
%     colorpod = (zpod-zpodmin)*(Dblue-Lblue)/(zpodmax-zpodmin)+ones(size(xpod,1),1)*Lblue;
    colorpodnum = (zpod-zpodmin)*(len_color-1)/(zpodmax-zpodmin)+ones(size(xpod,1),1);
    colorpod = interp1(1:len_color, jet, colorpodnum);
    scatter(xpod, ypod, sizepod, colorpod, 'filled');
    colormap('jet'); colorbar 
    caxis([300 900])
    %axis([-0.3*d0, (Na+0.1)*d0, -0.3*d0, (Na-0.5)*d0])
    %surf(X,Z)  % Z is the color matrix to fix the color at certain position
    drawnow
    M(j) = getframe(gcf);
    writeVideo(video_podo, M(j))
end

close(video_podo);
h.Visible = 'on';

%movie(M,1,3);
%movie2avi(M, 'Podo_collect_rand3.avi', 'compression', 'None', 'fps', 20);
end

% figure(1)
% hold on
% for istart = 1:Na^2
%     neighbors = Podoconnect{istart};
%     for j = 1:length(neighbors)
%         iend = neighbors(j);
%         plot([xpod(istart), xpod(iend)],[ypod(istart), ypod(iend)],'k:','linewidth', 1)
%     end
% end