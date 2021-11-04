
function [Podoconnect, xpod, ypod, Npod] = PodoConnectivity_hexagon(Na, d0)


Npod = 3*Na*(Na-1)+1; % total podosome number
xpod = zeros(Npod,1);
ypod = zeros(Npod,1);
Podoconnect = cell(Npod,1);

% define the coordinates
Angle_sides = [pi/3, 0, -pi/3, -2*pi/3, pi, 2*pi/3];
for i = 2:Na   % hexagon boundary number
    Nfirst = 3*(i-1)*(i-2)+2;
    xpod(Nfirst) = -d0*(i-1);
    ypod(Nfirst) = 0;
    
    for j=1:6*(i-1)-1  % six sides of hexagon
        Ang = Angle_sides(floor((j-1)/(i-1))+1);
        xpod(Nfirst+j) = xpod(Nfirst+j-1)+d0*cos(Ang);
        ypod(Nfirst+j) = ypod(Nfirst+j-1)+d0*sin(Ang);
    end
end

% define the connectivity, note here they can be overlapped
for i = 1:Npod
    neighbors = [];
    for j = 1:Npod
        if abs((xpod(i)-xpod(j))^2+(ypod(i)-ypod(j))^2-d0^2)<0.1*d0
            neighbors = [neighbors, j];
        end
    end
    [~, orderY] = sort(ypod(neighbors));
    if (length(neighbors)==6)&&(xpod(neighbors(orderY(1)))-xpod(neighbors(orderY(5)))~=0)
        nn = orderY(1);
        orderY(1) = orderY(2);
        orderY(2) = nn;
    end
    Podoconnect{i} = neighbors(orderY);
end

% % Plot the results
% figure(1)
% hold on
% for istart = 1:Npod
%     neighbors = Podoconnect{istart};
%     for j = 1:length(neighbors)
%         iend = neighbors(j);
%         plot([xpod(istart), xpod(iend)],[ypod(istart), ypod(iend)],'k:','linewidth', 1)
%     end
% end
% scatter(xpod, ypod);
% grid on;
% axis equal

end

