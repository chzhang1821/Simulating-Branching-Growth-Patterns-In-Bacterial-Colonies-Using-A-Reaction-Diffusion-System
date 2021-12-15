% Bacteria diff eq: B' = B + (Db*del2(B) + B*N*heaviside(B-beta) - mu*B)*dt
% Nutrient diff eq: N' = N + (del2(B) - B*N*heaviside(B-beta))*dt
% "frozen" bacterial diff eq: S' = S + mu*B*dt

clear all;
close all;

% defining the mesh
width= 500; %size of grid in each dimension x & y
h = 1/width; % step size in x and y


% initial condition
n0 = 1; % initial concentration of nutrient
b0 = 1; % initial concentration of bacteria


% parameters
Db = 0.01;
beta = 0.25;
mu = 0.01;


% dt = .001*h*h; % time step - usually small
dt = 0.25;
stoptime = 5000;
count_t=0;

[t,B,N,S] = initial_conditions(width,n0,b0);

tic
nframes = 1;
while t<stoptime
    % Euler algorithm
    nnew = N + (my_laplacian(N,width) - B.*N.*heaviside(B-beta))*dt;
    bnew = B + (Db*my_laplacian(B,width) + B.*N.*heaviside(B-beta) - mu*B)*dt;
    snew = S + mu*B*dt;
    N = nnew;
    B = bnew;
    S = snew;

 if mod(count_t,100)==0
    t
    nframes
    figure(1)
    h1=heatmap(N,'GridVisible','off');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    title('The density of nutrient')
    caxis([0,1])
    filename1 = 'nutrient.gif';
    
    drawnow;
    F1 = getframe(gcf);
    Im1 = frame2im(F1);
    [A1,map1] = rgb2ind(Im1,256);

    if nframes == 1
        imwrite(A1,map1,filename1,'gif','Loopcount',Inf,'DelayTime',0.05);
    else 
        imwrite(A1,map1,filename1,'gif','WriteMode','append','DelayTime',0.05);
    end

    figure(2)
    h2=heatmap(B,'GridVisible','off');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    title('The density of active bacteria')
    filename2 = 'active_bacteria.gif';

    drawnow;
    F2 = getframe(gcf);
    Im2 = frame2im(F2);
    [A2,map2] = rgb2ind(Im2,256);
    if nframes == 1
        imwrite(A2,map2,filename2,'gif','Loopcount',Inf,'DelayTime',0.05);
    else 
        imwrite(A2,map2,filename2,'gif','WriteMode','append','DelayTime',0.05);
    end

    figure(3)
    h3=heatmap(S,'GridVisible','off');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    title('The density of "frozen" bacteria')
    filename3 = 'frozen_bacteria.gif';

    drawnow;
    F3 = getframe(gcf);
    Im3 = frame2im(F3);
    [A3,map3] = rgb2ind(Im3,256);
    if nframes == 1
        imwrite(A3,map3,filename3,'gif','Loopcount',Inf,'DelayTime',0.05);
    else 
        imwrite(A3,map3,filename3,'gif','WriteMode','append','DelayTime',0.05);
    end
%     pause(0.0001);
end
    t = t+dt;
    count_t=count_t+1;
    nframes = nframes+1;
end


% % plot image
% figure(1)
% heatmap(N,'GridVisible','off');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% caxis([0,1]);
% title('The density of nutrient')
% 
% figure(2)
% heatmap(B,'GridVisible','off');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% title('The density of active bacteria')
% 
% figure(3)
% heatmap(S,'GridVisible','off');
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% title('The density of "frozen" bacteria')

delta = toc;
disp([num2str(nframes) ' frames in ' num2str(delta) ' seconds'])


%% functon of initial conditions
function [t,B,N,S] = initial_conditions(n,b0,n0)
    t = 0;

    % one colony
    m = n/2;

    % Initialize B to zero with a clump of b0
    B = zeros(n);

    % Location of each bacteria colony
    % Location of one colony
    [x,y] = meshgrid((1:n)-m,(1:n)-m);
    C = sqrt(x.^2+y.^2)<=15;
    B(C) = b0;
    for i = 1:n
        for j = 1:n
            B(i,j) = B(i,j)-B(i,j)*rand(1,1);
        end
    end

    % Locations of two colonies
%     m1 = n*4/12; % Location of colony 1
%     m2 = n*8/12; % Location of colony 2
% 
%     [x1,y1] = meshgrid((1:n)-m1,(1:n)-m);
%     C1 = sqrt(x1.^2+y1.^2)<=10;
%     [x2,y2] = meshgrid((1:n)-m2,(1:n)-m);
%     C2 = sqrt(x2.^2+y2.^2)<=10;
%     B(C1) = b0;
%     B(C2) = b0;
%     for i = 1:n
%         for j = 1:n
%             B(i,j) =B(i,j)-B(i,j)*rand(1,1);
%         end
%     end

    % Initialize N to n0
    N = ones(n);
    for i = 1:n
        for j = 1:n
            N(i,j) = n0+0.05*rand(1,1)*(-1)^randi(2);
        end 
    end

    % Initialize S to zero
    S = zeros(n);
end


%% laplacian function
% function out = my_laplacian(in,N)
%   out = 0.25*(in(:,[2:N+1 N])+in(:,[2 1:N])+in([2 1:N],:)+in([2:N+1 N],:))...
%       + 0.5*(in([2:N+1 N],[2:N+1 N])+in([2:N+1 N],[2 1:N])...
%       + in([2 1:N],[2 1:N])+in([2 1:N],[2:N+1 N]))...
%       - 3*in;
% end


% function out = my_laplacian(in)
%   out = -in ...
%       + .20*(circshift(in,[ 1, 0]) + circshift(in,[-1, 0])  ...
%       +      circshift(in,[ 0, 1]) + circshift(in,[ 0,-1])) ...
%       + .05*(circshift(in,[ 1, 1]) + circshift(in,[-1, 1])  ...
%       +      circshift(in,[-1,-1]) + circshift(in,[ 1,-1]));
% end

% function out = my_laplacian(in,N)
%   out = in(:,[2:N+1 N])+in(:,[2 1:N])+in([2 1:N],:)+in([2:N+1 N],:)...
%       - 4*in;
% end

function out = my_laplacian(in,N)
    out = zeros(N);
    for x = 2:N-1
        for y = 2:N-1
            out(x,y) = -3*in(x,y)...
                + 0.5*(in(x+1,y)+in(x-1,y)+in(x,y+1)+in(x,y-1))...
                + 0.25*(in(x+1,y+1)+in(x-1,y+1)+in(x+1,y-1)+in(x-1,y-1)); 
        end
    end
    out(:,N) = out(:,N-1);
    out(:,1) = out(:,2);
    out(N,:) = out(N-1,:);
    out(1,:) = out(2,:);
end




