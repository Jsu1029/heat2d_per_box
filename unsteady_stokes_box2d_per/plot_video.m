writerObj=VideoWriter('us_test3_moving2.avi','MPEG-4');
%  writerObj.FrameRate=50;
writerObj.FrameRate=10;
writerObj.Quality=95;
open(writerObj);
% 
  xf=load('xf.dat');
  yf=load('yf.dat');
  pot=load('velocity.dat');
  nleaves=load('nleaves.dat');



nt=100; %time steps
nn=0;

for i=1:nt
  nleaves_t=nleaves(i);
  i
  for j=1:nleaves_t
      xt=zeros(8,8);
      yt=zeros(8,8);
      p=zeros(8,8);
      xt(:,:)=xf((nn+j-1)*8+1:(nn+j)*8,1:8);
      yt(:,:)=yf((nn+j-1)*8+1:(nn+j)*8,1:8);
      p(:,:)=pot((nn+j-1)*8+1:(nn+j)*8,1:8);
      surf(xt,yt,p);
      hold on
      % 手动绘制边界线为黑色
        [m, n] = size(xt);
        % 绘制 xt, yt 边界上的线条
        % 设置透明度为1，调整线条粗细为0.5
        plot3(xt(1,:), yt(1,:), p(1,:), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 顶边
        plot3(xt(m,:), yt(m,:), p(m,:), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 底边
        plot3(xt(:,1), yt(:,1), p(:,1), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 左边
        plot3(xt(:,n), yt(:,n), p(:,n), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 右边
  end
  
  view([0,90]);
  shading interp;
  colormap(jet(256));
%  caxis([-60,60])
  colorbar
  hold off;
  nn=nn+nleaves_t;
  frame=getframe(gcf);
  writeVideo(writerObj,frame);
end


close(writerObj);
