writerObj=VideoWriter('eps10_order4.avi','MPEG-4');
%writerObj.FrameRate=50;
writerObj.FrameRate=10;
writerObj.Quality=95;
open(writerObj);

xf=load('xf.dat');
yf=load('yf.dat');
pot=load('velocity.dat');
nleaves=load('nleaves.dat');



nt=size(nleaves,1); %time steps
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
        [m, n] = size(xt);
        plot3(xt(1,:), yt(1,:), p(1,:), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 顶边
        plot3(xt(m,:), yt(m,:), p(m,:), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 底边
        plot3(xt(:,1), yt(:,1), p(:,1), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 左边
        plot3(xt(:,n), yt(:,n), p(:,n), 'Color', [0, 0, 0, 1], 'LineWidth', 0.5); % 右边
  end
  view([0,90]);
  shading interp;
  colormap(jet(256));
  caxis([-25,25])
  colorbar
  hold off;
  nn=nn+nleaves_t;
  frame=getframe(gcf);
  writeVideo(writerObj,frame);
end


close(writerObj);
