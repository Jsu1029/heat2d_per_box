function []=plottree(fileno)
filename='fort.';
filename=[filename,num2str(fileno)]
tree=load(filename);
[n,m]=size(tree);

for i=1:n
  rectangle('Position',tree(i,:));
  hold on;
end
axis equal;
xlim([-0.65,0.65]);
ylim([-0.65,0.65]);

end
