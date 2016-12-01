%% �����½���ͼʾ
% ���ò���Ϊ0.1��f_changeΪ�ı�ǰ���yֵ�仯����������һ���˳�������
syms x;f=x^2;
step=0.1;x=2;k=0;         %���ò���,��ʼֵ,������¼��
f_change=x^2;             %��ʼ����ֵ
f_current=x^2;            %���㵱ǰ����ֵ
ezplot(@(x,f)f-x^2)       %��������ͼ��
axis([-2,2,-0.2,3])       %�̶�������
hold on
while f_change>0.000000001                %�������������μ����ֵ֮��С��ĳ����������ѭ��
    x=x-step*2*x;                         %-2*xΪ�ݶȷ�����stepΪ�������������½�����
    f_change = f_current - x^2;           %�������κ���ֵ֮��
    f_current = x^2 ;                     %���¼��㵱ǰ�ĺ���ֵ
    plot(x,f_current,'ro','markersize',7) %��ǵ�ǰ��λ��
    drawnow;pause(0.2);
    k=k+1;
end
hold off
fprintf('�ڵ���%d�κ��ҵ�������СֵΪ%e����Ӧ��xֵΪ%e\n',k,x^2,x)

%% contour for gradient
[x,y] = meshgrid(0:10, 0:10);
z = x.^2 + y.^2;
contour(z)