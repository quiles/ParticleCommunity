clear all;


    a = load('./position0.txt');
    [N C] = size(a);

    cores = hsv(max(a(:,4)));
    
    for j=1:N
        p = plot3(a(j,1),a(j,2),a(j,3),'.');
        hold on;
        set(p,'Color',cores(a(j,4),:), 'MarkerSize',20);

    end;
    box on;
    hold off;



