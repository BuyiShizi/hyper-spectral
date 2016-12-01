function plot_arrow(x, y, string, linewidth)

    len_x = length(x);
   
    for i=1:len_x-1
        plot ([x(i), x(i+1)], [y(i), y(i+1)], string, 'linewidth', linewidth);
        pause(1)
    end
end