function[] = tofplot(input, pixelX, pixelY, output)
    fileID = fopen(input, 'r');
    if (fileID == -1)
        error('Could not open %s', input);
    end

    pixels = textscan(fileID, '%s', 'delimiter', '#');
    fclose(fileID);

    for j = 2:numel(pixels{1})
       [pixel, pos] = textscan(pixels{1}{j}, '%d%d', 1);
       data = textscan(pixels{1}{j}(pos+1:end),'%f%f');

       if (pixel{1,1} == pixelX && pixel{1,2} == pixelY)
           D = [data{1,1}];
           L = [data{1,2}];

           h = figure;
           plot(D, L);

           saveas(h, output, 'fig');
           saveas(h, output, 'jpeg');
       end
    end

    %close all;

end