function[] = tofplot_batch(name, root)

if (strcmp(name, '') == 1 || strcmp(root, '') == 1)
    error('Usage: tofplot_batch name root');
end

[histpixel, nsteps] = getparams(name, root);

for idx = 0:nsteps-1
    filename = strcat(root, '/', name, '/histograms/output/', name, '_', num2str(idx), '.txt');
    fileID = fopen(filename, 'r');
    if (fileID == -1)
        error('Could not open %s', filename);
    end
   
    pixels = textscan(fileID, '%s', 'delimiter', '#');
    fclose(fileID);
    
    for j = 2:numel(pixels{1})
       [pixel, pos] = textscan(pixels{1}{j}, '%d%d', 1);
       data = textscan(pixels{1}{j}(pos+1:end),'%f%f');
       
       for k = 1:length(histpixel{1})
           if (pixel{1,1} == histpixel{1,1}(k) && pixel{1,2} == histpixel{1,2}(k))
               D = data{1,1};
               L = data{1,2};

               fig_dir = strcat(root, '/', name, '/histograms/images/figures/');
               if ~(exist(fig_dir, 'dir') == 7)
                   error('%s is not a directory', outdir);
               end

               jpeg_dir = strcat(root, '/', name, '/histograms/images/jpeg/');
               if ~(exist(jpeg_dir, 'dir') == 7)
                   error('%s is not a directory', outdir);
               end

               h = figure;
               plot(D, L);

               saveas(h, strcat(fig_dir, name, '_', num2str(pixel{1,1}), ...
                   '_', num2str(pixel{1,2}), '_', num2str(idx)), 'fig');
               saveas(h, strcat(jpeg_dir, name, '_', num2str(pixel{1,1}), ...
                   '_', num2str(pixel{1,2}), '_', num2str(idx)), 'jpeg');
           end
       end
    end
end

close all;

end