function[histogram_pixels, nsteps] = ...
    getparams(name, root)

    filename = strcat(root, '/', name, '/', name, '_template.pbrt');
    fileID = fopen(filename, 'r');
    if (fileID == -1)
        error('Could not open %s', filename);
    end
    
    lines = textscan(fileID, '# %s', 'delimiter', '\n');
    
    nsteps = textscan(lines{1}{2}, '%f');
    hpixel_x = textscan(lines{1}{5}, '%f');
    hpixel_y = textscan(lines{1}{6}, '%f');
    
    histogram_pixels = [hpixel_x, hpixel_y];
    nsteps = nsteps{1};
    
    fclose(fileID);
end