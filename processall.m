function [M_Gmean, M_Gmax, M_GHs, M_GH10] = processall(path)
	addpath(path);
	files = dir(fullfile(path));
	n = size(files);
	n = n(1) - 2;

	M_Gmean = [0,0];
	M_Gmax = [0,0];
	M_GHs = [0,0];
	M_GH10 = [0,0];

	for i = 1:n
		m = tdfread(files(i+2).name, ',');
		E = m.elevation;
		[M_Gmean(i,:), M_Gmax(i,:), M_GHs(i,:), M_GH10(i,:)] = zerodown(E,0.5);
	end

        [success,message] = xlswrite('hasil.xls', M_Gmean, 1);
        [success,message] = xlswrite('hasil.xls', M_Gmax, 2);
        [success,message] = xlswrite('hasil.xls', M_GHs, 3);
        [success,message] = xlswrite('hasil.xls', M_GH10, 4);
    
