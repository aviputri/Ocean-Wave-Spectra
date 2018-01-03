function fx = lagrpol(x, y, n, xx)
	fx = 0.0;
	l = 1.0;

	for i = 1:n
		l = 1.0;
		for j = 1:n
			if (j ~= i)
				l = l * (xx - x(j)) / (x(i) - x(j));
			end
		end
		fx = fx + (l * y(i));
	end