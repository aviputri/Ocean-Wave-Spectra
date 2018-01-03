#!/usr/bin/python
import sys, os

if __name__ == '__main__':
	# change the content of '' into the input file folder ! THIS IS JUST AN EXAMPLE !
	src_path = 'C:\\Users\\Administrator\\Desktop\\'
	# change the content of '' into the output file folder ! THIS IS JUST AN EXAMPLE !
	dest_path = 'C:\\Users\\Administrator\\Desktop\\data\\'

	# change the content of '' into the name of the file ! THIS IS JUST AN EXAMPLE !
	with open(os.path.join(source_path, '1-NBU-Sensor1.txt'), 'r') as fread:
		raw = [line.lstrip() for line in fread]

	header = raw[0]
	raw = raw[1:]

	# positions of each value
	# :8 --> date (0-8)
	# 9:19 --> time (9-19)
	# 20: --> elevation (20-end)
	lines = [[line[:8], line[9:19], line[20:].lstrip()] for line in raw]
	del raw

	nline = 0
	counter = 0
	filecounter = 0

	while (nline <len(lines)):
		file_path = ''
		# naming files from 0000 - 0009
		if filecounter <10:
			file_path = '000' + str(filecounter)
		# naming files from 0010 - 0099
		if filecounter <100:
			file_path = '00' + str(filecounter)
		# naming files from 0100 - 0999
		if filecounter <1000:
			file_path = '0' + str(filecounter)
		# naming files from 1000 - 2330 <--in my case
		else:
			file_path = str(filecounter)
		with open(os.path.join(folder_path, file_path), 'w') as fwrite:
			fwrite.write('time,elevation\n')
			for n in xrange(nline, nline+7200):
				fwrite.write(str(counter) + ',' + lines [n][2])
				counter += 0.5 # interval of 0.5 seconds
		nline += 7200 # the number of lines in each file is 7200
		filecounter += 1