
# -*- coding: utf-8 -*-

prediction_list = ['diffCDC']
split_flags = ['1']
#superior_bitplanes = ['0', '2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '22', '24', '26', '28']
superior_bitplanes = ['30']
inferior_bitplanes = ['10', '12', '14', '16', '18', '20']

for prediction in prediction_list:
	for split in split_flags:
		for superior_bitplane in superior_bitplanes:
			for inferior_bitplane in inferior_bitplanes:
				with open('/home/douglascorrea/light-field/resultados_prediction/Bikes/decoding/bitplanes/summary/' + prediction + '_1_' + split + '_' + superior_bitplane + '_' + inferior_bitplane + '.txt') as fp:
					lines = fp.readlines()
					print(lines[-1].split(' ')[1].rstrip('\n') + ', ', end = '')