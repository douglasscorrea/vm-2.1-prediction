import sys
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

prediction_list = ['diffRDC']#['diffR', 'diffC', 'mule']
lambda_list = ['322_Bmin8']#['322', '8750', '127000', '1250000']

# prediction = sys.argv[1]
# _lambda = sys.argv[2]
for prediction in prediction_list:
	for _lambda in lambda_list:
		print(prediction, _lambda)
		partition_code_list = []
		with open('../partition_code_file_' + prediction + '_' + _lambda + '.txt', 'r') as fp:
			file_content = fp.readlines()

			for partition_code in file_content:
				partition_code = partition_code.rstrip('\n')
				#print(partition_code)

				partition_code_list.append(partition_code)

		#counts = dict(Counter(partition_code_list).most_common(10))
		counts = Counter(partition_code_list)
		partition_code, frequency = zip(*counts.items())
		index_sort = np.argsort(frequency)[::1]
		partition_code = np.array(partition_code)[index_sort]
		frequency = np.array(frequency)[index_sort]
		indexes = np.arange(len(partition_code))

		print(partition_code)

		if prediction == 'diffR':
			prediction_type = 'Differential Raster'
		elif prediction == 'diffC':
			prediction_type = 'Differential Central'
		elif prediction == 'mule':
			prediction_type = 'MuLE'
		elif prediction == 'diffRDC':
			prediction_type = 'Differential Raster DC Ref Plane'

		#plt.figure(figsize=(30, 10))
		plt.bar(indexes, frequency, tick_label=indexes)
		plt.xticks(indexes, partition_code)
		plt.title(prediction_type + ' - λ = ' +  _lambda + ' - Occurrences of each partition code', fontsize=14)
		plt.xlabel('Partition codes', fontsize=12)
		plt.ylabel('Occurrences', fontsize=12)
		plt.yscale('log')
		plt.xticks(rotation=90, fontsize=6)
		#plt.tight_layout()
		#ax = plt.gca()
		# ax.set_xticks(ax.get_xticks()[::2])
		plt.show()
		plt.savefig(prediction + '_' + _lambda + '.png', bbox_inches='tight', dpi=400)
		plt.close()