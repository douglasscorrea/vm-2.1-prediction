import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
import glob
import statistics
import os

figures_path = '/home/douglascorrea/light-field/resultados_prediction/Bikes/coding/rate_per_block_figures/'
inferior_bitplanes = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
mule_average_list = []
diffR_average_list = []
diffC_average_list = []
mule_std_list = []
diffR_std_list = []
diffC_std_list = []

for inferior_bitplane in inferior_bitplanes:
	print('Bmin = ' + str(inferior_bitplane))
	with open('/home/douglascorrea/light-field/resultados_prediction/Bikes/coding/ratePerBlock/ratePerBlock_mule_1_1_30_' + str(inferior_bitplane) + '.txt') as fp:
		content = fp.readlines()

		strip_list = [item.strip() for item in content]
		float_list = [float(item) for item in strip_list]
		
		mu = statistics.mean(float_list)
		mule_average_list.append(mu)
		variance = statistics.variance(float_list)
		sigma = math.sqrt(variance)
		mule_std_list.append(sigma)
		print('\tMuLE - Average: ' + str(mu) + ' - Std: ' + str(sigma))

		x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
		plt.plot(x, stats.norm.pdf(x, mu, sigma), label='MuLE')

	with open('/home/douglascorrea/light-field/resultados_prediction/Bikes/coding/ratePerBlock/ratePerBlock_diffR_1_1_30_' + str(inferior_bitplane) + '.txt') as fp:
		content = fp.readlines()

		strip_list = [item.strip() for item in content]
		float_list = [float(item) for item in strip_list]
		
		mu = statistics.mean(float_list)
		diffR_average_list.append(mu)
		variance = statistics.variance(float_list)
		sigma = math.sqrt(variance)
		diffR_std_list.append(sigma)
		print('\tdiffR - Average: ' + str(mu) + ' - Std: ' + str(sigma))

		x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
		plt.plot(x, stats.norm.pdf(x, mu, sigma), label='diffR')

	with open('/home/douglascorrea/light-field/resultados_prediction/Bikes/coding/ratePerBlock/ratePerBlock_diffC_1_1_30_' + str(inferior_bitplane) + '.txt') as fp:
		content = fp.readlines()

		strip_list = [item.strip() for item in content]
		float_list = [float(item) for item in strip_list]
		
		mu = statistics.mean(float_list)
		diffC_average_list.append(mu)
		variance = statistics.variance(float_list)
		sigma = math.sqrt(variance)
		diffC_std_list.append(sigma)
		print('\tdiffC - Average: ' + str(mu) + ' - Std: ' + str(sigma))

		x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
		plt.plot(x, stats.norm.pdf(x, mu, sigma), label='diffC')

	plt.legend()
	plt.title('Rate per block - Bmin = ' + str(inferior_bitplane))
	plt.savefig('rate_per_block.png', bbox_inches='tight', dpi=400)

	cmd = 'mv rate_per_block.png ' + figures_path + 'ratePerBlock_' + str(inferior_bitplane) + '.png'
	os.system(cmd)

	#plt.show()
	plt.close()

print('mule_rate_per_block_average = ' + str(mule_average_list))
print('diffR_rate_per_block_average = ' + str(diffR_average_list))
print('diffC_rate_per_block_average = ' + str(diffC_average_list))
# print('mule_rate_per_block_std = ' + str(mule_std_list))
# print('diffR_rate_per_block_std = ' + str(diffR_std_list))
# print('diffC_rate_per_block_std = ' + str(diffC_std_list))
