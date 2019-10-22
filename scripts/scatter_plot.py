import sys
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import os

prediction_list = ['diffRDC', 'diffR', 'mule']

inferior_bitplanes = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
mule_partition_numbers = [263, 265, 267, 274, 297, 256, 128, 684, 920, 1585, 1706, 2840, 3616, 3618, 3619]
diffRDC_partition_numbers = [420, 419, 425, 428, 454, 476, 629, 2229, 1755, 2124, 3374, 3584, 3594, 3600, 3610]
diffR_partition_numbers = [423, 423, 428, 435, 458, 484, 671, 2350, 2315, 2037, 2723, 3615, 3625, 3631, 3636]

#plt.figure(figsize=(30, 10))
plt.plot(inferior_bitplanes, mule_partition_numbers, '-or', label='MuLe')
for i, txt in enumerate(inferior_bitplanes):
	plt.annotate(txt, (inferior_bitplanes[i], mule_partition_numbers[i]))

plt.plot(inferior_bitplanes, diffRDC_partition_numbers, '-ob', label='diffRDC')
for i, txt in enumerate(inferior_bitplanes):
	plt.annotate(txt, (inferior_bitplanes[i], diffRDC_partition_numbers[i]))


plt.plot(inferior_bitplanes, diffR_partition_numbers, '-og', label='diffR')
for i, txt in enumerate(inferior_bitplanes):
	plt.annotate(txt, (inferior_bitplanes[i], diffR_partition_numbers[i]))

plt.xlabel('Bmin', fontsize=12)
plt.legend(fontsize=16)
plt.tight_layout()
plt.ylabel('Number of partitionings', fontsize=12)
plt.show()
#