lambdas = [322, 5800, 127000, 1250000]
exec_names = ['./mule-encoder-bin_MuLE', './mule-encoder-bin_DC', './mule-encoder-bin_DR']

for name in exec_names:
	for _lambda in lambdas:
		alg_name = name.split('_')[1]
		cmd1 = name + ' -lf /home/douglascorrea/light-field/dataset/Lenslet/Danger_de_Mort/ '
		cmd2 = '-o /home/douglascorrea/light-field/resultados_prediction/Danger_de_Mort/coding/output_' + alg_name + '_' + str(_lambda) + '.LF '
		cmd3 = '-cf /home/douglascorrea/GitHub/vm-2.1-prediction/MULE_CFGs/Bikes/I01_Bikes_22016.json '
		cmd4 = '-lambda ' + str(_lambda)
		print(cmd1 + cmd2 + cmd3 + cmd4)

