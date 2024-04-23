import os
config_file = "./blislab/bl_config.h"
f = open(config_file,'r')
f_str = f.read()
g = f_str.split('\n')
f.close()

def replace_para(list_str,para_list, val_list):
	for index_p, paras in enumerate(para_list):
		for index, items in enumerate(list_str):
			if(items.find(paras) != -1 and not(items.strip().startswith('//'))):
				list_str[index] = "#define "+str(paras)+" "+str(val_list[index_p])

def print_para(list_str,para_list):
	for index_p, paras in enumerate(para_list):
		for index, items in enumerate(list_str):
			if(items.find(paras) != -1 and not(items.strip().startswith('//'))):
				print(items)

para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC','DGEMM_MR','DGEMM_NR']
N = 511
name = str(N) + 'KC'
os.system('echo "" > result_autotune.'+str(name))
for i in range(1, 128, 4):
	f = open(config_file,'r')
	f_str = f.read()
	g = f_str.split('\n')
	f.close()

	val_list = [1920,i,256,8,4]
	replace_para(g,para_list,val_list)

	f = open(config_file,'w')
	ans = ""
	for items in g:
		ans = ans + items + '\n'
	f.write(ans)
	f.close()

	os.system('make')
	os.system('echo "PARA1 '+str(i)+'" >> result_autotune.'+str(name))
	cmd = './benchmark-blislab -n '+str(N)+' >> result_autotune.'+str(name)
	os.system(cmd)
