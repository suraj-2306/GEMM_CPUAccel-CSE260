import os
import copy
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

def reset_config():
	f = open(config_file,'r')
	f_str = f.read()
	g = f_str.split('\n')
	f.close()

	para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC']
	val_list = [1024, 332, 314, 16, 4]
	replace_para(g,para_list,val_list)
	
	f = open(config_file,'w')
	ans = ""
	for items in g:
		ans = ans + items + '\n'
	f.write(ans)
	f.close()

#para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC','DGEMM_MR','DGEMM_NR']
para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC']
mat_opt = [1024, 332, 314, 16, 4]
for index, params in enumerate(para_list):
	N = 513
	name = str(N) + "_" + params
	os.system('echo "" > result_autotune.'+str(name))
	reset_config()
	for i in range(-256, 256, 16):
		f = open(config_file,'r')
		f_str = f.read()
		g = f_str.split('\n')
		f.close()
	
		val_list = copy.deepcopy(mat_opt)
		val_list[index] = val_list[index] + i;
		replace_para(g,para_list,val_list)
	
		f = open(config_file,'w')
		ans = ""
		for items in g:
			ans = ans + items + '\n'
		f.write(ans)
		f.close()
	
		os.system('make > mf 2>&1')
		print(params+" "+str((i+256+16)/5.12)+"%")
		os.system('echo "'+params+' '+str(val_list[index])+'" >> result_autotune.'+str(name))
		cmd = './benchmark-blislab -n '+str(N)+'| grep Geo >> result_autotune.'+str(name)
		os.system(cmd)
