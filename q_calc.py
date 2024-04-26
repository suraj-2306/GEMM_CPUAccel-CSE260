import os
import copy
config_file = "./blislab/bl_config.h"
#f = open(config_file,'r')
#f_str = f.read()
#g = f_str.split('\n')
#f.close()

ori_val = []

def replace_para(para_list, val_list):
	f = open(config_file,'r')
	f_str = f.read()
	list_str = f_str.split('\n')
	f.close()
	for index_p, paras in enumerate(para_list):
		for index, items in enumerate(list_str):
			if(items.find(paras) != -1 and not(items.strip().startswith('//'))):
				list_str[index] = "#define "+str(paras)+" "+str(val_list[index_p])
	f = open(config_file,'w')
	ans = ""
	for items in list_str:
		ans = ans + items + '\n'
	f.write(ans)
	f.close()

def print_para(para_list):
	f = open(config_file,'r')
	f_str = f.read()
	list_str = f_str.split('\n')
	f.close()
	for index_p, paras in enumerate(para_list):
		for index, items in enumerate(list_str):
			if(items.find(paras) != -1 and not(items.strip().startswith('//'))):
				print(items)

def save_ori_config():
	para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC', 'DGEMM_MR', 'DGEMM_NR']
	f = open(config_file,'r')
	f_str = f.read()
	list_str = f_str.split('\n')
	f.close()
	global ori_val
	for index_p, paras in enumerate(para_list):
		for index, items in enumerate(list_str):
			if(items.find(paras) != -1 and not(items.strip().startswith('//'))):
				#print(items.strip().split()[-1])
				ori_val.append(items.strip().split()[-1])
save_ori_config()

def reset_config():
	para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC', 'DGEMM_MR', 'DGEMM_NR']
	val_list = [64, 64, 64, 4, 4]
	replace_para(para_list,val_list)

def restore_config():
	para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC', 'DGEMM_MR', 'DGEMM_NR']
	val_list = ori_val
	replace_para(para_list,val_list)

def write_para(para_list,params):
	val_list = copy.deepcopy(params)
	replace_para(para_list,val_list)


#para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC','DGEMM_MR','DGEMM_NR']
para_list = ['DGEMM_MC', 'DGEMM_KC', 'DGEMM_NC', 'DGEMM_MR', 'DGEMM_NR']
N = 1024
TAG = os.getcwd().split('/')[-1].replace('pa1-kdivij-ssathyaprakash.','')
name = str(N) + "_q_calc."+TAG+".txt" 
mat_opt = [1024, 332, 314, 16, 4]
config_list = [
    [10, 10, 10, 16, 4],
    [20, 17, 17, 16, 4],
    [38, 29, 27, 16, 4],
    [73, 47, 45, 16, 4],
    [141, 79, 73, 16, 4],
    [273, 132, 119, 16, 4],
    [529, 220, 196, 16, 4],
    [1025, 368, 320, 16, 4]
]
os.system('echo "" > result_autotune.'+str(name))
for index, params in enumerate(config_list):
	print(params)
	write_para(para_list,params)
	os.system('make > mf 2>&1')
	os.system('echo "'+str(params)+'" >> result_autotune.'+str(name))
	cmd = './benchmark-blislab -n '+str(N)+'| grep Geo >> result_autotune.'+str(name)
	os.system(cmd)

restore_config()
