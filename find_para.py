L3 = pow(2,12) * pow(2,10)
L2 = pow(2,7) * pow(2,10)
L1 = pow(2,3) * pow(2,10)
REGS =  32

N = 4096
#mc, kc, nc, mr, nr

def obj_func(mc,kc,nc,mr,nr):
	return mc*kc + nc*kc + mr*kc + mr*nr
DIV = 64
max_val = 0
for mc in range(int(N/DIV),N):
	for nc in range(int(N/DIV),N):
		for kc in range(int(N/DIV),N):
			for mr in range(int(REGS/2),REGS):
				for nr in range(int(REGS/2),REGS):
					#if(mc*kc <= L3 and nc*kc <= L2 and mr*kc <= L1 and mr*nr <= 32 and obj_func(mc,kc,nc,mr,nr) >= max_val/2 ):
					if(mc*kc <= L3 and nc*kc <= L2 and mr*kc <= L1 and mr*nr <= 32 ):
						max_val = obj_func(mc,kc,nc,mr,nr)
						print(mc,kc,nc,mr,nr,max_val,L3+L2+L1+REGS)

