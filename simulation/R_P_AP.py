import os
from math import ceil
from numpy import mean

Dir = os.getcwd()
NF = Dir + '/navigation_footprint/'
cut = 10

AP_r = {} # AP_r = {flag:[run1, run2, ...], ...}
P_r = {} # P_r = {flag:[run1, run2, ...], ...}
R_r = {} # R_r = {flag:[run1, run2, ...], ...}
		
fList = os.listdir(NF)
fList.sort()
for filename in fList:
	(run,flag) = filename.split('.')
	run = int(run[3:])
	flag = int(flag[4:])
	
	filename = NF + filename
	fp = open(filename,'r')
	line = fp.readline()
	hit_order = eval(line)
	total_hit = len(hit_order)
	increment = ceil(total_hit/cut) + 1	
	
	i = 0
	c = 0
	SUM = 0
	total_sofar = 0
	r_sofar = 0
	total_r = total_hit - hit_order.count(0)
	
	print '~~~~~~~~~~~~~~ %s ~~~~~~~~~~~~~~~' % filename
	print '%d hits in total, %d relevant' % (total_hit,total_r)
	print 'cutoff\trecall\t\tprecision\tavg_p\t\ttotal_sofar\trelevant_sofar'
	

	
	while i < total_hit:
		r = hit_order[i]
		i += 1
		c += 1
		if r == 0:
			total_sofar += 1
		else:
			r_sofar += 1
			total_sofar += 1
			SUM += 1.0*r_sofar/total_sofar
		if c == increment:
			c = 0
			cutoff = int(i/increment)
			AP = 1.0*SUM/total_r
			R = 1.0*r_sofar/total_r
			P = 1.0*r_sofar/total_sofar
			
			print '%d\t\t%f\t%f\t%f\t%d\t\t%d' % (cutoff,R,P,AP,total_sofar,r_sofar)
			
	cutoff = int(i/increment) + 1
	AP = 1.0*SUM/total_r
	R = 1.0*r_sofar/total_r
	P = 1.0*r_sofar/total_sofar
	
	print '%d\t\t%f\t%f\t%f\t%d\t\t%d' % (cutoff,R,P,AP,total_sofar,r_sofar)
			
	if AP_r.has_key(flag):
		AP_r[flag].append(AP)
		P_r[flag].append(P)
		R_r[flag].append(R)
	else:
		AP_r[flag] = []
		P_r[flag] = []
		R_r[flag] = []

print '\n'
print '---------------------------------- AVG OF 20 RUN ---------------------------------'
print 'flag\trecall\t\tprecision\tavg_p'
for flag in AP_r:
	avg_AP = mean(AP_r[flag])
	avg_P = mean(P_r[flag])
	avg_R = mean(R_r[flag])
	print '%d\t\t%f\t%f\t%f' % (flag, avg_R, avg_P, avg_AP)	
