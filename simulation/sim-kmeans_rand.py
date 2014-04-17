# Simulate user navigation based on SOM clustering result. 
# tab width:4
# external files necessary: cluster distribtion, distance, line2ID map, ECT.hits
# docID range(1~29381) in database

import operator
import _mysql
import MySQLdb as mdb

"""
functionc in this script:

read_as_mt		(input_csv,sep,head):										return data
som				(data,xdim,ydim,topo,seed=None):							return sm
sort_sm			(sm,NO_cluster,clusterID = None,NO_FILE=False):				return sorted_result
list_hits		(sorted_result,clusterID=None,NO_FILE=False):				return hit_o_in_cluster
revise_order	(hit_order_navg):											return revised_order
simulate		(hit_o_in_cluster,flag,p,NO_cluster,NO_FILE=False):			return hit_order_navg
re_som			(clusterID,clusters,data,xdim,ydim):						return new_cluster,hit_o_in_cluster
count_R_P_AP	(hit_order,flag):											return (R,P,AP)
log				(hit_o_in_cluster,sorted_result,NO_FILE=False)				return (D,monster,cluster_d,total_hit,cluster_log)
"""

#NO_cluster = 200		# initial cluster number
ID = 'kmeans.rand4.docF2'	# identifier of each result



def print_time():
	"""
	TO: record run time
	"""
	from datetime import datetime
	return datetime.now()

	
def readDocID(filename):
	line2docid = {}			# line2docid = {line number:docid, ...}
	fp = open(filename,'r')
	lineID = 1
	for line in fp:
		docid = eval(line)
		line2docid[lineID] = docid
		lineID += 1
	fp.close()
	return line2docid


def readClusterDistribution(filename,line2docid):
	cluster = {}			# cluster = {clusterID:[docid,...],...}		here clusterID: 0~199
	fp = open(filename,'r')
	lineID = 1
	for line in fp:
		clusterID = eval(line)
		docid = line2docid[lineID]
		lineID += 1
		if cluster.has_key(clusterID):
			cluster[clusterID].append(docid)
		else:
			cluster[clusterID] = [docid]
	total = 0
	for i in cluster:
		print 'cluster %d has %d docs' % (i,len(cluster[i]))
		total += len(cluster[i])
	print '%d docs in total' % total
	fp.close()
	return cluster


def sort_sm(filename,cluster,line2docid):
	sorted_cluster = {}		# sorted_cluster = {clusterID:[(docid,dist),...],...} sort by distance
	fp = open(filename,'r')
	docid2dist = {}			# docid2dist = {docid:distance,...}
	lineID = 1
	for line in fp:
		try:
			dist = eval(line)
		except:
			dist = 100000000000000
		docid = line2docid[lineID]
		docid2dist[docid] = dist
		lineID += 1
	fp.close()
	for clusterID in cluster.keys():
		docs = cluster[clusterID]
		tmp = {}		# tmp = {docid:distance,...} for docs in a cluster
		for docid in docs:
			dist = docid2dist[docid]
			tmp[docid] = dist
		tmp = sorted(tmp.iteritems(), key=operator.itemgetter(1))
		sorted_cluster[clusterID] = tmp
	return sorted_cluster


def list_hits(sorted_cluster):
	try:
		conn = mdb.connect('localhost', 'root', 'fdsarewq', 'PATENTS_ECT');	# connect to database
		cursor = conn.cursor()
		output = 'hit&relevance_table-'+ID+'.txt'		# list docID, hitID, hit_offset, hits_relevant group by cluster
		fp2 = open(output,'w')
		hit_o_in_cluster = {}		# hit_o_in_cluster = {clusterID:[relevance,...]} order by docid & offset

		for clusterID in sorted_cluster.keys():	# each cluster
			hit_o_in_cluster[clusterID] = []
			ls = sorted_cluster[clusterID]	# ls = [(docid,distance),...]
			fp2.write('-------------------'+str(clusterID)+'-------------------------\n')
			fp2.write('docID\t\thitID\t\thit_offset\t\thits_relevant\n')
			hit_o_in_cluster[clusterID] = []
			for (docid,dist) in ls:
				cursor.execute('SELECT docid, CAST(FileOffset as UNSIGNED), Relevance FROM hits \
				WHERE documenttabledocid = %d ORDER BY CAST(FileOffset as UNSIGNED)' % docid) # fetch hits in the doc
				rows = cursor.fetchall()
				for row in rows: 	# row = each hit 
					(hitid, offset, relevance) = row
					fp2.write('%d\t%d\t%d\t%d\n' % (docid,hitid,offset,relevance))
					hit_o_in_cluster[clusterID].append(relevance)						

		conn.commit()
		cursor.close()
		conn.close()
		fp2.close()
		
	except mdb.Error, e: 
		print "Error %d: %s" % (e.args[0],e.args[1])
		sys.exit(1)

	return hit_o_in_cluster


def count_R_P_AP(hit_order):
	"""
	TO: count R,P and AP in a hit relevance order
	INPUT:	hit_order	list	[0,0,0,1,0,0,1,0,1,0,0,0,1,...] 
	OUTPUT:	R,P,AP		tuple	
	"""	
	SUM = 0
	total_sofar = 0
	r_sofar = 0
	total_hit = len(hit_order)
	if total_hit == 0:
		return (0,0,0)
	else:
		total_r = total_hit-hit_order.count(0)
		for i in hit_order:
			if i == 0:
				total_sofar += 1
			else:
				r_sofar += 1
				total_sofar += 1
				SUM += 1.0*r_sofar/total_sofar
		if total_r == 0:
			(R,P,AP) = (0,0,0)
		else:
			AP = 1.0*SUM/total_r
			R = 1.0*r_sofar/total_r
			P = 1.0*total_r/total_hit
		return (R,P,AP)


def log(hit_o_in_cluster,sorted_cluster,NO_FILE=False):
	"""
	To: output clustet log
	INPUT:	hit_o_in_cluster	dict		{clusterID:[hit_r (in order)]} after some been re_clustered and sorted
	OUTPUT:	cluster_log			dict		{clusterID:[NO_DOC,NO_HIT,DEVIATION,PERCENTAGE,PRECISION,AP,DISTANCE(MIN),DISTANCE(MAX),DISTANCE(MEAN),DISTANCE(VAR)]}
			D					float		standard deviation on NO_HIT			
			cluster_d			dict		{clusterID:NO_HIT,...,'STD':D}
			monster				list		[clusterID,...] clusters whose NO_HIT >= 3D
			total_hit			int			total number of hits in the cluster set
	FILE:	cluster_log			txt			table of NO_DOC,NO_HIT,DEVIATION,PERCENTAGE,PRECISION,AP,DISTANCE(MIN),DISTANCE(MAX),DISTANCE(MEAN),DISTANCE(VAR)			
	"""

	from numpy import mean,var,std

	cluster_log = {}	# cluster_log = {clusterID:[NO_DOC,NO_HIT,DEVIATION,PERCENTAGE,PRECISION,AP,DISTANCE(MIN),DISTANCE(MAX),DISTANCE(MEAN),DISTANCE(VAR)]}
	cluster_d = {}		# cluster_d = {clusterID:NO_HIT,...,STD:DEVIATIOn}
	monster = []		# monster = [clusterID,clusterID]

	total_hit = 0
	devi = []
	for clusterID in hit_o_in_cluster:
		hit_order = hit_o_in_cluster.get(clusterID)
		(R,P,AP) = count_R_P_AP(hit_order)
		s_P = format(P,'%')
		s_AP = format(AP,'%')
		NO_HIT = len(hit_o_in_cluster.get(clusterID))
		total_hit += NO_HIT
		cluster_log[clusterID] = [NO_HIT,s_P,s_AP]
		devi.append(NO_HIT)		# get NO_HIT in each cluster to count deviation
		cluster_d[clusterID] = NO_HIT

	D = std(devi)
	cluster_d['STD'] = D
	del devi
	
	for clusterID in sorted_cluster:
		docs = sorted_cluster.get(clusterID)
		NO_DOC = len(docs)
		dt = []
		for (docID,dist) in docs:	# get distance list in the cluster
			dt.append(dist)
		MAX = max(dt)
		MIN = min(dt)
		MEAN = mean(dt)
		VAR = var(dt)
		
		NO_HIT = cluster_log.get(clusterID)[0]
		DEVIATION = NO_HIT/D
		if DEVIATION > 3:
			monster.append(clusterID)
		p = 1.0*NO_HIT/total_hit
		s_p = format(p,'%')

		cluster_log[clusterID].insert(0,NO_DOC)		# [NO_DOC,NO_HIT,P,AP]
		cluster_log[clusterID].insert(2,DEVIATION)		# [NO_DOC,NO_HIT,DEVIATION,P,AP]		
		cluster_log[clusterID].insert(3,s_p)		# [NO_DOC,NO_HIT,DEVIATION,PERCENTAGE,P,AP]
		cluster_log[clusterID].extend([MIN,MAX,MEAN,VAR])		# [NO_DOC,NO_HIT,DEVIATION,PERCENTAGE,P,AP,MIN,MAX,MEAN,VAR]	
	
	if NO_FILE == False:
		filename = 'cluster_log_'+ID+'.txt'
		fp = open(filename,'w')
		fp.write('ID\tNO_DOC\tNO_HIT\tSTD(times)\tPERCENTAGE\tPrecesion\tAvgPrecesion\tDIST(MIN)\tDIST(MAX)\tDIST(MEAN)\tDIST(VAR)\n')
		for k in cluster_log:
			[NO_DOC,NO_HIT,DEVI,PERCENTAGE,P,AP,MIN,MAX,MEAN,VAR] = cluster_log.get(k)	
			fp.write('%s\t%d\t%d\t%f\t%s\t%s\t%s\t%f\t%f\t%f\t%f\n' % (k,NO_DOC,NO_HIT,DEVI,PERCENTAGE,P,AP,MIN,MAX,MEAN,VAR))
		fp.close()
	else:
		pass

	return (D,monster,cluster_d,total_hit,cluster_log)


def simulate(hit_o_in_cluster,flag,p,NO_cluster,NO_FILE=False):
	"""
	INPUT:	hit_o_in_cluster	dict	{clusterID:[hit_r (in order)]}
			flag				int		threashold of least number of hits when leave a cluster
			p					float	threashold of precision when leave a cluster
	OUTPUT:	hit_order_navg	list	[hit_r (in navigation order)]
	FILES:	out					txt		simulation log: time, hit_r in each cluster in navigation order, total number of hits, number of relevant & non-relevant hits
	"""

	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for i in complete:
			fp.write('cluster:%s\n' % str(i))
			ls = hit_o_in_cluster.get(i)
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			for r in ls:
				total_sofar += 1
				counter += 1
				tmp.append(r)
				if r == 0:
					hit_order_navg.append(0)
				else:	# r == 1
					NO_r += 1
					total_r_sofar += 1
					hit_order_navg.append(total_r_sofar)
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
	else:
		for i in complete:
			ls = hit_o_in_cluster.get(i)
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			for r in ls:
				total_sofar += 1
				counter += 1
				tmp.append(r)
				if r == 0:
					hit_order_navg.append(0)
				else:	# r == 1
					NO_r += 1
					total_r_sofar += 1
					hit_order_navg.append(total_r_sofar)
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					break
				else:
					continue		

	return hit_order_navg

def read_hit_order(filename):
	total = 0
	ct_docs = 0
	hit_o_in_cluster = {}		# hit_o_in_cluster = {clusterID:{docid:[relevance,...]},...} order by docid & offset
	fp = open(filename,'r')
	clusterID = 0
	for line in fp:
		if line[0] == '-':
			clusterID += 1
#			print clusterID
			hit_o_in_cluster[clusterID] = {}
		elif line[0] == 'd':
			continue
		else:
			(docID,hitID,hit_offset,hits_relevant) = line.split('\t')
			hits_relevant = eval(hits_relevant)
			if hit_o_in_cluster[clusterID].has_key(docID):				
				hit_o_in_cluster[clusterID][docID].append(hits_relevant)
			else:
				hit_o_in_cluster[clusterID][docID] = [hits_relevant]
				ct_docs += 1
				
			total += 1
	print 'total hits: %d' % total
	print 'total docs: %d' % ct_docs
	print 'total clusters: %d' % len(hit_o_in_cluster)
	return hit_o_in_cluster
	
def simulate_random(hit_o_in_cluster,flag,p,NO_cluster,NO_FILE=False):	# randomly select a clusterID, then randomly select a document, navigate from the first hit in that doc, if relevant, keep on going, otherwise go to another doc
	
	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key/clusterID
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for clusterID in complete:
			fp.write('cluster:%s\n' % str(clusterID))
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
#			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
			
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						break	# go to another doc
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
		
	else:		# no file
		for clusterID in complete:
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
			
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						break	# go to another doc
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue

	return hit_order_navg


def simulate_random2(hit_o_in_cluster,flag,p,NO_cluster,doc_flag, NO_FILE=False):	# randomly select a clusterID, then randomly select a document, navigate from the first hit in that doc, if relevant, keep on going, otherwise, keep on going untill doc_flag of non-relevant hits has reached.
	
	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key/clusterID
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for clusterID in complete:
			fp.write('cluster:%s\n' % str(clusterID))
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
				non_re = doc_flag
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						if non_re == 0:
							break	# go to another doc
						else:
							non_re -= 1	# keep on going
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag*doc_flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
	return hit_order_navg	

def simulate_random3(hit_o_in_cluster,flag,p,NO_cluster,NO_FILE=False):	# randomly select a clusterID, then randomly select a document, navigate from the a random hit in that doc, if relevant, keep on going to next random hit, otherwise go to another doc
	
	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key/clusterID
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for clusterID in complete:
			fp.write('cluster:%s\n' % str(clusterID))
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
				shuffle(re_ls)
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						break	# go to another doc
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
	return hit_order_navg


def simulate_random4(hit_o_in_cluster,flag,p,NO_cluster,doc_flag,NO_FILE=False):	# randomly select a clusterID, then randomly select a document, navigate from the a random hit in that doc, if relevant, keep on going to next random hit, otherwise exit doc if doc_flag non-relevant hits found
	
	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key/clusterID
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for clusterID in complete:
			fp.write('cluster:%s\n' % str(clusterID))
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
				non_re = doc_flag
				shuffle(re_ls)
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						if non_re == 0:
							break	# go to another doc
						else:
							non_re -= 1	# keep on going
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
	return hit_order_navg
	

def simulate_random5(hit_o_in_cluster,flag,p,NO_cluster,doc_flag,NO_FILE=False):	# randomly select a clusterID, then randomly select a document, navigate from the a random hit in that doc, if relevant, keep on going sequenctially, otherwise go to next randmo hit. exit doc if doc_flag non-relevant hits found
	
	from random import shuffle
	complete = hit_o_in_cluster.keys()	# get a full list of key/clusterID
	shuffle(complete)	# sort the list randomly

	total_r_sofar = 0
	hit_order_navg = []

	if NO_FILE == False:
		filename = 'sim_out_'+ID+'.txt'
		fp = open(filename, 'a')

		for clusterID in complete:
			fp.write('cluster:%s\n' % str(clusterID))
			tmp = []
			total_sofar = 0
			NO_r = 0
			counter = 0
			P = 0
			
			ls = hit_o_in_cluster.get(clusterID)
			docs = ls.keys()
			shuffle(docs)
			for docID in docs:
				re_ls = hit_o_in_cluster.get(clusterID).get(docID)
				non_re = doc_flag
				shuffle(re_ls)
				for r in re_ls:
					total_sofar += 1
					counter += 1
					tmp.append(r)
					if r == 0:		# if non-relevant
						hit_order_navg.append(0)
						if non_re == 0:
							break	# go to another doc
						else:
							non_re -= 1	# keep on going
					else:			# if relevant
						NO_r += 1
						total_r_sofar += 1
						hit_order_navg.append(total_r_sofar)
						
				P = 1.0*NO_r/total_sofar
				if counter >= flag and P < p:
					fp.write('%d hits navigated through\n' % counter)
					fp.write('hit nagivated in this cluster: %s' % tmp)
					fp.write('\n')
					break
				else:
					continue
		fp.write('**************************************************************************************\n')
	return hit_order_navg
						
def main():
	ap = {}	# ap = {flag:[avg ap,...]}
	i = 0
	doc_flag = 2
#	line2docid = readDocID('line2docid.txt')
#	cluster = readClusterDistribution('cluster',line2docid)
#	sorted_cluster = sort_sm('distance',cluster,line2docid)
#	hit_o_in_cluster = list_hits(sorted_cluster)
#	(D,monster,cluster_d,total_hit,cluster_log) = log(hit_o_in_cluster,sorted_cluster,NO_FILE=False)
	
	hit_o_in_cluster = read_hit_order('hit&relevance_table-10.18.kmeans.txt')
	Run_ct = 20
	while i < Run_ct:
		i += 1

		NO_cluster = INI_NO_cluster
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~run %d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' % i
		print '*****************************************************************'
		print '*START TIME!\t\t\t%s\t*' % print_time()
		
	
		p = 0.5
		Flag_ct = 10
		for flag in range(10,50*doc_flag):
			f = 0
			s = 0
			print '***********flag == %d*************' % flag
			while f < Flag_ct:
				hit_order = simulate_random4(hit_o_in_cluster,flag,p,NO_cluster,doc_flag,NO_FILE=False)
				total = len(hit_order)
				Non_r = hit_order.count(0)
				NO_r = total-Non_r	
				(R,P,AP) = count_R_P_AP(hit_order)
				s += AP
				print '%d hits intotal,%d non-relevant,%d relevant' % (total,Non_r,NO_r)
				f += 1
			avg = s/Flag_ct
			print 'flag: %d\tAvg Ap: %f' % (flag,avg)

			if ap.has_key(flag):
				ap[flag].append(avg)
			else:
				ap[flag] = [avg]
	
	print '!!!!!!!!!!!!!!!FINAL AVG AP!!!!!!!!!!!!!!'
	print ap
	for flag in ap:
		avg_avg_ap = sum(ap.get(flag))/Run_ct
		print 'flag: %d\tAvg Avg Ap: %f' % (flag,avg_avg_ap)
	
if __name__ == '__main__':
	INI_NO_cluster = 200
#	NO_FILE = True
	main()
