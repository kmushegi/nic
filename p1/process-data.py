import sys
import os

statsdir = "stats-ga-extended-2"

def findMaxIndex(pcs):
	maxIndex = 0
	maxValue = pcs[0]

	for i in range(1,len(pcs)):
		if pcs[i] > maxValue:
			maxIndex = i
			maxValue = pcs[i]
		else:
			continue
	return maxIndex

def findMinIndex(pcs):
	minIndex = 0
	minValue = pcs[0]

	for i in range(1,len(pcs)):
		if pcs[i] < minValue:
			minIndex = i
			minValue = pcs[i]
		else:
			continue
	return minIndex

def findAverage(pcs):
	sum = 0.0
	for i in range(0,len(pcs)):
		sum += float(pcs[i])

	avg = sum/len(pcs)
	return avg

def analyzedata(sd):

	files = []
	percents = []
	times = []
	iterations = []

	for filename in os.listdir(sd): #for every file in the directory
		if filename.endswith(".txt"): #with extension .txt
			files.append(filename)

			temp_percents = []
			temp_times = []
			temp_iterations = []

			with open(sd+"/"+filename) as f: #read file line by line
				lineCount = 0
				for line in f:
					if line != "\n": #if its not just an empty line
						line = line.replace("\n","")
						#print line
						lineCount += 1 #count number of lines of error checking
						tokens = line.split(" ")

						temp_percents.append(tokens[2])
						temp_iterations.append(tokens[3])
						temp_times.append(tokens[4])
						#print tokens[2]

			if lineCount != 48: #there are 48 problems per experiment
				print "Number of Lines Looks Sketchy. != 48. Exiting.."
				sys.exit()

			curr_avg_p = findAverage(temp_percents)
			curr_avg_t = findAverage(temp_times)
			curr_avg_i = findAverage(temp_iterations)

			percents.append(curr_avg_p)
			times.append(curr_avg_t)
			iterations.append(curr_avg_i)
			#print curr_avg_p
			#print curr_avg_t
			#print curr_avg_i

	print ("Number of Files " + str(len(files)))
	max_index_p = findMaxIndex(percents)
	min_index_t = findMinIndex(times)
	min_index_i = findMinIndex(iterations)

	print ("Best Sat. Clauses: \t "+ files[max_index_p] + "\twith " + str(percents[max_index_p]) + "%")
	print ("Least Avg. Time: \t " + files[min_index_t] + "\twith " + str(times[min_index_t]) + "s")
	print ("Least Avg. Its to Best:  " + files[min_index_i] + "\twith " + str(iterations[max_index_p]) + " iterations")

analyzedata(statsdir);




