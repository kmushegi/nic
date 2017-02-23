import sys
import os

statsdir = "stats-ga"

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

	for filename in os.listdir(sd):
		if filename.endswith(".txt"):
			files.append(filename)

			temp_percents = []
			temp_times = []
			temp_iterations = []
			with open(sd+"/"+filename) as f:
				for line in f:
					if line != "\n":
						line = line.replace("\n","")

						tokens = line.split(" ")

						temp_percents.append(tokens[2])
						temp_iterations.append(tokens[3])
						temp_times.append(tokens[4])
						#print tokens[2]

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




