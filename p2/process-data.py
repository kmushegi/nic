import sys
import os

statsdir = "stats"

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

def analyzeData(sd):
	files = []
	times = []
	values = []

	for filename in os.listdir(sd):
		if filename.endswith(".txt"):
			files.append(filename)

			temp_times = []
			temp_values = []

			with open(sd+"/"+filename) as f:
				lineCount = 0
				for line in f:
					if line != "\n":
						line = line.replace("\n","")
						#print line
						lineCount += 1

						tokens = line.split(" ")
						print(tokens[0] + " " + tokens[1])
						temp_times.append(tokens[0])
						temp_values.append(tokens[1])

			curr_avg_t = findAverage(temp_times)
			curr_avg_v = findAverage(temp_values)

			times.append(curr_avg_t)
			values.append(curr_avg_v)

	print("Number of Files " + str(len(files)))
	min_time_index = findMinIndex(times)
	min_value_index = findMinIndex(values)

	print("Least Avg. Time: \t" + files[min_time_index] + "\t with " + str(times[min_time_index]) + "s")
	print("Best Avg. Value: \t" + files[min_value_index] + "\t with" + str(values[min_value_index]))

