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
	errors = []

	for filename in os.listdir(sd):
		if filename.endswith(".txt"):
			files.append(filename)

			temp_errors = []

			with open(sd + "/" + filename) as f:
				next(f)
				for line in f:
					tokens = line.split(" ")

					temp_errors.append(float(tokens[3]))

			curr_avg_error = findAverage(temp_errors)
			errors.append(curr_avg_error)

	print ("Number of Files " + str(len(files)))
	min_error_index = findMinIndex(errors)

	print("Least Avg. Error:\t" + files[min_error_index] + "\twith\t" + str(errors[min_error_index]))

analyzeData(statsdir)