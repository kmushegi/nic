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
	
	w, h = 20, 11;
	medians = [[0 for x in range(w)] for y in range(h)] 

	for filename in os.listdir(sd):
		if filename.endswith(".txt"):
			files.append(filename)

			print(filename)

			temp_times = []
			temp_values = []

			with open(sd+"/"+filename) as f:
				lineCount = 0
				for line in f:
					if line != "\n":
						line = line.replace("\n","")
						#print line

						tokens = line.split(" ")
						temp_times.append(tokens[0])
						temp_values.append(float(tokens[1]))

						for x in range(2, len(tokens)-1):
							medians[x-2][lineCount] = float(tokens[x])

						lineCount += 1

			counter = 0
			for i in medians:
				medianRow = i
				medianRow = sorted(medianRow)
				print("Median at " + str(counter) + ": " + str((medianRow[10] + medianRow[9])/2))
				counter+=1

			curr_avg_v = findAverage(temp_values)
			print("Average Best Value: " + str(curr_avg_v))
			temp_values = sorted(temp_values)
			print("Median Best Value : " + str((temp_values[10] + temp_values[9])/2))
			print("\n")

analyzeData(statsdir)

