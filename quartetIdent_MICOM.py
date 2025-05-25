import cobra
import copy
import math
import numpy as np
import os
import pandas as pd
import time

from cobra.io.mat import load_matlab_model
from micom import Community, data
from multiprocessing import Pool
from os.path import join


st = time.time()


# Calculating the community growth
def HidCooperation(modelFile, commonRichness):
	com = Community(modelFile)
	medium = com.medium
	medium.update(
		(key,value*commonRichness) for key, value in medium.items()
		)
	com.medium = medium
	soln = com.cooperative_tradeoff(fraction=1.0)
	return soln

# Calculate individual growth of microbes for varying richness
def calcMicrobeAloneGr(model, richnessVal):
	modelNew = copy.deepcopy(model)
	medium = modelNew.medium
	medium.update(
		(key, value*richnessVal) for key, value in medium.items()
		)
	modelNew.medium = medium
	solution = modelNew.optimize()
	microbeAloneGr = solution.objective_value
	return microbeAloneGr



# Call the calcMicrobeAloneGr function
def runMicrobeAloneGr(modelLocation, commonRichness, modelNums):
	model = []
	endVal = 100
	microbeAloneGr = np.zeros([2,endVal])
	richness = [3*commonRichness*val/endVal for val in range(endVal)]
	for iterVal in range(2):
		model = load_matlab_model(join(modelLocation,"models","model_{}.mat".format(modelNums[iterVal])))
		with Pool(8) as pool:
			microbeAloneGr[iterVal] = pool.starmap(calcMicrobeAloneGr, [(model, richness[iterNum]) for iterNum in range(endVal)])
	index = np.argmin([abs(richnessVal - commonRichness) for richnessVal in richness])
	return microbeAloneGr, richness, index

# Invoke the callMicrobeAlobeGr to simulate individual growth for varied richness
if __name__ == "__main__":
	commonRichness = 2

	# Create the list of models 
	modelLocation = os.getcwd()
	
	files = []
	modelNums = [1,2]
	file1 = join(modelLocation,"models","model_1.mat")
	file2 = join(modelLocation,"models","model_2.mat")
	files.append(file1)
	files.append(file2)
	ids = ["model_1", "model_2"]	
	modelFiles = pd.DataFrame({"id":ids, "file": files})

	microbeAloneGr, richness, index = runMicrobeAloneGr(modelLocation, commonRichness, modelNums)

# Call HidCooperation function to calculate the community biomass

	comSoln = HidCooperation(modelFiles, commonRichness)

	commGr = comSoln.members.growth_rate * comSoln.members.abundance

# Identifying the left-over resource
	leftOverRs = np.zeros([2])

	val1 = abs(commGr[2] - microbeAloneGr[1]).argmin()

	if commonRichness > richness[val1]:
		leftOverRs[0] = commonRichness - richness[val1]
	else:
		leftOverRs[0] = 0

	val2 = abs(commGr[1] - microbeAloneGr[0]).argmin()

	if commonRichness > richness[val2]:
		leftOverRs[1] = commonRichness - richness[val2]
	else:
		leftOverRs[1] = 0

	NetInteraction = np.zeros([2])
	netPos = np.zeros([2])
	netNeg = np.zeros([2])


	for modelVal in range(2):

		indGr = microbeAloneGr[modelVal,index]
		leftOverGr = microbeAloneGr[modelVal, np.array((abs(richness-leftOverRs[modelVal])).argmin())]
		if leftOverGr < 0:
			leftOverGr = 0
		
		if commGr[modelVal+1] !=0:

			NetInteraction[modelVal] = commGr[modelVal+1] - indGr
			netPos[modelVal] = commGr[modelVal+1] - leftOverGr
			netNeg[modelVal] = NetInteraction[modelVal] - netPos[modelVal]
			
	Quartet = [netPos[0],netNeg[0], netPos[1], netNeg[1]]
	print('Net Interaction:',NetInteraction)
	print('Quartet:',Quartet)
	# print(indGr)
	# print(commGr)
	# print(leftOverGr)
	with open('results/output_1_2.txt','w') as f:
		print(NetInteraction, file = f)
		print(Quartet, file = f)
		# print(indGr, file = f)
		# print(commGr, file = f)
		# print(leftOverGr, file = f)

	et = time.time()
	elapsed_time = et - st
	print('Execution time:', elapsed_time, 'seconds')

