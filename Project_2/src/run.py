import os
import sys
from fileReader import readImage
from fileReader import readTextFile

optimalFolder = "..\Test Images GT" # you may have to specify the complete path
studentFolder = "..\Solution Images" # you may have to specify the complete path
colorValueSlackRange = 40
blackValueThreshold = 100 # colors below 100 is black
pixelRangeCheck = 4
checkEightSurroundingPixels = True

def readFilesFromFolder(directory, folderNum):
	allFiles = []
	path = directory + "\\" + folderNum
	for filename in os.listdir(path):
		if filename.endswith(".jpg") or filename.endswith(".png"):
			filename = os.path.join(path, filename)
			allFiles.append(readImage(filename))
	return allFiles


def comparePics(studentPic, optimalSegmentPic):
	# for each pixel in studentPic, compare to corresponding pixel in optimalSegmentPic 
	global colorValueSlackRange
	global checkEightSurroundingPixels
	global pixelRangeCheck
	width, height= studentPic.shape
	counter = 0 #counts the number of similar pics
	numberOfBlackPixels = 0
	for w in range(width):
		for h in range(height):
			#if any pixel nearby or at the same position is within the range, it counts as correct
			color1 = studentPic[w][h]
			color2 = optimalSegmentPic[w][h]
			if color1 < blackValueThreshold:
				#black color
				numberOfBlackPixels +=1
				if(int(color1) == int(color2)):
					counter +=1
					continue
				elif checkEightSurroundingPixels:
					#check surroundings
					correctFound = False
					for w2 in range(w-pixelRangeCheck, w + pixelRangeCheck):
						if(correctFound):
							break
						for h2 in range(h - pixelRangeCheck, h + pixelRangeCheck):
							if(w2 >=0 and h2 >= 0 and w2 < width and h2 < height):

								color2 = optimalSegmentPic[w2][h2]
								if( color1 - colorValueSlackRange< color2  and color2 < colorValueSlackRange + color1):
									correctFound = True
									counter +=1
									break
	
	return float(counter)/float(max(numberOfBlackPixels,1))


def main(argv):
	optimalFiles = readFilesFromFolder(optimalFolder, argv[0])
	studentFiles = readFilesFromFolder(studentFolder, '')
	totalScore = 0
	for student in studentFiles:
		highestScore = 0
		for opt in optimalFiles:
			result1 = comparePics(opt,student)
			result2 = comparePics(student,opt)
			result = min(result1,result2)
			highestScore = max(highestScore,result)
		totalScore += highestScore
		a = highestScore*100
		print("Score: %.2f" % a + "%")
	a = totalScore/len(studentFiles)*100
	print("Total Average Score: %.2f" % a + "%")


if __name__ == '__main__':
	print("Folder num", sys.argv[1:])
	main(sys.argv[1:])

