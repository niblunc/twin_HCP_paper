import os
import glob
import sys
import fnmatch

#basedir='/projects/niblab/scripts/milkshake'
basedir='/projects/niblab/data/eric_data/design_files/imagine'

#list_name=os.path.join(basedir, 'feat_scripts', 'list_files')
#slope_name=os.path.join(basedir,'race.txt')
ev_name=os.path.join(basedir,'EDxRate','cov_files','appEDxRating.txt')

#this file has the columns for the mat file
#0 is input
#1 is mean
#2 is appEDxRating score, demeaned
#3 is the BMI intercept, demeaned
#4 is the BMI slope, demeaned

g=open('EDxRate_input.txt','w')
h=open('mean.txt','w')
j=open('appEDxRate_cov.txt','w')
k=open('BMIintercept_cov.txt','w')
l=open('BMIslope_cov.txt','w')
m=open('interactionappBMIslope_cov.txt', 'w')
i=0 


#for dir in glob.glob('/projects/niblab/data/eric_data/W1/milkshake/level2_grace_edit/cs*.gfeat'):
#start the function 
for dir in glob.glob('/projects/niblab/data/eric_data/W1/imagine/level1_grace_edit/cs*.feat'):	
	sub0=dir.split('/')
	sub=sub0[8].strip('++.feat')
	print(sub)
#	now sub is a list of the subject numbers csXXX
#	this is saying: open the file and read as search2
	with open(ev_name, 'r') as search2:
		for line2 in search2:
			line2 = line2.split('\t')
#			this is opening the file and reading in starting at line 2, splitting at the tab
			name2 = line2[0]
#			print(name2)
			mean = line2[1]
			mean = mean.strip('\n')
			appEDxRate2 = line2[2]
			appEDxRate = appEDxRate2.strip('\n')
#			print(appEDxRate)
			BMIintercept = line2[3]
			BMIintercept = BMIintercept.strip('\n')
			BMIslope = line2[4]
			BMIslope = BMIslope.strip('\n')
			interaction = line2[5]
			interaction= interaction.strip('\n')
			print(name2+' has a combination score of '+appEDxRate+' ,BMI intercept of '+BMIintercept+' , and a BMI slope of '+BMIslope)
			if fnmatch.fnmatch(name2, sub):
				print(name2)
				print(sub)
				print(line2)
				i=i+1
				print('set feat_files('+str(i)+') '+'"'+dir+'/COPE.feat"')
				g.write('set feat_files('+str(i)+') '+'"/projects/niblab/data/eric_data/W1/imagine/level1_grace_edit/'+name2+'++.feat"'+'\n')
				h.write('set fmri(evg'+str(i)+'.1) '+mean+'\n')  
				j.write('set fmri(evg'+str(i)+'.2) '+appEDxRate+'\n')  
				k.write('set fmri(evg'+str(i)+'.3) '+BMIintercept+'\n')
				l.write('set fmri(evg'+str(i)+'.4) '+BMIslope+'\n')
				m.write('set frmi(evg'+str(i)+'.5) '+interaction+'\n')  
#				print('set fmri(evg'+str(i)+'.2) '+risk[0])
				#print('set fmri(evg'+str(i)+'.3) '+line2[4])
				#print('set fmri(evg3.'+str(i)+') '+line2[4])
				#print('set fmri(evg1.'+str(i)+') 1')
g.close()
h.close()
