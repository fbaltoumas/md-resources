#!/usr/bin/env python
from __future__ import division;
from sys import argv;
import math;


help="""
gmx_jarzynski.py: Potential of Mean Force (PMF) calculations with the Jarzynski equality for pulling simulations with GROMACS.


The Jarzynski Equality estimates PMF (DeltaF) as:

exp(DeltaF) = < k*T*exp(-W/(k*T))>,  or  DeltaF = <-k*T*ln(-W/(k*T))> 

where k is the Boltzmann constant, T is the temperature, W is the Work produced by the simulation and the '< >' denote the ensemble average.


BASIC USAGE:
   python gmx_jarzynski.py -i [--input] pull-files.txt -t [--temperature] temperature -v [--velocity] velocity [-x [--xvg]] [-h [--help]]


Options:
---------
* -i or --input: The input file.  A text containing a list of pull xvg files to be processed.  An example is shown below.
* -t or --temperature:  The simulation temperature (in K).  Default value is 300 K.
* -v or --velocity:  The pulling velocity or pulling rate (in nm/ps).  Please input the same value you used in your *mdp file.  Default value is 0.001 nm/ps.
* -x or --xvg: Enable XVG format for the output.  If not set, a simple text is produced.
* -e or --error: also calculate the error estimate.
* -h or --help: This help text.



Results are given in the standard output (STDOUT).  The first column is the reaction coordinate (in Angstroms), while the second column is the PMF curve (in kcal/mol).

FORMAT for the pull-files list text.

The list should be a two column, tab or space separated file.  The first column should have the force vs time xvg file (pullf.xvg) and the second should have the extension vs time xvg file (pullx.xvg).  An example for a set of ten runs (run 0 to run 9) is shown below:

run0_pullf.xvg    run0_pullx.xvg
run1_pullf.xvg    run1_pullx.xvg
run2_pullf.xvg    run2_pullx.xvg
run3_pullf.xvg    run3_pullx.xvg
run4_pullf.xvg    run4_pullx.xvg
run5_pullf.xvg    run5_pullx.xvg
run6_pullf.xvg    run6_pullx.xvg
run7_pullf.xvg    run7_pullx.xvg
run8_pullf.xvg    run8_pullx.xvg
run9_pullf.xvg    run9_pullx.xvg


ATTENTION: DO NOT leave a final blank line in the text, or the script will fail.
""";


temp=300;
vel=0.001;
xvg=False;
error=False;

if len(argv)==1:
	print """
gmx_jarzynski.py: Potential of Mean Force (PMF) calculations with the Jarzynski equality for pulling simulations with GROMACS.

No arguments given.  Use '-h' or '--help' to display the help text.
	""";
	exit();


for i in range(len(argv)):
	if argv[i]=="-i" or argv[i]=="--input":
		pull_files=argv[i+1].rstrip("\n");
	if argv[i]=="-t" or argv[i]=="--temperature":
		temp=float(argv[i+1]);
	if argv[i]=="-v" or argv[i]=="--velocity":
		vel=float(argv[i+1]);
	if argv[i]=="-x" or argv[i]=="--xvg":
		xvg=True;
	if argv[i]=="-w" or argv[i]=="--work":
		printWork=True;
	if argv[i]=="-e" or argv[i]=="--error":
		error=True;
	if argv[i]=="-h" or argv[i]=="--help":
		print help;
		exit();



class Work():
	def __init__(self):
		self.velocity=float();
		self.force=list();
		self.time=list();
		self.work=list();
		self.extension=list();
		self.dt=float();
	def calcWork(self):
		#print "to be done";
		dt=0;
		for i in range(0, len(self.time)-1, 1):
			dt+=self.time[i+1]-self.time[i];
		self.dt=dt/len(self.time);
		#print dt;
		fsum=0;
		for f in self.force:
			fsum+=(f/4.184/10)*self.velocity*self.dt;  #conversion from kj/mol/nm to kcal/mol/angstroms
			self.work.append(fsum);		

class XVG():
	def __init__(self):
		self.time=list();
		self.measurement=list();
	def parseXvg(self,xvg):
		#print xvg;
		inp=open(xvg, "r");
		for line in inp:
			if line[0] not in ["#", "@"]:
				spl=line.split();
				self.time.append(float(spl[0]));
				self.measurement.append(float(spl[1]));
		inp.close();


class PMF():
	def __init__(self):
		self.pmf=list();
		self.stdev=list();
		self.stderr=list();
		self.coord=list();
	def calcJarzynski(self,workList, temperature):
		pmf=list();
		z_coord=list();
		boltzmann=0.0019872041; #boltzmann constant in kcal/(mol*K)
		kT=boltzmann*temperature;

		for i in range(len(workList[0].time)):
			#print i;
			Fexp=0;
			z=0;
			stdev=0;
			for work in workList:
				Fexp+=math.exp(-1*work.work[i]/kT);
				z+=work.extension[i];
			Fexp=Fexp/len(workList); #mean
			z=z/len(workList); #mean
			z=z*10; #conversion to angstroms;
			if Fexp!=0:
				fes=-1*kT*math.log(Fexp);
			else:
				fes=0;
			self.coord.append(z);
			self.pmf.append(fes);
			for work in workList:
				stdev+=(Fexp-(-1*kT*math.log(math.exp(-1*work.work[i]/kT))))**2;
			stdev=math.sqrt(stdev/(len(workList)-1));
			stderr=stdev/math.sqrt(len(workList));
			self.stdev.append(stdev);
			self.stderr.append(stderr);
				
	




workList=list();
file_list=open(pull_files, "r");
for line in file_list:
	files=line.split();
	pullf=XVG();
	pullf.parseXvg(files[0]);
	pullx=XVG();
	pullx.parseXvg(files[1]);
	work=Work();
	work.time=pullx.time;
	work.velocity=vel;
	work.force=pullf.measurement;
	work.extension=pullx.measurement;
	work.calcWork();
	workList.append(work);
file_list.close();



pmf=PMF();
pmf.calcJarzynski(workList, temp);		



#finally, print the output

if xvg==True:
	print """@    title "Potential of Mean Force"
@    xaxis  label "\\xz\\f{} (\\cE\\C)"
@    yaxis  label "\\xD\\f{}F (kcal/mol)"
@ legend on
@ legend box on
@ legend loctype view
@ legend 0.78, 0.8
@ legend length 2
@ s0 legend "Jarzynski PMF" """;
	if error==True:	
		print "@TYPE xydy";
	else:
		print "@TYPE xy";

for i in range(len(pmf.pmf)):
	if error==True:
		print "%f\t%f\t%f" %(pmf.coord[i], pmf.pmf[i], pmf.stderr[i]);
	else:
		print "%f\t%f" %(pmf.coord[i], pmf.pmf[i]);

