from pymol import cmd,stored;
import numpy;
import math;

def norm(a):
	return math.sqrt(numpy.sum(a*a));
	
	
def cremer_pople(selection, state=-1, quiet='no'):
		"""
		cremer_pople
		Calculates Puckering parameters for carbohydrate rings using the
		Cremer - Pople (CP) definitions (Phi, Theta and Q).  The CP parameters are calculated
		using triangular decomposition, as described by Hill and Reilly
		in J. Chem. Inf. Model. 2007, 47(3):1031-5 (PubMed: 17367125).
		This function is based on the original implementation of the Hill-Reilly method
		introduced in the Supporting Information of the aforementioned paper.
		
		USAGE:
		cremer_pople selection, state=-1, quiet=no
		
		selection:	A selection (or object) corresponding to a carbohydrate ring.
		state:		The number of the state that is to be analyzed. Default is "-1", i.e. current state.
		quiet:		Choice of output. Default is "no", which prints the output.  If set to 'yes', no output is printed,
					but the results are returned as a list object instead.  This enables the use of 'cremer_pople' in 
					scripting.
		
		cremer_pople performs CP puckering calculations for a single state.  Setting the 'state' option to '-1'
		(default) calculates CP for the current state.  In multi-state objects/selections (e.g. NMR structures or
		Molecular Dynamics trajectories), each CP calculations for each frame can be performed simply by setting 'state'
		to the appropriate value (e.g. for frame No. 10, use state=10).
		Do NOT (!!!) set state to '0', as the script will fail.  If you want to perform calculations for all states of an object,
		use the 'cremer_pople_multi' function instead.
		
		EXAMPLES:
		
		For a selection of a ring, in a single state object, or for a multi-state object and its current state:
		cremer_pople ring
		
		For a selection of a ring in a multi-state object, for a specific state, e.g. state 15
		cremer_pople ring, state=15
		
		For saving the results in a list object, which can then be used programmatically:
		cp_list=[]
		cp_list=cremer_pople ring, quiet=yes
		
		For average CP calculations in multi-state objects, use cremer_pople_multi instead.
		"""
		atoms=numpy.zeros((6,3),dtype="float64");
		stored.atms=[];
		cmd.iterate_state(state, selection, 'stored.atms.append([name,x,y,z])');
		for i in range(len(stored.atms)):
			atoms[i]=stored.atms[i][1:4];
		center=numpy.add.reduce(atoms)/6.0;
		atoms=atoms-center;
		
		r1a=numpy.zeros((3),dtype="float64");
		r2a=numpy.zeros((3),dtype="float64");
		for j,i in enumerate(atoms[0:6]):
			r1a+=i*math.sin(2.0*math.pi*j/6);
			r2a+=i*math.cos(2.0*math.pi*j/6);

		n=numpy.cross(r1a,r2a);
		n=n/norm(n);

		z=numpy.dot(atoms,n);
		q2cosphi=0;
		q2sinphi=0;
		q1cosphi=0;
		q1sinphi=0;
		q3=0;
		bigQ=0;
		sqrt_2=math.sqrt(2);
		inv_sqrt_6=math.sqrt(1/6);

		for j,i in enumerate(z):
			q2cosphi+=i*math.cos(2*math.pi*2*j/6);
			q2sinphi-=i*math.sin(2*math.pi*2*j/6);
			q1cosphi+=i*math.cos(2*math.pi*j/6);
			q1sinphi-=i*math.sin(2*math.pi*j/6);
			q3+=i*math.cos(j*math.pi);
			bigQ+=i*i;
		q2cosphi=sqrt_2*inv_sqrt_6*q2cosphi;
		q2sinphi=sqrt_2*inv_sqrt_6*q2sinphi;
		q3 = inv_sqrt_6 * q3;
		q2=math.sqrt(q2cosphi* q2cosphi + q2sinphi*q2sinphi);
		q1=math.sqrt(q1cosphi* q1cosphi + q1sinphi*q1sinphi);
		bigQ=math.sqrt(bigQ);

		if q2cosphi>0.0:
			if q2sinphi>0.0:
				phi=math.degrees(math.atan(q2sinphi/q2cosphi));

			else:
				phi=240.0-abs(math.degrees(math.atan(q2sinphi/q2cosphi)));

		else:
	 
			if q2sinphi>0.0:
				phi=60-abs(math.degrees(math.atan(q2sinphi/q2cosphi)));
				if phi<0:
					phi=360-abs(phi);

			else:
				phi=60+abs(math.degrees(math.atan(q2sinphi/q2cosphi)));


		theta=abs(math.degrees(math.atan(q2/q3)));
		if quiet=="yes":
			return [phi, theta, bigQ];
		else:
			print("Phi Theta Q");
			print("%.3f %.3f %.3f" %(phi,theta,bigQ));


		
cmd.extend(cremer_pople);


def cremer_pople_multi(selection, quiet='yes'):
	"""
	cremer_pople_multi
	Performs Cremer-Pople (CP) puckering calculations for cabohydrate rings in multi-state objects,
	such as NMR ensembles or Molecular Dynamics trajectories.  The function essentially calls the 'cremer_pople'
	function for each state of the object and, by default, returns as a result the average values (Mean +/- Standard Error)
	of the three CP parameters (Phi, Theta and Q).
	For more details on the implementation, check the help text of the 'cremer_pople' function (type 'help cremer_pople').
	
	USAGE:
	cremer_pople_multi selection, quiet='yes'
	
	where selection is the selection (or object) of a multi-state ring structure.
	
	'quiet' determines whether the results of individual frames are printed or not.  If set to 'yes', only the final result
	(average Phi, Theta and Q)is printed.  If set to 'no', the values of each individual frame are also printed.
	
	EXAMPLES:
	Obtain the averages for a multi-state object (ring):
	cremer_pople_multi ring
	
	Print all frame values for the same object, as well as the averages:
	cremer_pople_multi ring, quiet=no
	
	This function is intended for use with mult-state objects.  For single state objects, use the 'cremer_pople' function instead.
	"""
	states=cmd.count_states(selection);
	cp=[];
	for i in range(1, states+1, 1):
		cp.append(cremer_pople(selection, state=i, quiet="yes"));
	phis=[];
	thetas=[];
	Qs=[];
	if quiet=="no":
		print("Phi Theta Q");
	for cremer in cp:
		if quiet=="no":
			print("%.3f %.3f %.3f" %(cremer[0],cremer[1],cremer[2]));
		phis.append(cremer[0]);
		thetas.append(cremer[1]);
		Qs.append(cremer[2]);
	avgPhi=numpy.mean(phis);
	stderrPhi=numpy.std(phis)/math.sqrt(states);
	avgTheta=numpy.mean(thetas);
	stderrTheta=numpy.std(thetas)/math.sqrt(states);
	avgQ=numpy.mean(Qs);
	stderrQ=numpy.std(Qs)/math.sqrt(states);
	print("\nAverage Values (mean +/- stdErr):\nAv. Phi Av. Theta Av. Q");
	print("%.3f+/-%.3f %.3f+/-%.3f %.3f+/-%.3f" %(avgPhi,stderrPhi, avgTheta,stderrTheta, avgQ,stderrQ));

cmd.extend(cremer_pople_multi);


def dihe(selection, a1, a2, a3, a4, state=-1):
	atm1=selection + " and name %s" %a1;
	atm2=selection + " and name %s" %a2;
	atm3=selection + " and name %s" %a3;
	atm4=selection + " and name %s" %a4;
	return cmd.get_dihedral(atm1,atm2,atm3,atm4, state=state);	

def aver_dihe(selection, a1, a2, a3, a4):
	states=cmd.count_states(selection);
	dihedrals=list();
	states=cmd.count_states(selection);
	for i in range(1,states+1, 1):
		dihedrals.append(dihe(selection, a1, a2, a3, a4, state=i));
	print(dihedrals);

cmd.extend(aver_dihe);

def rao_dihedrals(selection, quiet='yes'):
		"""
		rao_dihedrals
		Calculates the virtual torsion angles a1, a2 and a3 of a hexopyranose ring, using the definitions of Rao and co-workers,
		as in the book.
		"""
		stored.names=[];
		states=cmd.count_states(selection);
		print(states);
		a1=[];
		a2=[];
		a3=[];
		for i in range(1, states+1, 1):
			cmd.iterate_state(i, selection, 'stored.names.append(name)');
			a1.append(dihe(selection, stored.names[3], stored.names[4], stored.names[1], stored.names[0], state=i));
			a2.append(dihe(selection, stored.names[5], stored.names[1], stored.names[3], stored.names[4], state=i));
			a3.append(dihe(selection, stored.names[1], stored.names[3], stored.names[5], stored.names[4], state=i));
		av1=numpy.mean(a1);
		stdErr1=numpy.std(a1)/math.sqrt(states);
		av2=numpy.mean(a2);
		stdErr2=numpy.std(a2)/math.sqrt(states);
		av3=numpy.mean(a3);
		stdErr3=numpy.std(a3)/math.sqrt(states);
		if quiet=="no":
			print("a1 a2 a3");
			for i in range(len(a1)):
				print("%.3f %.3f %.3f" %(a1[i],a2[i],a3[i]));
		print("Averages(Mean +/- Standard Error):\na1 a2 a3\n %.3f+/-%.3f %.3f+/-%.3f %.3f+/-%.3f" %(av1,stdErr1, av2,stdErr2, av3,stdErr3));

cmd.extend(rao_dihedrals);




print("""CREMER - POPLE Puckering calculations.
This script implements functions for performing Cremer-Pople (CP) puckering calculations for carbohydrate rings.
The method is based on the approach presented by Hill and Reilly in the following publication:

Hill, A.D., Reilly, P.J. (2007) Puckering coordinates of monocyclic rings by triangular decomposition.
J. Chem. Inf. Model. 47(3):1031-5, doi: 10.1021/ci600492e, PubMed: 17367125 

Two functions are implemented:

cremer_pople: performs CP calculations for single-state objects, or the specific state of a multi-state object.
cremer_pople_multi: calculates CP averages for the states of a multi-state object (NMR Ensemble, Molecular Dynamics trajectory etc).

For each function, you can find more on its use by typing:
help cremer_pople
help cremer_pople_multi


Written by Fotis A. Baltoumas, as part of his PhD thesis.  Contact: fbaltoumas @ biol.uoa.gr
""");
