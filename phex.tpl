# template for beta-D-Glucopyranose
# d  - distance to d_atom
# a  - angle calculated acording "distance path"
# tn - torsion angle calculated acording "distance path"
# te - torsion exocyclic (Cn-1_Cn_On_Hn)-in disacharide one is psi ,and phi (O5_C1_O1_H1)
# ch - improper diherdal for chiral center
# o  - out of plane dihedral - ring puckering Pickett and Strauss
# C  - one or two or tree Cartesian coordinate on place of torsion angle(1), angle(2) and distance(3)

# atom	   d_atom	distance	angle  	torsion	note
  O5 	      C1               C2               C3
  C5 	      O5               C1               C2 	
  C4 	      C5               O5               C1 
  C3 	      C4               C5               O5       
  C2 	      C3               C4               C5
  C1 	      C2               C3               C4
  C6 	      C5               O5               C4
  O6 	      C6               C5               O5
  HO6	      O6               C6               C5
  O4 	      C4               C5               C3
  HO4	      O4               C4               C3
  O3 	      C3               C4               C2     
  HO3	      O3               C3               C2
  O2 	      C2               C3               C1
  HO2	      O2               C2               C1
  O1 	      C1               C2               O5
  HO1	      O1               C1               O5
