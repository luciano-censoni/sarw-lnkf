# AUTHOR: Luciano Censoni, 2018.
# Institute of Chemistry - University of Campinas
# luciano.censoni@gmail.com

# This script calculates the self-information content of each constraint
# based on a self-avoiding random walk (SARW) polymer model.
# The input format is the one defined in constraint_generator.py


from sys import argv, exit, stdout
from numpy import arange, exp, log2
from collections import Counter



write = stdout.write
def progress_bar(index, length, est_index, est_length, fill):
    write("\r|")
    for k in range(100):
      if k/100.0 <= est_index/est_length: write(fill)
      else: write("-")
    write("|i= " + str(index) + " len= " + str(length) +" {0:.2f}".format(100.0*est_index/est_length)+"%")


def flatten(lista):
  result = []
  extend = result.extend
  for sub in lista:
      extend(sub)
  return result


def saw_expectation(loop1, loop2, dlow, dupp):
  #likelihood assessed by SAW polymer model
  #zero to infinity is guaranteed to sum to 1.0
  #probability is the integral from dmin to dmax
  step = 0.01
  x = arange(dlow, dupp, step)
  two_nu = 1.15
  l2nu = abs(loop1 - loop2) **two_nu
  Al2nu = 16.58 * l2nu
  c = 2.4343 / (Al2nu**1.5)
  f = c * exp( -1.0 * ( (0.8555/Al2nu) * x**2.0  )**1.6 ) * x**2.0
  prob = sum(f*step)
  return -1.0 * log2(prob) #self-information of observed event


#I/O
try:
  restr_file = argv[1]
except:
  restr_file = raw_input("Input constraint list file: ")

try:
  pdb_file = argv[2]
except:
  pdb_file = raw_input("Input pdb file: ")

try:
  redund_dist = float(argv[3])
except:
  redund_dist = float( raw_input("Input residue redundancy distance: "))

try:
  min_represent = int(argv[4])
except:
  min_represent = int(raw_input("Input minimum residue representation in restriction set: "))

try:
  range_dmin = float(argv[5])
except:
  range_dmin = float(raw_input("Input minimum distance for restriction consideration: "))

try:
  range_dmax = float(argv[6])
except:
  range_dmax = float(raw_input("Input maximum distance for restriction consideration: "))
assert range_dmin <= range_dmax
dist_range = [range_dmin, range_dmax]
del range_dmin, range_dmax

try:
  monomer = float(argv[7])
except:
  monomer = float(raw_input("Input monomer diameter for information calculation: "))
monomer14 = monomer**1.4 #CURRENTLY UNUSED

FILTER = ("filter" in argv)


try:
  with open(restr_file) as f:
    restr_set = []
    for line in f:
      if not line.strip() or line[0]=="#": continue
      l = line.split()
      restr_set.append( [ int(l[0]), int(l[1]), float(l[2]), float(l[3]), float(l[8]), float(l[9]), float(l[10]) ] ) #atom1, atom2, dmin, dmax, act_x, act_y, act_z
except:
  print "Problem reading", restr_file, ". Terminating."
  exit()

try:
  with open(pdb_file) as f:
    residues = {}
    last_res = None
    for line in f:
      if not line.strip(): continue
      l = line.split()
      if l[0]=="TER": break
      if not (l[0]=="ATOM" or l[0]=="HETATM"): continue
      if not (line[13:15] == "CA"): continue #residues without CA ignored
      this_res = int(line[22:26])
      if last_res is None: last_res = this_res #first accepted residue
      else: assert (this_res==last_res or this_res==(last_res+1)) #no jumps!
      residues[ int(line[6:11]) ] = int(line[22:26]) #atom index, residue index
      last_res = this_res
except Exception as e:
  print "Problem reading", pdb_file, ". Terminating.", type(e)
  exit()


#PROCESSING
for i, restr in enumerate(restr_set):
  restr_set[i].append( saw_expectation( residues[ restr[0] ], residues[ restr[1] ], restr[2], restr[3]))


#BEGIN PRUNING; criterion = higher info given distance range
initial_num_restr = len(restr_set)
if FILTER:
  print "\nBEGINNING GREEDY LOW_INFO FILTER"

  print "begin filter len(restr_set)=", initial_num_restr
  resrepresent = dict( Counter(                                         #count
                   map( lambda atom: residues[atom],                    #transform
                     flatten(                                           #flatten
                       map( lambda restr: restr[:2], restr_set) ) ) ) ) #slice
  print "initial max resrepresent =", max(resrepresent.values()), "min =", min(resrepresent.values())
  #for testing purposes:
  assert max(resrepresent.values()) == min(resrepresent.values()) #assert restrictions are all-on-all
  assert set(residues.values()) == set(resrepresent.keys()) #making sure all pdb residues are represented at least once
                                                            #in the set of restrictions

  if dist_range[0] >= 0.0: #filter by distance range
    print "removing distant restrictions"
    restr_set = filter(lambda x: dist_range[0] <= 0.5*(x[2]+x[3]) <= dist_range[1], restr_set)
    print "after removal len(restr_set)=", len(restr_set)
    resrepresent = dict( Counter(                                         #count
                     map( lambda atom: residues[atom],                    #transform
                       flatten(                                           #flatten
                         map( lambda restr: restr[:2], restr_set) ) ) ) ) #slice
    print "current max resrepresent = ", max(resrepresent.values()), "min =", min(resrepresent.values())
    print "current length = ", len(restr_set)
    if min(resrepresent.values()) < min_represent: raise ValueError

  temp_n_contacts = len(restr_set)
  print "number of contacts:", temp_n_contacts

  for k in resrepresent.keys(): resrepresent[k] = 0
  restr_set.sort( key=lambda item: item[-1], reverse=True) #decreasing #FIRST SORT

  for i in range(len(restr_set)): #establishing min_represent
    resrepresent[ residues[ restr_set[i][0] ] ] += 1
    resrepresent[ residues[ restr_set[i][1] ] ] += 1
    if len( filter(lambda x: x < min_represent, resrepresent.values() )) > 0: continue
    del restr_set[i+1:] #deletes from here on to the end
    break

  print "current max resrepresent = ", max(resrepresent.values()), "min =", min(resrepresent.values())
  print "current length = ", len(restr_set)
  if min(resrepresent.values()) < min_represent: raise ValueError


  current_len = len(restr_set)
  print "begin pruning len(restr_set)=", current_len
  restr_set.sort( key=lambda item: item[-1]) #increasing #SECOND SORT

  i = -1
  removed = False
  while True: #lowering max_represent
    i += 1

    progress_bar(i, current_len, (i+1)**0.5, (current_len+1)**0.5, u"\u25A3")

    if i == current_len: #will try to access beyond limits
      if not removed: #no change in last pass
        write("\n")
        break
      else: #one more pass needed
        i = -1
        removed = False
        continue

    res1 = residues[ restr_set[i][0] ]
    res2 = residues[ restr_set[i][1] ]
    cnt1 = 0
    cnt2 = 0
    for j in range( i+1, current_len):
      if (residues[ restr_set[j][0] ] == res1 or residues[ restr_set[j][1] ] == res1): cnt1 += 1
      if (residues[ restr_set[j][0] ] == res2 or residues[ restr_set[j][1] ] == res2): cnt2 += 1
      if cnt1 >= min_represent and cnt2 >= min_represent:
        #minimum number of higher-info restrictions for each residue
        removed = True
        del restr_set[i]
        current_len -= 1
        i -= 1
        break
    #end while

  resrepresent = dict( Counter(                                         #count
                   map( lambda atom: residues[atom],                    #transform
                     flatten(                                           #flatten
                       map( lambda restr: restr[:2], restr_set) ) ) ) ) #slice
  print "current max resrepresent =", max(resrepresent.values()), "min =", min(resrepresent.values())
  if min(resrepresent.values()) < min_represent: raise ValueError

  print "FINISHED GREEDY LOW_INFO FILTER"

  #REDUNDANCY FILTER
  #decreasing is the output format
  restr_set.sort( key=lambda item: item[-1], reverse=True) #decreasing #THIRD SORT

  #redundancy pruning
  if redund_dist >= 0: #zero will prune only constraints about the same pair
    print "\nCommencing pairwise redundancy pruning."

    current_len = len(restr_set)
    print "len(restr_set)=", current_len

    i = -1
    while True: #removing equivalent pairs
      i += 1

      progress_bar(i, current_len, (i+1)**0.5, (current_len+1)**0.5, u"\u2588")

      if i >= current_len:
        write('\n')
        break #last i is len-1

      to_remove = []
      for j in range( i+1, current_len):
        if (0 < redund_dist < 1): #overlap
          overlap = len( set(range( residues[restr_set[i][0]], residues[restr_set[i][1]]+1 )) &
                         set(range( residues[restr_set[j][0]], residues[restr_set[j][1]]+1 )) )
          if overlap >= (redund_dist * len(set(residues.values()))):
                        # min( len(set(range( residues[restr_set[i][0]], residues[restr_set[i][1]]+1 ))),
                        #      len(set(range( residues[restr_set[j][0]], residues[restr_set[j][1]]+1 ))) ) ):
            to_remove.insert(0, j)
        else: #(redund_dist >= 1); endpoints
          delta = ( abs( residues[restr_set[j][0]] - residues[restr_set[i][0]]) +
                    abs( residues[restr_set[j][1]] - residues[restr_set[i][1]]) )
          if delta <= redund_dist:
            to_remove.insert(0, j)
      for item in to_remove:
        del restr_set[item]
        current_len -= 1

    print "Finished pairwise redundancy pruning."
  else: print "Skipping pairwise redundancy pruning."
else:
  print "SKIPPING GREEDY LOW_INFO_FILTER"
  print "Skipping pairwise redundancy pruning."
  temp_n_contacts = -1.0

#decreasing is the output format; might be repeating but there's no harm
restr_set.sort( key=lambda item: item[-1], reverse=True) #decreasing

final_num_restr = len(restr_set)
print
#I/O
try:
  with open(restr_file.rpartition(".")[0]+".sorted.dat", 'w') as f:
    f.write("#restr_frac="+str(float(final_num_restr)/float(initial_num_restr))+"   ")
    f.write("#contact_frac="+str(float(final_num_restr)/float(temp_n_contacts))+"\n")
    f.write(reduce( lambda x,y: str(x)+"     "+str(y), ("ATOM1","ATOM2","DMIN","DMAX","INFORMATION\n")))
    for restr in restr_set:
      f.write(reduce(lambda x,y: str(x)+"     "+str(y), restr))
      f.write("\n")

except:
  print "Problem writing output file. Terminating."
  exit()

print "Constraint file written successfully."
