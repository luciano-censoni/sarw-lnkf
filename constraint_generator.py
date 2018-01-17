# AUTHOR: Luciano Censoni, 2018
# Institute of Chemistry - University of Campinas
# luciano.censoni@gmail.com

# The routine defined in this module generates a
# ".dat" constraint file with a set of geometric
# constraints corresponding to all the CA-CA
# distances in the given pdb.

# A few of the arguments and behaviours are due to
# legacy code requirements.


def generate_constraints(pdb_file, constrn_k=None, padding=1.0):
  from math import sqrt
  import mypdb as pdb

  constrn_k = str(constrn_k) #legacy
  p = pdb.PDB(pdb_file)

  lista = filter( lambda item: item[1] == 'CA', p.get_atomos_lista())
  num_c_alpha = len(lista)

  restr = []
  ind = -1
  for i in range(num_c_alpha):
    for j in range(i+1, num_c_alpha):
      #all pairs of CA
      #x,y,z: 3,4,5
      #6: resid, 7: resname
      #0: atomid
      x1, x2 = lista[i][3], lista[j][3]
      y1, y2 = lista[i][4], lista[j][4]
      z1, z2 = lista[i][5], lista[j][5]

      act_x = (x2-x1)
      act_y = (y2-y1)
      act_z = (z2-z1)

      dist  = sqrt( act_x**2 + act_y**2 + act_z**2 )
      cycle = lista[j][6] - lista[i][6] #length of cycle to be closed; always positive

      ind += 1
      restr.append([lista[i][0], lista[j][0], dist-padding, dist+padding, "#True", ind, cycle, act_x, act_y, act_z])
  #i,j

  _file_restr = pdb_file.rpartition(".")[0]+".dat"
  file_restr = open(_file_restr, "w")

  if restr == []: arq_restr.write("# EMPTY FILE") #legacy
  else:
    file_restr.write("# total = " + str(len(restr)) + "\n")
    file_restr.write("#ind1    ind2    dmin        dmax        k    bool    rest_index    loop   act_x        act_y        act_z\n")

    for item in restr:
      file_restr.write(str(item[0]) + '    ' + str(item[1]) + '    ' + str(item[2]) + '    ' + str(item[3]) + '    ' + constrn_k + '    ' + item[4] + '    ' + str(item[5]) + '    ' + str(item[6]) + '    ' + str(item[7]) + '    ' + str(item[8]) + '    ' + str(item[9]) + '\n') #reduce( lambda x, y: x+'    '+y, item)

  file_restr.close()
  return _file_restr
