# AUTHOR: Luciano Censoni, 2018
# Institute of Chemistry - University of Campinas
# luciano.censoni@gmail.com


from subprocess import check_output, call, check_call
from copy import deepcopy
from json import load, dump
from os.path import exists
from os import rename, mkdir
import sys
import fileinput
from sys import argv, stdout
def write(a):
  stdout.write("\r\033[K"+a)
  stdout.flush()
from numpy.random import choice
from matplotlib import pyplot as pp
pp.rcParams.update({'font.size': 18})


#ugly hack, I know
if len(argv) == 1: argv = ['read-folding-rate.py', 'None', '0.0', '9.5', '8', 'None', '0.01', '8.0', 'fetch', 'parse', 'calc', 'plot']


_fetch_pdb  = ("fetch"  in argv)
_parse_data = ("parse"  in argv)
_calc_info  = ("calc"   in argv)
_plot       = ("plot"   in argv)


try: min_represent = int(argv[1])
except: min_represent = None #legacy; unused

try: range_dmin = float(argv[2])
except: range_dmin = 0.0

try: range_dmax = float(argv[3])
except: range_dmax = 9.5

try: redund_dist = float(argv[4])
except: redund_dist = 8

try: min_loop = int(argv[5]) #legacy; unused
except: min_loop = None

try: padding  = float(argv[6])
except: padding = 0.01

try: monomer = float(argv[7])
except: monomer = 8.0

run = lambda x: check_output(x, shell=True)
path = lambda st: "./no-fragment/"+st


if _parse_data:
  #will parse and save
  with open("proteins_17_Dec_14.csv") as f: fratedata = f.readlines()

  data = []
  skip = False
  for line in fratedata:
    l = line.split(",")
    prot = l[0]
    if prot == "PDB Id": continue
    for linha in data:
      if linha['id'] == prot: #this id has been seen already
        skip = True
        break
    if skip:
      skip = False
      continue
    #new protein
    try:
      thisp = {}
      thisp['id']       = prot
      thisp['length']   = int(l[3])
      thisp['aresid']   = 0.0 if l[9] == "NA" else (float(l[9])/thisp['length'])
      thisp['bresid']   = 0.0 if l[11] == "NA" else (float(l[11])/thisp['length'])
      thisp['ftype']    = l[17]
      thisp['corder']   = float(l[18])
      thisp['frate']    = float(l[20])
      thisp['frgstart'] = l[23]
      thisp['frgend']   = l[24]
      thisp['consist_length'] = True
      if thisp['frgstart'] != "NA":
        if thisp['length'] != int(thisp['frgend']) - int(thisp['frgstart']) + 1:
          #print "INCONSISTENT LENGTH!"
          thisp['consist_length'] = False
      print 'Finished parsing protein', thisp['id']
      data.append(deepcopy(thisp))
    except:
      print "Failed parsing protein", prot
      continue

  del fratedata, thisp
  with open("formatted-data.json","w") as f: dump(data, f)
else:
  #just load json data
  with open("formatted-data.json") as f: data = load(f)

#filter out problematic pdb files
filt_data = []
if _fetch_pdb: mkdir(path(''))
for prot in data:
  if not prot['consist_length']: continue
  if prot['frgstart'] != "NA": continue #this second if probably contains the previous
  if prot['id'] in ['1AVZ', '1B9C', '1FMK', '1JON', '1L8W', '1M9S', '2BLM']: continue #bad pdb files
  #1avz: unexplained gap 148-179   #1b9c: unexplained gap 64-66-68   #1fmk: unexplained gap 409-424
  #1jon: unexplained gap 301-307   #l8w: unexplained gap 92-113   #1m9s: unexplained gap 319-391
  #2blm: unexplained gap 83-86
  if prot['id'] in ['1CUN','1DIV','1HRC','1LMB','1PGB','1PHP','1QOP','1RA9','3CHY']: continue #duplicated in the database

  filt_data.append(deepcopy(prot))
  if _fetch_pdb:
    run("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"+prot['id'].lower()+".ent.gz")
    run("gunzip pdb"+prot['id'].lower()+".ent.gz")
    rename("pdb"+prot['id'].lower()+".ent", path(prot['id'])+".pdb")
    #automated processing for proteins with alternate conformations and leading sequences:
    #1PNJ, 1BF4, 1BNZ, 1URN (leading sequences)
    #1YCC, 1YEA (leading sequences and single heteroresidue at position 72)
    if prot['id'] == '1PNJ':
      for line_number, line in enumerate(fileinput.input(path("1PNJ.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
        elif line_number < 202: continue
        else: #from line 202 onwards
          sys.stdout.write(line)

    if prot['id'] == '1BF4':
      for line_number, line in enumerate(fileinput.input(path("1BF4.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
        elif line_number < 636: continue
        else: #from line 636 onwards
          sys.stdout.write(line)

    if prot['id'] == '1BNZ':
      for line_number, line in enumerate(fileinput.input(path("1BNZ.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
        elif line_number < 643: continue
        else: #from line 643 onwards
          sys.stdout.write(line)

    if prot['id'] == '1URN':
      for line_number, line in enumerate(fileinput.input(path("1URN.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
        elif line_number < 1756: continue
        else: #from line 1756 onwards
          sys.stdout.write(line)

    if prot['id'] == '1YCC':
      for line_number, line in enumerate(fileinput.input(path("1YCC.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
          sys.stdout.write("REMARK HETATM transformed residue 72\n")
        elif line_number < 472: continue
        else: #from line 472 onwards
          if line[:6] == "HETATM":
            sys.stdout.write("ATOM  " ) #fix for residue 72
            sys.stdout.write( line[6:])
          else: sys.stdout.write(line)

    if prot['id'] == '1YEA':
      for line_number, line in enumerate(fileinput.input(path("1YEA.pdb"), inplace=1)):
        if line_number == 0:
          sys.stdout.write("REMARK modified protein file with leading sequences removed\n")
          sys.stdout.write("REMARK HETATM transformed residue 72\n")
        elif line_number < 439: continue
        else: #from line 439 onwards
          if line[:6] == "HETATM":
            sys.stdout.write("ATOM  " ) #fix for residue 72
            sys.stdout.write( line[6:])
          else: sys.stdout.write(line)
#end for prot in data

del data #this should probably be dumped one more time
files = run("ls "+path("")+" | grep pdb").split()
for prot in filt_data: assert prot['id']+".pdb" in files
#if the assertion passes then we have structures for all proteins in filt_data
del files

if _calc_info:
  from constraint_generator import generate_constraints
  to_remove = []
  call("ls " + path("*.dat") + " | xargs rm", shell=True) #clean-up!

  for linha in filt_data:
    pid = linha['id']
    try:
      generate_constraints(path(pid)+".pdb", padding=padding)
      if not exists( path(pid)+".dat"): raise IOError
      with open(path(pid)+".dat") as f:
        if "EMPTY FILE" in f.readline():
          remove(path(pid)+".dat")
          raise IOError
      print "Constraints generated for", pid
    except:
      to_remove.append(pid)
      print "Constraint generation failed for", pid
      continue

    for list_min_rep in [2,3,4]:
      try:
        check_call("python constrn_information.py "+ path(pid)+".dat "+ path(pid)+".pdb " + " ".join([str(redund_dist), str(list_min_rep), str(range_dmin), str(range_dmax), str(monomer), "filter"]), shell=True) #run
        if not exists( path(pid)+".sorted.dat"): raise IOError #probably will never encounter this clause
        rename( path(pid)+".sorted.dat", path(pid)+"-filter.sorted.dat")

        check_call("python constrn_information.py "+ path(pid)+".dat "+ path(pid)+".pdb " + " ".join([str(redund_dist), str(list_min_rep), str(range_dmin), str(range_dmax), str(monomer), ""]), shell=True) #run nofilter
        if not exists( path(pid)+".sorted.dat"): raise IOError #probably never will encounter this clause
        rename( path(pid)+".sorted.dat", path(pid)+"-nofilter.sorted.dat")

        print "Information calculated for", pid, "list_min_rep =", list_min_rep
      except Exception as e:
        print "Error calculating information for", pid, "list_min_rep =", list_min_rep
        print '\x1b[5;37;41m' + str(type(e)) + '\x1b[0m'
        break #breaks out of for list_
    #end for list_min_rep

    with open(path(pid)+"-nofilter.sorted.dat") as info_f: info = info_f.readlines()
    assert "frac" in info[0]
    del info[0:2] #legacy; frac = 1.0 always
    tot_inf = sum( map( lambda item: float(item.split()[-1]), info) ) / len(info)
    print "avg information for", pid, "=", tot_inf, '\n'
    linha['avinf'] = tot_inf

    with open(path(pid)+"-filter.sorted.dat") as info_f: info = info_f.readlines()
    assert "frac" in info[0]
    linha['restr_frac']   = float( info[0].split()[0].rpartition("=")[2] )
    linha['contact_frac'] = float( info[0].split()[1].rpartition("=")[2] )
    del info[0:2]
    tot_inf = sum( map( lambda item: float(item.split()[-1]), info) ) / len(info)
    print "avg information for", pid, "=", tot_inf, '\n'
    linha['favinf'] = tot_inf

  print "removing:"
  print to_remove
  for item in to_remove:
    for prot in filt_data:
      if prot['id'] == item:
        filt_data.remove( prot)
        print "removido", item
        break #breaks out of for prot

  print "Final size:", len(filt_data)
  with open("filtered-data.json","w") as f: dump(filt_data, f)
  print "Information data saved."
else: #do not calculate information
  with open("filtered-data.json") as f: filt_data = load(f)
  print "Information data loaded."

print "Calculating correlations."
from scipy.stats import pearsonr


ftype  = map( lambda x: x['ftype'],        filt_data)
co     = map( lambda x: x['corder'],       filt_data)
fr     = map( lambda x: x['frate'],        filt_data)
inf    = map( lambda x: x['avinf'],        filt_data)
finf   = map( lambda x: x['favinf'],       filt_data)
rfrac  = map( lambda x: x['restr_frac'],   filt_data)
cfrac  = map( lambda x: x['contact_frac'], filt_data)
aresid = map( lambda x: x['aresid'],       filt_data)
bresid = map( lambda x: x['bresid'],       filt_data)

cotwo     = []
frtwo     = []
inftwo    = []
finftwo   = []
comulti   = []
frmulti   = []
infmulti  = []
finfmulti = []

tworef = []
for i in range(len(ftype)):
  if ftype[i] == "Two":
    cotwo.append(  co[i] )
    inftwo.append( inf[i] )
    finftwo.append( finf[i] )
    frtwo.append(  fr[i] )
    tworef.append( i)
  else: #"Multi"
    comulti.append(  co[i] )
    infmulti.append( inf[i] )
    finfmulti.append( finf[i] )
    frmulti.append(  fr[i] )

rfrac = sum(rfrac)/(1.0* len(rfrac))
cfrac = sum(cfrac)/(1.0* len(cfrac))


if _plot:
  print 'bootstrapping co and inf'
  co_boot = []
  inf_boot = []
  finf_boot = []
  cotwo_boot = []
  inftwo_boot = []
  finftwo_boot = []
  comulti_boot = []
  infmulti_boot = []
  finfmulti_boot = []
  bootstrapN = 1000000
  for k in range(bootstrapN):
    write("k = "+ str(k))
    temp_co = []
    temp_inf = []
    temp_finf = []
    temp_fr = []
    temp_cotwo = []
    temp_inftwo = []
    temp_finftwo = []
    temp_frtwo = []
    temp_comulti = []
    temp_infmulti = []
    temp_finfmulti = []
    temp_frmulti = []
    temp = choice( len(co), len(co), replace=True) #this is numpy.random.choice
    temptwo = choice( len(cotwo), len(cotwo), replace=True)
    tempmulti = choice( len(comulti), len(comulti), replace=True)
    for i in temp:
      temp_co.append( co[i])
      temp_inf.append( inf[i])
      temp_finf.append( finf[i] )
      temp_fr.append( -1.0 * fr[i])
    co_boot.append( pearsonr(temp_co, temp_fr)[0])
    inf_boot.append( pearsonr(temp_inf, temp_fr)[0])
    finf_boot.append( pearsonr(temp_finf, temp_fr)[0])
    for i in temptwo:
      temp_cotwo.append( cotwo[i] )
      temp_inftwo.append( inftwo[i] )
      temp_finftwo.append( finftwo[i] )
      temp_frtwo.append( -1.0 * frtwo[i] )
    cotwo_boot.append( pearsonr(temp_cotwo, temp_frtwo)[0])
    inftwo_boot.append( pearsonr(temp_inftwo, temp_frtwo)[0])
    finftwo_boot.append( pearsonr(temp_finftwo, temp_frtwo)[0])
    for i in tempmulti:
      temp_comulti.append( comulti[i] )
      temp_infmulti.append( infmulti[i] )
      temp_finfmulti.append( finfmulti[i] )
      temp_frmulti.append( -1.0 * frmulti[i] )
    comulti_boot.append( pearsonr(temp_comulti, temp_frmulti)[0])
    infmulti_boot.append( pearsonr(temp_infmulti, temp_frmulti)[0])
    finfmulti_boot.append( pearsonr(temp_finfmulti, temp_frmulti)[0])
  write("k = "+str(k)+'\n')
  co_boot.sort()
  inf_boot.sort()
  finf_boot.sort()
  cotwo_boot.sort()
  inftwo_boot.sort()
  finftwo_boot.sort()
  comulti_boot.sort()
  infmulti_boot.sort()
  finfmulti_boot.sort()


  pp.hist(co_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(co, fr)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(inf_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(inf, fr)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-mixed-nofilter.pdf"); pp.clf()

  pp.hist(co_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(co, fr)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(finf_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(finf, fr)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-mixed-filter.pdf"); pp.clf()


  pp.hist(cotwo_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(cotwo, frtwo)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(inftwo_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(inftwo, frtwo)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-two-nofilter.pdf"); pp.clf()

  pp.hist(cotwo_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(cotwo, frtwo)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(finftwo_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(finftwo, frtwo)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-two-filter.pdf"); pp.clf()


  pp.hist(comulti_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(comulti, frmulti)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(infmulti_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(infmulti, frmulti)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-multi-nofilter.pdf"); pp.clf()

  pp.hist(comulti_boot, 200, histtype="step", color='k', linestyle='dotted')
  pp.xlim((-0.02, 1.02))
  pp.plot(-1.0*pearsonr(comulti, frmulti)[0], 0.0, 'ko', clip_on=False, zorder=10)
  pp.hist(finfmulti_boot, 200, histtype="step", color='r')
  pp.plot(-1.0*pearsonr(finfmulti, frmulti)[0], 0.0, 'rs', clip_on=False, zorder=10)
  pp.grid()
  pp.xlabel(r"$\rho$")
  pp.savefig("bootstrap-multi-filter.pdf"); pp.clf()

  #plot bar charts
  print "CORRELATIONS IN THE BAR CHARTS"
  barwidth = 0.20
  barindices = [1.0, 2.0, 3.0]
  _, ax = pp.subplots()
  final_rho_co    =  map(lambda x, y: -1.0* pearsonr(x, y)[0], [co, cotwo, comulti], [fr, frtwo, frmulti])
  print "co: ", final_rho_co
  conf_range_co   = [map(lambda x, y: y - x[int(0.025*bootstrapN)-1], [co_boot, cotwo_boot, comulti_boot], final_rho_co),
                     map(lambda x, y: x[int(0.975*bootstrapN)-1] - y, [co_boot, cotwo_boot, comulti_boot], final_rho_co)]
  print "co range:", conf_range_co

  final_rho_inf   =  map(lambda x, y: -1.0* pearsonr(x, y)[0], [inf, inftwo, infmulti], [fr, frtwo, frmulti])
  print "inf: ", final_rho_inf
  conf_range_inf  = [map(lambda x, y: y - x[int(0.025*bootstrapN)-1], [inf_boot, inftwo_boot, infmulti_boot], final_rho_inf),
                     map(lambda x, y: x[int(0.975*bootstrapN)-1] - y, [inf_boot, inftwo_boot, infmulti_boot], final_rho_inf)]
  print "inf range:", conf_range_inf

  final_rho_finf  =  map(lambda x, y: -1.0* pearsonr(x, y)[0], [finf, finftwo, finfmulti], [fr, frtwo, frmulti])
  print "finf: ", final_rho_finf
  conf_range_finf = [map(lambda x, y: y - x[int(0.025*bootstrapN)-1], [finf_boot, finftwo_boot, finfmulti_boot], final_rho_inf),
                     map(lambda x, y: x[int(0.975*bootstrapN)-1] - y, [finf_boot, finftwo_boot, finfmulti_boot], final_rho_inf)]
  print "finf range:", conf_range_finf

  rects1 = ax.bar(barindices, final_rho_co, barwidth, color='b', yerr=conf_range_co, error_kw=dict(ecolor='k', lw=2, capsize=5, capthick=2) )
  rects2 = ax.bar( [x+barwidth for x in barindices], final_rho_inf, barwidth, color='y', yerr=conf_range_inf, hatch="////", error_kw=dict(ecolor='k', lw=2, capsize=5, capthick=2) )
  rects3 = ax.bar( [x+(2*barwidth) for x in barindices], final_rho_finf, barwidth, color='r', yerr=conf_range_finf, hatch="xxxx", error_kw=dict(ecolor='k', lw=2, capsize=5, capthick=2) )

  ax.set_ylabel("Corr. Coefficient")
  ax.set_xticks([x+barwidth for x in barindices])
  ax.set_xticklabels(("Full Set", "Two-state", "Multistate"))
  ax.set_xlim((1.0-0.2, 3.0+3.0*barwidth+0.2))
  ax.legend((rects1[0], rects2[0], rects3[0]), ("CO", "I", "I$_r$"), loc='lower left')
  pp.grid()
  pp.savefig("barchart-both.pdf"); pp.clf()
#end if _plot

from scipy.stats import linregress
mreg, breg, rreg, preg, ereg      = linregress(inf, fr)
fmreg, fbreg, frreg, fpreg, fereg = linregress(finf, fr)

if _plot:
  _ = pp.figure(figsize=(14,6))
  pp.rcParams.update({'font.size': 18})
  ax = pp.subplot(121)
  pp.scatter(cotwo, frtwo, color="r", marker="o")
  pp.scatter(comulti, frmulti, color="b", marker="^")
  pp.xlabel("Contact Order")
  pp.ylabel("$ln(k_f)$")
  pp.grid()
  ax.yaxis.set_label_coords(-0.055555, 0.5)
  ax = pp.subplot(122)
  pp.scatter(inftwo, frtwo, color="r", marker="o")
  pp.scatter(infmulti, frmulti, color="b", marker="^")
  pp.xlabel("Avg. Topological Information")
  pp.ylabel("$ln(k_f)$")
  pp.grid()
  ax.yaxis.set_label_coords(-0.055555, 0.5)
  pp.savefig("scatterplot-nofilter.pdf"); pp.clf()

  _ = pp.figure(figsize=(6,6))
  pp.rcParams.update({'font.size': 18})
  ax = pp.subplot(111)
  pp.scatter(finftwo, frtwo, color="r", marker="o")
  pp.scatter(finfmulti, frmulti, color="b", marker="^")
  pp.xlabel("Red. Topological Information")
  pp.ylabel("$ln (k_f)$")
  pp.grid()
  ax.yaxis.set_label_coords(-0.08, 0.5)
  pp.savefig("scatterplot-filter.pdf"); pp.clf()

  pp.close("all")


dev = []
devr = []
fdev = []
fdevr = []
delta = []
deltar = []
gain = []
gainr = []
for i in range(len(inf)):
  dev.append(   abs( fr[i] - (breg  + mreg  * inf[i]) ) )
  devr.append(  abs( fr[i] - (breg  + mreg  * inf[i]) ) /abs(fr[i]) )
  fdev.append(  abs( fr[i] - (fbreg + fmreg * finf[i]) ) )
  fdevr.append( abs( fr[i] - (fbreg + fmreg * finf[i]) ) /abs(fr[i]) )
  delta.append(  dev[-1] - fdev[-1] )
  deltar.append( devr[-1] - fdevr[-1] )
  gain.append(  (finf[i] - inf[i]) )
  gainr.append( (finf[i] - inf[i]) /inf[i] )


print "DELTA I"
print "ALPHA:", pearsonr(aresid, gain)[0], "BETA:", pearsonr(bresid, gain)[0]
print "RELATIVE:"
print "ALPHA:", pearsonr(aresid, gainr)[0], "BETA:", pearsonr(bresid, gainr)[0]


#pp.subplot(121)
#pp.scatter(aresid, gain)
#pp.subplot(122)
#pp.scatter(bresid, gain)
#pp.show()