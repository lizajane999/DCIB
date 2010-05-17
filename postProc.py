#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Analyse output from dcib """
"""Lisa Miller 5/2010 """
###############################################################################
# Program information
###############################################################################
__author__ = "Lisa Miller"
__date__ = "1 May 2010"
__version__ = "$Revision: 1.0$"
__credits__ = """ """

###############################################################################
# Imports
###############################################################################
import re
from operator import itemgetter, attrgetter
###############################################################################
# Main Program
###############################################################################

if __name__ == '__main__':
  
  docInd = {}
  docCats = {}
  pdgc = {}
  clustAsgn = {}
  catAsgn = {} # will list cluster(s) assigning this cat
  wordInd = {}
  clustWords = {}

  i = 0
  with open("data/docindices.txt") as ind_file:
    for index in ind_file:
      index = index.strip()
      docInd[str(i)] = index
      i += 1
      
  print "indices: "
  print docInd
  
  with open("data/categories.txt") as cat_file:
    for line in cat_file:
      S = line.partition("\t") # take the doc id off
      line = S[2]
      line = re.sub(r'\s+', r' ', line)
      words = line.split()
      docCats[S[0]] = words
      
  print "categories:"
  print docCats
  
 # make a list of the category names
  for d,c in docCats.iteritems():
    for cat in c:
      if not catAsgn.has_key(cat):
	catAsgn[cat] = []
  
  i = 0
  with open("assignments.txt") as asgn_file:
    for line in asgn_file:
      S = line.partition(":") # take the doc id off
      probs = S[2].split()
      key = docInd[S[0]]
      pdgc[key] = probs
      
  print "assignments:"
  print pdgc
  
  # now look at assignment probs
  for doc in pdgc:
    p = 0.0
    i = 0
    allEqual = True
    for prob in pdgc[doc]:
      if float(prob) > p:
	p = float(prob)
	index = str(i)
	if i > 0:
	  allEqual = False
      if float(prob) < p:
	allEqual = False	
      i += 1
    if allEqual == False:
      if not clustAsgn.has_key(index):
	clustAsgn[index] = []
      clustAsgn[index].append(doc)
    else:
      if not clustAsgn.has_key("equal"):
	clustAsgn["equal"] = []
      clustAsgn["equal"].append(doc)
  print "+++++++++++++++++++++++++++++++++++++++++++++++"
  print clustAsgn
  print "==============================================="
  fn = "clustAsgn.txt"
  asgn_file = open(fn, "w")
  asgn_file.write("*** Documents assigned to clusters, with categories ***\n")
  for cat in clustAsgn:
    print "cluster:"+str(cat)
    asgn_file.write("cluster:"+str(cat) + '\n')
    for doc in sorted(clustAsgn[cat]):
      print str(doc) + " " + str(docCats[doc])
      asgn_file.write(str(doc))
      asgn_file.write("\t")
      for c in docCats[doc]:
	catAsgn[c].append(cat)
	asgn_file.write(str(c))
	asgn_file.write("\t")
      asgn_file.write('\n')
    asgn_file.write('\n\n')

  
 # print catAsgn  
  
  asgn_file.write("*** Categories with # of appearances in clusters ***\n")
  for cat in catAsgn:
    count = 0
    temp = []
    t2 = []
    asgn_file.write(cat + '\n')
    for k in sorted(catAsgn[cat]):
      if k not in temp:
	temp.append(k)
	c = catAsgn[cat].count(k)
	t2.append((k,c))
    temp = sorted(t2, key=itemgetter(1), reverse=True)
    for k in temp:
      asgn_file.write(str(k[0])+'\t'+str(k[1])+'\n')
    asgn_file.write('\n')  
  asgn_file.close()
  i = 0
  with open("data/wordindices.txt") as word_file:
    for w in word_file:
      w = w.strip()
      wordInd[i] = w
      i += 1
  
  with open("centers.txt") as cent_file:
    for line in cent_file:
      i = 0
      S = line.partition(":") # take the cluster id off
      temp = []
      probs = S[2].split()
      for p in probs:
	if float(p) > 0.0001:
	  word = wordInd[i]
	  temp.append((word,p))
	i += 1
	#sort by probability    
      clustWords[S[0]] = sorted(temp, key=itemgetter(1), reverse=True)
      
  fn = "clustWds.txt"
  wd_file = open(fn, "w")
  
  for clust in sorted(clustWords.keys()):
    wd_file.write("cluster: "+str(clust)+"\n")
    i = 0;
    for word in clustWords[clust]:
     # print word
      if i < 100:
	wd_file.write(str(word[0])+"\t"+ str(word[1])+"\n")
      i += 1
    wd_file.write("\n")
    
  wd_file.close()  




