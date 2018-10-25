import matplotlib.pyplot as plt
plt.switch_backend('agg')

def end_ratios(average_proteins, min_proteins, max_proteins):
  avg_end_BamA = average_proteins[0][239]
  avg_end_OmpA = average_proteins[1][239]
  avg_end_lgA = average_proteins[2][239]

  min_end_BamA = min_proteins[0][239]
  min_end_OmpA = min_proteins[1][239]
  min_end_lgA = min_proteins[2][239]

  max_end_BamA = max_proteins[0][239]
  max_end_OmpA = max_proteins[1][239]
  max_end_lgA = max_proteins[2][239]

  print('Avg relative protein lvl after 2hrs: ' + str(avg_end_BamA) + ' : ' + str(avg_end_OmpA) + ' : ' + str(avg_end_lgA))
  print('Min relative protein lvl after 2hrs: ' + str(min_end_BamA) + ' : ' + str(min_end_OmpA) + ' : ' + str(min_end_lgA))
  print('Max relative protein lvl after 2hrs: ' + str(max_end_BamA) + ' : ' + str(max_end_OmpA) + ' : ' + str(max_end_lgA))

def plt_runs_of_interest(runs_of_interest, time):
  from itertools import product

  dna_labels = ['25ng DNA ', '250ng DNA ', '1000ng DNA ']
  rna_labels = ['0microg RNA', '1microg RNA', '2microg RNA', '3microg RNA', '4microg RNA', '5microg RNA',]

  fig = plt.figure(figsize=(30, 30))
  labels = [each[0] + each[1] for each in list(product(dna_labels, rna_labels))]
  for i in range(len(runs_of_interest)):
      plt.subplot(9,2,i+1)
      plt.plot(time, runs_of_interest[i], label=labels[i])

      plt.ylabel('Protein Present (molecules)')
      plt.xlabel('time (minutes)')
      plt.legend()

  plt.savefig('runs_of_interest.png')

def plt_avg_min_max(average_mRNAs, average_proteins, min_mRNAs, min_proteins, max_mRNAs, max_proteins, time):
  fig = plt.figure(figsize=(15, 15))

  #mRNA
  plt.subplot(231)
  plt.plot(time, average_mRNAs[0], label='Avg mRNA BamA')
  plt.plot(time, average_mRNAs[1], label='Avg mRNA OmpA')
  plt.plot(time, average_mRNAs[2], label='Avg mRNA ')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(232)
  plt.plot(time, min_mRNAs[0], label='Min mRNA BamA')
  plt.plot(time, min_mRNAs[1], label='Min mRNA OmpA')
  plt.plot(time, min_mRNAs[2], label='Min mRNA lgA')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(233)
  plt.plot(time, max_mRNAs[0], label='Max mRNA BamA')
  plt.plot(time, max_mRNAs[1], label='Max mRNA OmpA')
  plt.plot(time, max_mRNAs[2], label='Max mRNA lgA')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  #protein
  plt.subplot(234)
  plt.plot(time, average_proteins[0], label='Avg BamA')
  plt.plot(time, average_proteins[1], label='Avg OmpA')
  plt.plot(time, average_proteins[2], label='Avg lgA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(235)
  plt.plot(time, min_proteins[0], label='Min BamA')
  plt.plot(time, min_proteins[1], label='Min OmpA')
  plt.plot(time, min_proteins[2], label='Min lgA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(236)
  plt.plot(time, max_proteins[0], label='Max BamA')
  plt.plot(time, max_proteins[1], label='Max OmpA')
  plt.plot(time, max_proteins[2], label='Max lgA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.savefig('model.png')

def plt_run(run_result, time):
  fig = plt.figure(figsize=(15, 15))

  #mRNA
  plt.subplot(2,3,1)
  plt.plot(time, run_result[0], label='mRNA BamA')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(2,3,2)
  plt.plot(time, run_result[2], label='mRNA OmpA')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(2,3,3)
  plt.plot(time, run_result[4], label='mRNA lgA')

  plt.ylabel('mRNA Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  #protein
  plt.subplot(2,3,4)
  plt.plot(time, run_result[1], label='BamA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(2,3,5)
  plt.plot(time, run_result[3], label='OmpA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.subplot(2,3,6)
  plt.plot(time, run_result[5], label='lgA')

  plt.ylabel('Protein Present (molecules)')
  plt.xlabel('time (minutes)')
  plt.legend()

  plt.savefig('run.png')

def plt_res_scatter(mRNAs, proteins, time):
  unnest_mRNAs_BamA = sum(mRNAs[0], [])
  unnest_protein_BamA = sum(proteins[0], [])

  unnest_mRNAs_OmpA = sum(mRNAs[1], [])
  unnest_protein_OmpA = sum(proteins[1], [])

  unnnest_mRNAs_lgA = sum(mRNAs[2], [])
  unnnest_mRNAs_lgA = sum(proteins[2], [])

  rep_time = time * 279936

  plt.subplot(231)
  plt.scatter(unnest_mRNAs_BamA, rep_time)

  plt.subplot(232)
  plt.scatter(unnest_mRNAs_OmpA, rep_time)

  plt.subplot(233)
  plt.scatter(unnest_mRNAs_lgA, rep_time)

  plt.subplot(234)
  plt.scatter(unnest_protein_BamA, rep_time)

  plt.subplot(235)
  plt.scatter(unnest_protein_OmpA, rep_time)

  plt.subplot(236)
  plt.scatter(unnest_protein_lgA, rep_time)

  plt.savefig('scatters.png')

def plt_BamA_mRNA_dists(mRNAs_per_time):
  fig = plt.figure(figsize=(15, 15))

  axis_end = max(mRNAs_per_time[11][0])

  plt.subplot(621)
  plt.hist(mRNAs_per_time[0][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(622)
  plt.hist(mRNAs_per_time[1][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(623)
  plt.hist(mRNAs_per_time[2][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(624)
  plt.hist(mRNAs_per_time[3][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(625)
  plt.hist(mRNAs_per_time[4][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(626)
  plt.hist(mRNAs_per_time[5][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(627)
  plt.hist(mRNAs_per_time[6][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(628)
  plt.hist(mRNAs_per_time[7][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(629)
  plt.hist(mRNAs_per_time[8][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(mRNAs_per_time[9][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(mRNAs_per_time[10][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(mRNAs_per_time[11][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('BamA_mRNA_dists.png')

def plt_OmpA_mRNA_dists(mRNAs_per_time):
  fig = plt.figure(figsize=(15, 15))

  axis_end = max(mRNAs_per_time[11][1])

  plt.subplot(621)
  plt.hist(mRNAs_per_time[0][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(622)
  plt.hist(mRNAs_per_time[1][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(623)
  plt.hist(mRNAs_per_time[2][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(624)
  plt.hist(mRNAs_per_time[3][1],bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(625)
  plt.hist(mRNAs_per_time[4][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(626)
  plt.hist(mRNAs_per_time[5][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(627)
  plt.hist(mRNAs_per_time[6][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(628)
  plt.hist(mRNAs_per_time[7][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(629)
  plt.hist(mRNAs_per_time[8][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(mRNAs_per_time[9][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(mRNAs_per_time[10][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(mRNAs_per_time[11][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('OmpA_mRNA_dists.png')

def plt_lgA_mRNA_dists(mRNAs_per_time):
  fig = plt.figure(figsize=(15, 15))

  axis_end = max(mRNAs_per_time[11][2])

  plt.subplot(621)
  plt.hist(mRNAs_per_time[0][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(622)
  plt.hist(mRNAs_per_time[1][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(623)
  plt.hist(mRNAs_per_time[2][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(624)
  plt.hist(mRNAs_per_time[3][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(625)
  plt.hist(mRNAs_per_time[4][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(626)
  plt.hist(mRNAs_per_time[5][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(627)
  plt.hist(mRNAs_per_time[6][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(628)
  plt.hist(mRNAs_per_time[7][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(629)
  plt.hist(mRNAs_per_time[8][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(mRNAs_per_time[9][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(mRNAs_per_time[10][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(mRNAs_per_time[11][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('mRNA Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('lgA_mRNA_dists.png')

def plt_BamA_protein_dists(proteins_per_time):
  fig = plt.figure(figsize=(20, 20))

  axis_end = max(proteins_per_time[11][0])

  plt.subplot(621)
  plt.hist(proteins_per_time[0][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(622)
  plt.hist(proteins_per_time[1][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(623)
  plt.hist(proteins_per_time[2][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(624)
  plt.hist(proteins_per_time[3][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(625)
  plt.hist(proteins_per_time[4][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(626)
  plt.hist(proteins_per_time[5][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(627)
  plt.hist(proteins_per_time[6][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(628)
  plt.hist(proteins_per_time[7][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(629)
  plt.hist(proteins_per_time[8][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(proteins_per_time[9][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(proteins_per_time[10][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(proteins_per_time[11][0], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Proteins Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('BamA_dists.png')

def plt_OmpA_protein_dists(proteins_per_time):
  fig = plt.figure(figsize=(20, 20))

  axis_end = max(proteins_per_time[11][1])

  plt.subplot(621)
  plt.hist(proteins_per_time[0][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(622)
  plt.hist(proteins_per_time[1][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(623)
  plt.hist(proteins_per_time[2][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(624)
  plt.hist(proteins_per_time[3][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(625)
  plt.hist(proteins_per_time[4][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(626)
  plt.hist(proteins_per_time[5][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(627)
  plt.hist(proteins_per_time[6][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(628)
  plt.hist(proteins_per_time[7][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(629)
  plt.hist(proteins_per_time[8][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(proteins_per_time[9][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(proteins_per_time[10][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(proteins_per_time[11][1], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('OmpA_dists.png')

def plt_lgA_protein_dists(proteins_per_time):
  fig = plt.figure(figsize=(20, 20))

  axis_end = max(proteins_per_time[11][2])

  plt.subplot(621)
  plt.hist(proteins_per_time[0][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(622)
  plt.hist(proteins_per_time[1][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(623)
  plt.hist(proteins_per_time[2][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(624)
  plt.hist(proteins_per_time[3][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(625)
  plt.hist(proteins_per_time[4][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(626)
  plt.hist(proteins_per_time[5][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(627)
  plt.hist(proteins_per_time[6][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(628)
  plt.hist(proteins_per_time[7][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(629)
  plt.hist(proteins_per_time[8][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,10)
  plt.hist(proteins_per_time[9][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,11)
  plt.hist(proteins_per_time[10][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  plt.subplot(6,2,12)
  plt.hist(proteins_per_time[11][2], bins=200, range=[0.0, axis_end])

  plt.ylabel('count')
  plt.xlabel('Protein Present (molecules)')

  fig.subplots_adjust(hspace=1)

  plt.savefig('lgA_dists.png')
