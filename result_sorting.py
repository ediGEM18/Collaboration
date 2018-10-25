def get_mRNAs(total_results):
  mRNAs_BamA = []
  mRNAs_OmpA = []
  mRNAs_lgA = []
  for res in total_results:
    mRNAs_BamA.append(res[1][:, 0])
    mRNAs_OmpA.append(res[1][:, 2])
    mRNAs_lgA.append(res[1][:, 4])
  return (mRNAs_BamA, mRNAs_OmpA, mRNAs_lgA)

def get_mRNAs_at_time(total_results, time):
  mRNAs = get_mRNAs(total_results)
  mRNAs_at_time_BamA = [res[time] for res in mRNAs[0]]
  mRNAs_at_time_OmpA = [res[time] for res in mRNAs[1]]
  mRNAs_at_time_lgA = [res[time] for res in mRNAs[2]]
  return (mRNAs_at_time_BamA, mRNAs_at_time_OmpA, mRNAs_at_time_lgA)

def get_proteins(total_results):
  protein_BamA = []
  protein_OmpA = []
  protein_lgA = []
  for res in total_results:
    protein_BamA.append(res[1][:, 1])
    protein_OmpA.append(res[1][:, 3])
    protein_lgA.append(res[1][:, 5])
  return (protein_BamA, protein_OmpA, protein_lgA)

def get_proteins_at_time(total_results, time):
  proteins = get_proteins(total_results)
  protein_BamA = [res[time] for res in proteins[0]]
  protein_OmpA = [res[time] for res in proteins[1]]
  protein_lgA = [res[time] for res in proteins[2]]
  return (protein_BamA, protein_OmpA, protein_lgA)

def avg_mRNAs(total_results):
  mRNAs = get_mRNAs(total_results)
  avg_mRNA_BamA = [sum(i)/float(len(total_results)) for i in zip(*mRNAs[0])]
  avg_mRNA_OmpA = [sum(i)/float(len(total_results)) for i in zip(*mRNAs[1])]
  avg_mRNA_lgA = [sum(i)/float(len(total_results)) for i in zip(*mRNAs[2])]
  return (avg_mRNA_BamA, avg_mRNA_OmpA, avg_mRNA_lgA)

def min_mRNAs(total_results):
  mRNAs = get_mRNAs(total_results)
  min_mRNA_BamA = [min(i) for i in zip(*mRNAs[0])]
  min_mRNA_OmpA = [min(i) for i in zip(*mRNAs[1])]
  min_mRNA_lgA = [min(i) for i in zip(*mRNAs[2])]
  return (min_mRNA_BamA, min_mRNA_OmpA, min_mRNA_lgA)

def max_mRNAs(total_results):
  mRNAs = get_mRNAs(total_results)
  max_mRNA_BamA = [max(i) for i in zip(*mRNAs[0])]
  max_mRNA_OmpA = [max(i) for i in zip(*mRNAs[1])]
  max_mRNA_lgA = [max(i) for i in zip(*mRNAs[2])]
  return (max_mRNA_BamA, max_mRNA_OmpA, max_mRNA_lgA)

def avg_proteins(total_results):
  proteins = get_proteins(total_results)
  avg_BamA = [sum(i)/float(len(total_results)) for i in zip(*proteins[0])]
  avg_OmpA = [sum(i)/float(len(total_results)) for i in zip(*proteins[1])]
  avg_lgA = [sum(i)/float(len(total_results)) for i in zip(*proteins[2])]
  return (avg_BamA, avg_OmpA, avg_lgA)

def min_proteins(total_results):
  proteins = get_proteins(total_results)
  min_BamA = [min(i) for i in zip(*proteins[0])]
  min_OmpA = [min(i) for i in zip(*proteins[1])]
  min_lgA = [min(i) for i in zip(*proteins[2])]
  return (min_BamA, min_OmpA, min_lgA)

def max_proteins(total_results):
  proteins = get_proteins(total_results)
  max_BamA = [max(i) for i in zip(*proteins[0])]
  max_OmpA = [max(i) for i in zip(*proteins[1])]
  max_lgA = [max(i) for i in zip(*proteins[2])]
  return (max_BamA, max_OmpA, max_lgA)

def find_run(params, total_results):
  i = 0
  for res in total_results:
      if params == res[0]:
          return i
      else:
          i = i + 1

def find_runs(param1, loc1, param2, loc2, total_results):
  i = 0
  runs = []
  for res in total_results:
      if (param1 == res[0][loc1]) and (param2 == res[0][loc2]):
          runs.append(i)
      i = i + 1
  return runs

def get_run(run, total_results):

  mRNAs_BamA = total_results[run][1][:, 0]
  mRNAs_OmpA = total_results[run][1][:, 2]
  mRNAs_lgA = total_results[run][1][:, 4]

  protein_BamA = total_results[run][1][:, 1]
  protein_OmpA = total_results[run][1][:, 3]
  protein_lgA = total_results[run][1][:, 5]

  return (mRNAs_BamA, mRNAs_OmpA, mRNAs_lgA, protein_BamA, protein_OmpA, protein_lgA)
