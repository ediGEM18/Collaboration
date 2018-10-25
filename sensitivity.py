from SALib.analyze import fast
import numpy as np
import params
import result_sorting
import pickler

def fast_sens(time, total_results):
    proteins_at_time = result_sorting.get_proteins_at_time(total_results, time)

    Y = np.array(proteins_at_time[0])
    problem = {'num_vars': 3,
               'names': ['plasmid_copies_BamA', 'mRNA_deg_rt_BamA', 'deg_rt_BamA'],
               'bounds': [[params.plasmid_copies_BamA[0], params.plasmid_copies_BamA[2]],
                          [params.mRNA_deg_rt_BamA[0], params.mRNA_deg_rt_BamA[3]],
                          [params.deg_rt_BamA[0], params.deg_rt_BamA[2]]]}

    print('FAST ' + str(time + 1) + ' minutes BamA')
    Si = fast.analyze(problem, Y, print_to_console=False)
    print('First order indices: ' + str(Si['S1']))
    print('Total order indices: ' + str(Si['ST']))

    Y = np.array(proteins_at_time[1])
    problem = {'num_vars': 3,
               'names': ['plasmid_copies_OmpA', 'mRNA_deg_rt_OmpA', 'deg_rt_OmpA'],
               'bounds': [[params.plasmid_copies_OmpA[0], params.plasmid_copies_OmpA[2]],
                          [params.mRNA_deg_rt_OmpA[0], params.mRNA_deg_rt_OmpA[3]],
                          [params.deg_rt_OmpA[0], params.deg_rt_OmpA[2]]]}

    print('FAST ' + str(time + 1) + ' minutes OmpA')
    Si = fast.analyze(problem, Y, print_to_console=False)
    print('First order indices: ' + str(Si['S1']))
    print('Total order indices: ' + str(Si['ST']))

    Y = np.array(proteins_at_time[2])
    problem = {'num_vars': 3,
               'names': ['plasmid_copies_lgA', 'mRNA_deg_rt_lgA', 'deg_rt_lgA'],
               'bounds': [[params.plasmid_copies_lgA[0], params.plasmid_copies_lgA[2]],
                          [params.mRNA_deg_rt_lgA[0], params.mRNA_deg_rt_lgA[3]],
                          [params.deg_rt_lgA[0], params.deg_rt_lgA[2]]]}

    print('FAST ' + str(time + 1) + ' minutes lgA')
    Si = fast.analyze(problem, Y, print_to_console=False)
    print('Total order indices: ' + str(Si['ST']))

print('reading...')
total_results = pickler.read('res_pickle')

for time in [(i* 20) - 1 for i in range(1, 13)]:
    fast_sens(time, total_results)
