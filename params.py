import itertools
from math import log

#Variables
plasmid_copies_BamA = [2.9580509e15, 2.9580509e16, 1.1832204e17]
transcription_rt_BamA = [(60./2430.)*60.]
mRNA_deg_rt_BamA = [log(2)/1., log(2)/5., log(2)/10., log(2)/15.]
translation_rt_BamA = [(20./810.)*60.]
deg_rt_BamA = [log(2)/(10. * 60.), log(2)/(20. * 60.), log(2)/(30. * 60.)]

#vars OmpA
plasmid_copies_OmpA = [5.9017229e15, 5.9017229e16, 2.3606892e17]
transcription_rt_OmpA = [(60./1038.)*60.]
mRNA_deg_rt_OmpA = [log(2)/1., log(2)/5., log(2)/10., log(2)/15.]
translation_rt_OmpA = [(20./346.)*60.]
deg_rt_OmpA = [log(2)/(10. * 60.), log(2)/(20. * 60.), log(2)/(30. * 60.)]

#vars lgA
plasmid_copies_lgA = [6.0442428e15, 6.0442428e16, 2.4176971e17]
transcription_rt_lgA = [(60./945.)*60.]
mRNA_deg_rt_lgA = [log(2)/1., log(2)/5., log(2)/10., log(2)/15.]
translation_rt_lgA = [(20./315.)*60.]
deg_rt_lgA = [log(2)/(10. * 60.), log(2)/(20. * 60.), log(2)/(30. * 60.)]

#Initial Reagents
kDaltons_per_milligram = 6.022e20
sense_strand_BamA_mRNA_kDa = 830.382
anti_sense_strand_BamA_mRNA_kDa = 820.8
init_mRNA_BamA = [0.0, 
                  (kDaltons_per_milligram/((sense_strand_BamA_mRNA_kDa + anti_sense_strand_BamA_mRNA_kDa)/2)),
                  (2 * kDaltons_per_milligram/((sense_strand_BamA_mRNA_kDa + anti_sense_strand_BamA_mRNA_kDa)/2)),
                  (3 * kDaltons_per_milligram/((sense_strand_BamA_mRNA_kDa  + anti_sense_strand_BamA_mRNA_kDa)/2)),
                  (4 * kDaltons_per_milligram/((sense_strand_BamA_mRNA_kDa  + anti_sense_strand_BamA_mRNA_kDa)/2)),
                  (5 * kDaltons_per_milligram/((sense_strand_BamA_mRNA_kDa  + anti_sense_strand_BamA_mRNA_kDa)/2))]
init_BamA = [0.]

init_mRNA_OmpA = [0.]
init_OmpA = [0.]

init_mRNA_lgA = [0.]
init_lgA = [0.]

settings = list(itertools.product(plasmid_copies_BamA,
    transcription_rt_BamA,
    mRNA_deg_rt_BamA,
    translation_rt_BamA,
    deg_rt_BamA,
    plasmid_copies_OmpA,
    transcription_rt_OmpA,
    mRNA_deg_rt_OmpA,
    translation_rt_OmpA,
    deg_rt_OmpA,
    plasmid_copies_lgA,
    transcription_rt_lgA,
    mRNA_deg_rt_lgA,
    translation_rt_lgA,
    deg_rt_lgA,
    init_mRNA_BamA,
    init_BamA,
    init_mRNA_OmpA,
    init_OmpA,
    init_mRNA_lgA,
    init_lgA
))
