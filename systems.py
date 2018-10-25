class Systems():

  def __init__(self, settings):
    #vars BamA
    self.plasmid_copies_BamA = settings[0]
    self.transcription_rt_BamA = settings[1]
    self.mRNA_deg_rt_BamA = settings[2]
    self.translation_rt_BamA = settings[3]
    self.deg_rt_BamA = settings[4]

    #vars OmpA
    self.plasmid_copies_OmpA = settings[5]
    self.transcription_rt_OmpA = settings[6]
    self.mRNA_deg_rt_OmpA = settings[7]
    self.translation_rt_OmpA = settings[8]
    self.deg_rt_OmpA = settings[9]

    #vars lgA
    self.plasmid_copies_lgA = settings[10]
    self.transcription_rt_lgA = settings[11]
    self.mRNA_deg_rt_lgA = settings[12]
    self.translation_rt_lgA = settings[13]
    self.deg_rt_lgA = settings[14]

  def sys(self, y, t):
      #reactants BamA
      mRNA_BamA = y[0]
      protein_BamA = y[1]

      #reactants OmpA
      mRNA_OmpA = y[2]
      protein_OmpA = y[3]

      #reactants lgA
      mRNA_lgA = y[4]
      protein_lgA = y[5]

      #equations BamA
      d_BamA_mRNA_dt = (self.plasmid_copies_BamA * self.transcription_rt_BamA) - (self.mRNA_deg_rt_BamA * mRNA_BamA)
      d_BamA_dt = (self.translation_rt_BamA * mRNA_BamA) - (self.deg_rt_BamA * protein_BamA)

      #equations OmpA
      d_OmpA_mRNA_dt = (self.plasmid_copies_OmpA * self.transcription_rt_OmpA) - (self.mRNA_deg_rt_OmpA * mRNA_OmpA)
      d_OmpA_dt = (self.translation_rt_OmpA * mRNA_OmpA) - (self.deg_rt_OmpA * protein_OmpA)

      #equations lgA
      d_lgA_mRNA_dt = (self.plasmid_copies_lgA * self.transcription_rt_lgA) - (self.mRNA_deg_rt_lgA * mRNA_lgA)
      d_lgA_dt = (self.translation_rt_lgA * mRNA_lgA) - (self.deg_rt_lgA * protein_lgA)
      
      return [d_BamA_mRNA_dt, d_BamA_dt, d_OmpA_mRNA_dt, d_OmpA_dt, d_lgA_mRNA_dt, d_lgA_dt]
