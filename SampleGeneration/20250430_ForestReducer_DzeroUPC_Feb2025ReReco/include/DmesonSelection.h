using namespace std;
#include <vector>

// ============================================================================
// CUT FUNCTIONS --------------------------------------------------------------

// General function to check any cut against a provided Dpt and Dy value
bool CheckCut(
  vector<vector<float>> DcutValue,
  float Dvalue,
  float Dpt,
  float Dy
) {
  bool pass = false;
  for (int iCut = 0; iCut < DcutValue.size(); iCut++) {
    float cutVal =  DcutValue[iCut][0];
    float cutSign = DcutValue[iCut][1];
    float Dptmin =  DcutValue[iCut][2];
    float Dptmax =  DcutValue[iCut][3];
    float Dymin =   DcutValue[iCut][4];
    float Dymax =   DcutValue[iCut][5];
    if (
      Dpt >= Dptmin && Dpt <= Dptmax &&
      Dy >= Dymin && Dy <= Dymax &&
      (Dvalue - cutVal) * cutSign >= 0.
    ) {
      pass = true;
      break;
    }
  }
  return pass;
}

// Flexible version of cut selection, accepts vector of cut selection matrixes
// that will be iterated through.
// New cut settings should use this!
bool DCutSelection(
  DzeroTreeMessenger &MDzero,
  int iD,
  vector<vector<vector<float>>> DCutValues
) {
  float pt = MDzero.Dpt[iD];
  float y = MDzero.Dy[iD];
  bool pass = true;
  for (int i = 0; i < DCutValues.size(); i++) {
    pass = pass * CheckCut(DCutValues[i], MDzero.Dchi2cl[iD], pt, y);
  }
  return pass;
}

// Old version of cut selection with explicit parameters.
bool DCutSelection(
  DzeroTreeMessenger &MDzero,
  int iD,
  vector<vector<float>> Dchi2clCutValue,
  vector<vector<float>> DalphaCutValue,
  vector<vector<float>> DdthetaCutValue,
  vector<vector<float>> DsvpvSigCutValue,
  vector<vector<float>> Dtrk1PtCutValue,
  vector<vector<float>> Dtrk2PtCutValue
) {
  float pt = MDzero.Dpt[iD];
  float y = MDzero.Dy[iD];
  bool pass = false;
  if (
    CheckCut(Dchi2clCutValue,   MDzero.Dchi2cl[iD],   pt, y) &&
    CheckCut(DalphaCutValue,    MDzero.Dalpha[iD],    pt, y) &&
    CheckCut(DdthetaCutValue,   MDzero.Ddtheta[iD],   pt, y) &&
    CheckCut(DsvpvSigCutValue,  MDzero.DsvpvDistance[iD] / MDzero.DsvpvDisErr[iD],  pt, y) &&
    CheckCut(Dtrk1PtCutValue,   MDzero.Dtrk1Pt[iD],   pt, y) &&
    CheckCut(Dtrk2PtCutValue,   MDzero.Dtrk2Pt[iD],   pt, y)
  ) {
    pass = true;
  }
  return pass;
}

// ============================================================================
// CUT SETTINGS ---------------------------------------------------------------
// Define cuts using:
//   {cut_value, pass_above/below, ptmin, ptmax, ymin, ymax}
//
// Notes:
// - pass_above/below indicates if values above the cut value are accepted (1)
//   or if values below the cut value are accepted (-1).
// - "_nom" is the nominal selection.
// - "_syst" are systematic settings. These should be looser than nominal".
// - "_loose" are looser-than-systematic settings, for cut validation
//   and/or BDT training

vector<vector<float>> Dchi2clCutValue_nom = {
  {0.1, 1,    0,  2,  -2.4, 2.4},
  {0.1, 1,    2, 12,  -2.4, 2.4}
};
vector<vector<float>> Dchi2clCutValue_syst = {
  {0.05, 1,   0,  2,  -2.4, 2.4},
  {0.05, 1,   2, 12,  -2.4, 2.4}
};
vector<vector<float>> Dchi2clCutValue_loose = {
  {0.03, 1,   0,  2,  -3, 3},
  {0.03, 1,   2, 12,  -3, 3},
};

vector<vector<float>> DalphaCutValue_nom = {
  {0.4,  -1,  0,  2,  -2.4, 2.4},
  {0.2,  -1,  2,  5,  -2.4,  -1},
  {0.4,  -1,  2,  5,    -1,   1},
  {0.2,  -1,  2,  5,     1, 2.4},
  {0.25, -1,  5,  8,  -2.4,  -1},
  {0.35, -1,  5,  8,    -1,   1},
  {0.25, -1,  5,  8,     1, 2.4},
  {0.25, -1,  8, 12,  -2.4,  -1},
  {0.4,  -1,  8, 12,    -1,   1},
  {0.25, -1,  8, 12,     1, 2.4}
};
vector<vector<float>> DalphaCutValue_syst = {
  {0.6,  -1,  0,  2,  -2.4, 2.4},
  {0.3,  -1,  2,  5,  -2.4,  -1},
  {0.6,  -1,  2,  5,    -1,   1},
  {0.3,  -1,  2,  5,     1, 2.4},
  {0.45, -1,  5,  8,  -2.4,  -1},
  {0.55, -1,  5,  8,    -1,   1},
  {0.45, -1,  5,  8,     1, 2.4},
  {0.45, -1,  8, 12,  -2.4,  -1},
  {0.6,  -1,  8, 12,    -1,   1},
  {0.45, -1,  8, 12,     1, 2.4}
};
vector<vector<float>> DalphaCutValue_loose = {
  {1.0,  -1,  0,  2,  -3,  3},
  {0.8,  -1,  2, 12,  -3,  3}
};

vector<vector<float>> DdthetaCutValue_nom = {
  {0.5,  -1,  0,  2,  -2.4, 2.4},
  {0.3,  -1,  2,  5,  -2.4,  -1},
  {0.5,  -1,  2,  5,    -1,   1},
  {0.3,  -1,  2,  5,     1, 2.4},
  {0.3,  -1,  5, 12,  -2.4, 2.4}
};
vector<vector<float>> DdthetaCutValue_syst = {
  {0.6,  -1,  0,  2,  -2.4, 2.4},
  {0.3,  -1,  2,  5,  -2.4,  -1},
  {0.6,  -1,  2,  5,    -1,   1},
  {0.3,  -1,  2,  5,     1, 2.4},
  {0.45, -1,  5,  8,  -2.4,  -1},
  {0.55, -1,  5,  8,    -1,   1},
  {0.45, -1,  5,  8,     1, 2.4},
  {0.45, -1,  8, 12,  -2.4,  -1},
  {0.6,  -1,  8, 12,    -1,   1},
  {0.45, -1,  8, 12,     1, 2.4}
};
vector<vector<float>> DdthetaCutValue_loose = {
  {1.0,  -1,  0,  2,  -3,  3},
  {0.8,  -1,  2, 12,  -3,  3}
};

vector<vector<float>> DsvpvSigCutValue_nom = {
  {2.5, 1,  0,  2, -2.4, 2.4},
  {2.5, 1,  2,  5, -2.4, 2.4},
  {3.5, 1,  5, 12, -2.4, 2.4}
};
vector<vector<float>> DsvpvSigCutValue_syst = {
  {2.0, 1,  0,  2, -2.4, 2.4},
  {2.0, 1,  2,  5, -2.4, 2.4},
  {2.5, 1,  5, 12, -2.4, 2.4}
};
vector<vector<float>> DsvpvSigCutValue_loose = {
  {1.5, 1,  0,  5, -3, 3},
  {2.0, 1,  5, 12, -3, 3},
};

vector<vector<float>> Dtrk1PtCutValue_nom = {
  {0.8, 1,  0,  2, -2.4, 2.4},
  {1.0, 1,  2, 12, -2.4, 2.4}
};
vector<vector<float>> Dtrk1PtCutValue_syst = {
  {0.5, 1,  0,  2, -2.4, 2.4},
  {0.7, 1,  2, 12, -2.4, 2.4}
};
vector<vector<float>> Dtrk1PtCutValue_loose = {
  {0.3, 1,  0,  2, -3, 3},
  {0.5, 1,  2, 12, -3, 3}
};
vector<vector<float>> Dtrk2PtCutValue_nom = Dtrk1PtCutValue_nom;
vector<vector<float>> Dtrk2PtCutValue_syst = Dtrk1PtCutValue_syst;
vector<vector<float>> Dtrk2PtCutValue_loose = Dtrk1PtCutValue_loose;

bool DpassCutNominal(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutNominal = {
    Dchi2clCutValue_nom,
    DalphaCutValue_nom,
    DdthetaCutValue_nom,
    DsvpvSigCutValue_nom,
    Dtrk1PtCutValue_nom,
    Dtrk2PtCutValue_nom
  };
  return DCutSelection(MDzero, iD, DCutNominal);;
}

bool DpassCutLoose(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutLoose = {
    Dchi2clCutValue_loose,
    DalphaCutValue_loose,
    DdthetaCutValue_loose,
    DsvpvSigCutValue_loose,
    Dtrk1PtCutValue_loose,
    Dtrk2PtCutValue_loose
  };
  return DCutSelection(MDzero, iD, DCutLoose);
}

bool DpassCutSystDsvpvSig(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutSystDsvpvSig = {
    Dchi2clCutValue_nom,
    DalphaCutValue_nom,
    DdthetaCutValue_nom,
    DsvpvSigCutValue_syst,
    Dtrk1PtCutValue_nom,
    Dtrk2PtCutValue_nom
  };
  return pass = DCutSelection(MDzero, iD, DCutSystDsvpvSig);
}

bool DpassCutSystDtrkPt(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutSystDtrkPt = {
    Dchi2clCutValue_nom,
    DalphaCutValue_nom,
    DdthetaCutValue_nom,
    DsvpvSigCutValue_nom,
    Dtrk1PtCutValue_syst,
    Dtrk2PtCutValue_syst
  };
  return DCutSelection(MDzero, iD, DCutSystDtrkPt);
}

bool DpassCutSystDalpha(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutSystDalpha = {
    Dchi2clCutValue_nom,
    DalphaCutValue_syst,
    DdthetaCutValue_syst,
    DsvpvSigCutValue_nom,
    Dtrk1PtCutValue_nom,
    Dtrk2PtCutValue_nom
  };
  return DCutSelection(MDzero, iD, DCutSystDalpha);
}

bool DpassCutSystDchi2cl(DzeroTreeMessenger &MDzero, int iD) {
  vector<vector<vector<float>>> DCutSystDchi2cl = {
    Dchi2clCutValue_syst,
    DalphaCutValue_nom,
    DdthetaCutValue_nom,
    DsvpvSigCutValue_nom,
    Dtrk1PtCutValue_nom,
    Dtrk2PtCutValue_nom
  };
  return DCutSelection(MDzero, iD, DCutSystDchi2cl);
}
