#include "classifier.h"

/* Autogenerated code from the script config2c.pl */

static const char *naccess_residue_name[] = {"A", "ALA", "ANY", "ARG", "ASN", "ASP", "C", "CYS", "DA", "DC", "DG", "DI", "DT", "DU", "G", "GLN", "GLU", "GLY", "HIS", "I", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SEC", "SER", "T", "THR", "TRP", "TYR", "U", "VAL", };
static const char *naccess_A_atom_name[] = {"C5", "C8", "N9", "N1", "N7", "N6", "C2", "N3", "C6", "C4", };
static double naccess_A_atom_radius[] = {1.80, 1.80, 1.60, 1.60, 1.60, 1.60, 1.80, 1.60, 1.80, 1.80, };
static int naccess_A_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_A_cfg = {
    10,
 "A",
    (char**) naccess_A_atom_name,
    (double*) naccess_A_atom_radius,
    (freesasa_atom_class*) naccess_A_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_ALA_atom_name[] = {"CB", };
static double naccess_ALA_atom_radius[] = {1.87, };
static int naccess_ALA_atom_class[] = {FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_ALA_cfg = {
    1,
 "ALA",
    (char**) naccess_ALA_atom_name,
    (double*) naccess_ALA_atom_radius,
    (freesasa_atom_class*) naccess_ALA_atom_class,
 {"ALA", 107.89, 43.94, 63.94, 36.71, 71.17, 0},
};

static const char *naccess_ANY_atom_name[] = {"OXT", "O5'", "C", "P", "N", "CB", "OP1", "C5'", "O", "CA", "C3'", "C4'", "C1'", "C2'", "O4'", "O3'", "OP3", "OP2", "O2'", };
static double naccess_ANY_atom_radius[] = {1.40, 1.40, 1.76, 1.90, 1.65, 1.87, 1.40, 1.80, 1.40, 1.87, 1.80, 1.80, 1.80, 1.80, 1.40, 1.40, 1.40, 1.40, 1.40, };
static int naccess_ANY_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_ANY_cfg = {
    19,
 "ANY",
    (char**) naccess_ANY_atom_name,
    (double*) naccess_ANY_atom_radius,
    (freesasa_atom_class*) naccess_ANY_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_ARG_atom_name[] = {"NH2", "CD", "NE", "CG", "CZ", "NH1", };
static double naccess_ARG_atom_radius[] = {1.65, 1.87, 1.65, 1.87, 1.76, 1.65, };
static int naccess_ARG_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_ARG_cfg = {
    6,
 "ARG",
    (char**) naccess_ARG_atom_name,
    (double*) naccess_ARG_atom_radius,
    (freesasa_atom_class*) naccess_ARG_atom_class,
 {"ARG", 238.33, 41.72, 196.61, 161.10, 77.23, 0},
};

static const char *naccess_ASN_atom_name[] = {"CG", "OD1", "ND2", };
static double naccess_ASN_atom_radius[] = {1.76, 1.40, 1.65, };
static int naccess_ASN_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_ASN_cfg = {
    3,
 "ASN",
    (char**) naccess_ASN_atom_name,
    (double*) naccess_ASN_atom_radius,
    (freesasa_atom_class*) naccess_ASN_atom_class,
 {"ASN", 143.97, 41.03, 102.94, 97.83, 46.14, 0},
};

static const char *naccess_ASP_atom_name[] = {"OD1", "OD2", "CG", };
static double naccess_ASP_atom_radius[] = {1.40, 1.40, 1.76, };
static int naccess_ASP_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_ASP_cfg = {
    3,
 "ASP",
    (char**) naccess_ASP_atom_name,
    (double*) naccess_ASP_atom_radius,
    (freesasa_atom_class*) naccess_ASP_atom_class,
 {"ASP", 140.48, 41.76, 98.72, 91.19, 49.29, 0},
};

static const char *naccess_C_atom_name[] = {"C6", "C4", "N3", "N4", "C2", "N1", "C5", "O2", };
static double naccess_C_atom_radius[] = {1.80, 1.80, 1.60, 1.60, 1.80, 1.60, 1.80, 1.40, };
static int naccess_C_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_C_cfg = {
    8,
 "C",
    (char**) naccess_C_atom_name,
    (double*) naccess_C_atom_radius,
    (freesasa_atom_class*) naccess_C_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_CYS_atom_name[] = {"SG", };
static double naccess_CYS_atom_radius[] = {1.85, };
static int naccess_CYS_atom_class[] = {FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_CYS_cfg = {
    1,
 "CYS",
    (char**) naccess_CYS_atom_name,
    (double*) naccess_CYS_atom_radius,
    (freesasa_atom_class*) naccess_CYS_atom_class,
 {"CYS", 134.24, 41.92, 92.33, 36.49, 97.75, 0},
};

static const char *naccess_DA_atom_name[] = {"N6", "N7", "N3", "C2", "C6", "C4", "C5", "C8", "N1", "N9", };
static double naccess_DA_atom_radius[] = {1.60, 1.60, 1.60, 1.80, 1.80, 1.80, 1.80, 1.80, 1.60, 1.60, };
static int naccess_DA_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_DA_cfg = {
    10,
 "DA",
    (char**) naccess_DA_atom_name,
    (double*) naccess_DA_atom_radius,
    (freesasa_atom_class*) naccess_DA_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_DC_atom_name[] = {"N1", "C5", "O2", "C6", "C4", "N3", "N4", "C2", };
static double naccess_DC_atom_radius[] = {1.60, 1.80, 1.40, 1.80, 1.80, 1.60, 1.60, 1.80, };
static int naccess_DC_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_DC_cfg = {
    8,
 "DC",
    (char**) naccess_DC_atom_name,
    (double*) naccess_DC_atom_radius,
    (freesasa_atom_class*) naccess_DC_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_DG_atom_name[] = {"N1", "N9", "C5", "C8", "O6", "C6", "C4", "N2", "N7", "N3", "C2", };
static double naccess_DG_atom_radius[] = {1.60, 1.60, 1.80, 1.80, 1.40, 1.80, 1.80, 1.60, 1.60, 1.60, 1.80, };
static int naccess_DG_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_DG_cfg = {
    11,
 "DG",
    (char**) naccess_DG_atom_name,
    (double*) naccess_DG_atom_radius,
    (freesasa_atom_class*) naccess_DG_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_DI_atom_name[] = {"N7", "N3", "C2", "O6", "C6", "C4", "C5", "C8", "N1", "N9", };
static double naccess_DI_atom_radius[] = {1.60, 1.60, 1.80, 1.40, 1.80, 1.80, 1.80, 1.80, 1.60, 1.60, };
static int naccess_DI_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_DI_cfg = {
    10,
 "DI",
    (char**) naccess_DI_atom_name,
    (double*) naccess_DI_atom_radius,
    (freesasa_atom_class*) naccess_DI_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_DT_atom_name[] = {"O4", "N1", "O2", "C5", "C4", "C7", "C6", "N3", "C2", };
static double naccess_DT_atom_radius[] = {1.40, 1.60, 1.40, 1.80, 1.80, 1.80, 1.80, 1.60, 1.80, };
static int naccess_DT_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_DT_cfg = {
    9,
 "DT",
    (char**) naccess_DT_atom_name,
    (double*) naccess_DT_atom_radius,
    (freesasa_atom_class*) naccess_DT_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_DU_atom_name[] = {"C2", "N3", "C6", "C4", "O2", "C5", "O4", "N1", };
static double naccess_DU_atom_radius[] = {1.80, 1.60, 1.80, 1.80, 1.40, 1.80, 1.40, 1.60, };
static int naccess_DU_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_DU_cfg = {
    8,
 "DU",
    (char**) naccess_DU_atom_name,
    (double*) naccess_DU_atom_radius,
    (freesasa_atom_class*) naccess_DU_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_G_atom_name[] = {"N1", "N9", "C8", "C5", "C4", "O6", "C6", "N3", "C2", "N2", "N7", };
static double naccess_G_atom_radius[] = {1.60, 1.60, 1.80, 1.80, 1.80, 1.40, 1.80, 1.60, 1.80, 1.60, 1.60, };
static int naccess_G_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_G_cfg = {
    11,
 "G",
    (char**) naccess_G_atom_name,
    (double*) naccess_G_atom_radius,
    (freesasa_atom_class*) naccess_G_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_GLN_atom_name[] = {"OE1", "CD", "CG", "NE2", };
static double naccess_GLN_atom_radius[] = {1.40, 1.76, 1.87, 1.65, };
static int naccess_GLN_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_GLN_cfg = {
    4,
 "GLN",
    (char**) naccess_GLN_atom_name,
    (double*) naccess_GLN_atom_radius,
    (freesasa_atom_class*) naccess_GLN_atom_class,
 {"GLN", 178.24, 41.72, 136.52, 126.35, 51.89, 0},
};

static const char *naccess_GLU_atom_name[] = {"OE2", "CG", "CD", "OE1", };
static double naccess_GLU_atom_radius[] = {1.40, 1.87, 1.76, 1.40, };
static int naccess_GLU_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_GLU_cfg = {
    4,
 "GLU",
    (char**) naccess_GLU_atom_name,
    (double*) naccess_GLU_atom_radius,
    (freesasa_atom_class*) naccess_GLU_atom_class,
 {"GLU", 172.09, 41.72, 130.37, 112.13, 59.96, 0},
};

static const char *naccess_GLY_atom_name[] = {"CA", };
static double naccess_GLY_atom_radius[] = {1.87, };
static int naccess_GLY_atom_class[] = {FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_GLY_cfg = {
    1,
 "GLY",
    (char**) naccess_GLY_atom_name,
    (double*) naccess_GLY_atom_radius,
    (freesasa_atom_class*) naccess_GLY_atom_class,
 {"GLY", 80.30, 80.30, 0.00, 42.62, 37.69, 0},
};

static const char *naccess_HIS_atom_name[] = {"CG", "CD2", "CE1", "NE2", "ND1", };
static double naccess_HIS_atom_radius[] = {1.76, 1.76, 1.76, 1.65, 1.65, };
static int naccess_HIS_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_HIS_cfg = {
    5,
 "HIS",
    (char**) naccess_HIS_atom_name,
    (double*) naccess_HIS_atom_radius,
    (freesasa_atom_class*) naccess_HIS_atom_class,
 {"HIS", 182.75, 38.76, 143.99, 85.61, 97.14, 0},
};

static const char *naccess_I_atom_name[] = {"C4", "O6", "C6", "N3", "C2", "N7", "N1", "N9", "C8", "C5", };
static double naccess_I_atom_radius[] = {1.80, 1.40, 1.80, 1.60, 1.80, 1.60, 1.60, 1.60, 1.80, 1.80, };
static int naccess_I_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_I_cfg = {
    10,
 "I",
    (char**) naccess_I_atom_name,
    (double*) naccess_I_atom_radius,
    (freesasa_atom_class*) naccess_I_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_ILE_atom_name[] = {"CG1", "CG2", "CD1", };
static double naccess_ILE_atom_radius[] = {1.87, 1.87, 1.87, };
static int naccess_ILE_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_ILE_cfg = {
    3,
 "ILE",
    (char**) naccess_ILE_atom_name,
    (double*) naccess_ILE_atom_radius,
    (freesasa_atom_class*) naccess_ILE_atom_class,
 {"ILE", 175.10, 41.16, 133.94, 36.10, 139.00, 0},
};

static const char *naccess_LEU_atom_name[] = {"CD2", "CG", "CD1", };
static double naccess_LEU_atom_radius[] = {1.87, 1.87, 1.87, };
static int naccess_LEU_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_LEU_cfg = {
    3,
 "LEU",
    (char**) naccess_LEU_atom_name,
    (double*) naccess_LEU_atom_radius,
    (freesasa_atom_class*) naccess_LEU_atom_class,
 {"LEU", 178.40, 39.50, 138.90, 36.45, 141.95, 0},
};

static const char *naccess_LYS_atom_name[] = {"NZ", "CG", "CE", "CD", };
static double naccess_LYS_atom_radius[] = {1.50, 1.87, 1.87, 1.87, };
static int naccess_LYS_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_LYS_cfg = {
    4,
 "LYS",
    (char**) naccess_LYS_atom_name,
    (double*) naccess_LYS_atom_radius,
    (freesasa_atom_class*) naccess_LYS_atom_class,
 {"LYS", 200.21, 41.72, 158.49, 84.31, 115.90, 0},
};

static const char *naccess_MET_atom_name[] = {"CE", "CG", "SD", };
static double naccess_MET_atom_radius[] = {1.87, 1.87, 1.85, };
static int naccess_MET_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_MET_cfg = {
    3,
 "MET",
    (char**) naccess_MET_atom_name,
    (double*) naccess_MET_atom_radius,
    (freesasa_atom_class*) naccess_MET_atom_class,
 {"MET", 193.72, 41.72, 152.00, 36.45, 157.27, 0},
};

static const char *naccess_PHE_atom_name[] = {"CD1", "CZ", "CG", "CD2", "CE2", "CE1", };
static double naccess_PHE_atom_radius[] = {1.76, 1.76, 1.76, 1.76, 1.76, 1.76, };
static int naccess_PHE_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_PHE_cfg = {
    6,
 "PHE",
    (char**) naccess_PHE_atom_name,
    (double*) naccess_PHE_atom_radius,
    (freesasa_atom_class*) naccess_PHE_atom_class,
 {"PHE", 199.40, 38.12, 161.27, 34.25, 165.14, 0},
};

static const char *naccess_PRO_atom_name[] = {"CD", "CG", };
static double naccess_PRO_atom_radius[] = {1.87, 1.87, };
static int naccess_PRO_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_PRO_cfg = {
    2,
 "PRO",
    (char**) naccess_PRO_atom_name,
    (double*) naccess_PRO_atom_radius,
    (freesasa_atom_class*) naccess_PRO_atom_class,
 {"PRO", 135.84, 27.09, 108.76, 15.17, 120.67, 0},
};

static const char *naccess_SEC_atom_name[] = {"SE", };
static double naccess_SEC_atom_radius[] = {1.80, };
static int naccess_SEC_atom_class[] = {FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_SEC_cfg = {
    1,
 "SEC",
    (char**) naccess_SEC_atom_name,
    (double*) naccess_SEC_atom_radius,
    (freesasa_atom_class*) naccess_SEC_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_SER_atom_name[] = {"OG", };
static double naccess_SER_atom_radius[] = {1.40, };
static int naccess_SER_atom_class[] = {FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_SER_cfg = {
    1,
 "SER",
    (char**) naccess_SER_atom_name,
    (double*) naccess_SER_atom_radius,
    (freesasa_atom_class*) naccess_SER_atom_class,
 {"SER", 116.56, 43.38, 73.18, 68.03, 48.53, 0},
};

static const char *naccess_T_atom_name[] = {"C5", "O2", "N1", "O4", "C2", "N3", "C6", "C7", "C4", };
static double naccess_T_atom_radius[] = {1.80, 1.40, 1.60, 1.40, 1.80, 1.60, 1.80, 1.80, 1.80, };
static int naccess_T_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_T_cfg = {
    9,
 "T",
    (char**) naccess_T_atom_name,
    (double*) naccess_T_atom_radius,
    (freesasa_atom_class*) naccess_T_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_THR_atom_name[] = {"OG1", "CG2", };
static double naccess_THR_atom_radius[] = {1.40, 1.87, };
static int naccess_THR_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_THR_cfg = {
    2,
 "THR",
    (char**) naccess_THR_atom_name,
    (double*) naccess_THR_atom_radius,
    (freesasa_atom_class*) naccess_THR_atom_class,
 {"THR", 139.25, 41.70, 97.55, 63.50, 75.75, 0},
};

static const char *naccess_TRP_atom_name[] = {"CZ3", "CZ2", "CE2", "CG", "CD1", "CE3", "CD2", "CH2", "NE1", };
static double naccess_TRP_atom_radius[] = {1.76, 1.76, 1.76, 1.76, 1.76, 1.76, 1.76, 1.76, 1.65, };
static int naccess_TRP_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_TRP_cfg = {
    9,
 "TRP",
    (char**) naccess_TRP_atom_name,
    (double*) naccess_TRP_atom_radius,
    (freesasa_atom_class*) naccess_TRP_atom_class,
 {"TRP", 248.97, 42.38, 206.59, 59.83, 189.14, 0},
};

static const char *naccess_TYR_atom_name[] = {"CE1", "CE2", "CG", "CD2", "OH", "CZ", "CD1", };
static double naccess_TYR_atom_radius[] = {1.76, 1.76, 1.76, 1.76, 1.40, 1.76, 1.76, };
static int naccess_TYR_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_TYR_cfg = {
    7,
 "TYR",
    (char**) naccess_TYR_atom_name,
    (double*) naccess_TYR_atom_radius,
    (freesasa_atom_class*) naccess_TYR_atom_class,
 {"TYR", 212.23, 38.13, 174.10, 76.34, 135.89, 0},
};

static const char *naccess_U_atom_name[] = {"N3", "C2", "C6", "C4", "O2", "C5", "N1", "O4", };
static double naccess_U_atom_radius[] = {1.60, 1.80, 1.80, 1.80, 1.40, 1.80, 1.60, 1.40, };
static int naccess_U_atom_class[] = {FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_APOLAR, FREESASA_ATOM_POLAR, FREESASA_ATOM_POLAR, };
static struct classifier_residue naccess_U_cfg = {
    8,
 "U",
    (char**) naccess_U_atom_name,
    (double*) naccess_U_atom_radius,
    (freesasa_atom_class*) naccess_U_atom_class,
 {NULL, 0, 0, 0, 0, 0},
};

static const char *naccess_VAL_atom_name[] = {"CG2", "CG1", };
static double naccess_VAL_atom_radius[] = {1.87, 1.87, };
static int naccess_VAL_atom_class[] = {FREESASA_ATOM_APOLAR, FREESASA_ATOM_APOLAR, };
static struct classifier_residue naccess_VAL_cfg = {
    2,
 "VAL",
    (char**) naccess_VAL_atom_name,
    (double*) naccess_VAL_atom_radius,
    (freesasa_atom_class*) naccess_VAL_atom_class,
 {"VAL", 151.40, 41.17, 110.23, 36.12, 115.28, 0},
};

static struct classifier_residue *naccess_residue_cfg[] = {
    &naccess_A_cfg, &naccess_ALA_cfg, &naccess_ANY_cfg, &naccess_ARG_cfg, &naccess_ASN_cfg, &naccess_ASP_cfg, &naccess_C_cfg, &naccess_CYS_cfg, &naccess_DA_cfg, &naccess_DC_cfg, &naccess_DG_cfg, &naccess_DI_cfg, &naccess_DT_cfg, &naccess_DU_cfg, &naccess_G_cfg, &naccess_GLN_cfg, &naccess_GLU_cfg, &naccess_GLY_cfg, &naccess_HIS_cfg, &naccess_I_cfg, &naccess_ILE_cfg, &naccess_LEU_cfg, &naccess_LYS_cfg, &naccess_MET_cfg, &naccess_PHE_cfg, &naccess_PRO_cfg, &naccess_SEC_cfg, &naccess_SER_cfg, &naccess_T_cfg, &naccess_THR_cfg, &naccess_TRP_cfg, &naccess_TYR_cfg, &naccess_U_cfg, &naccess_VAL_cfg, };

const freesasa_classifier freesasa_naccess_classifier = {
    34,    (char**) naccess_residue_name,
    "NACCESS",
    (struct classifier_residue **) naccess_residue_cfg,
};

