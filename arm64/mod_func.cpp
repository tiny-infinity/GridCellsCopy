#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _itheta_reg(void);
extern void _itheta_stell_reg(void);
extern void _kdr_reg(void);
extern void _naf_reg(void);
extern void _stellate_mech_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mod/itheta.mod\"");
    fprintf(stderr, " \"mod/itheta_stell.mod\"");
    fprintf(stderr, " \"mod/kdr.mod\"");
    fprintf(stderr, " \"mod/naf.mod\"");
    fprintf(stderr, " \"mod/stellate_mech.mod\"");
    fprintf(stderr, "\n");
  }
  _itheta_reg();
  _itheta_stell_reg();
  _kdr_reg();
  _naf_reg();
  _stellate_mech_reg();
}

#if defined(__cplusplus)
}
#endif
