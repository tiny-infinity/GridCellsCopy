/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__naf
#define _nrn_initial _nrn_initial__naf
#define nrn_cur _nrn_cur__naf
#define _nrn_current _nrn_current__naf
#define nrn_jacob _nrn_jacob__naf
#define nrn_state _nrn_state__naf
#define _net_receive _net_receive__naf 
#define _f_mh _f_mh__naf 
#define mh mh__naf 
#define states states__naf 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _p[0]
#define gmax_columnindex 0
#define i _p[1]
#define i_columnindex 1
#define g _p[2]
#define g_columnindex 2
#define h _p[3]
#define h_columnindex 3
#define ina _p[4]
#define ina_columnindex 4
#define Dh _p[5]
#define Dh_columnindex 5
#define _g _p[6]
#define _g_columnindex 6
#define _ion_ina	*_ppvar[0]._pval
#define _ion_dinadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpha(void);
 static void _hoc_beta(void);
 static void _hoc_mh(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_naf", _hoc_setdata,
 "alpha_naf", _hoc_alpha,
 "beta_naf", _hoc_beta,
 "mh_naf", _hoc_mh,
 0, 0
};
#define alpha alpha_naf
#define beta beta_naf
 extern double alpha( double , double );
 extern double beta( double , double );
 /* declare global and static user variables */
#define Inf Inf_naf
 double Inf[2];
#define Tau Tau_naf
 double Tau[2];
#define cai cai_naf
 double cai = 0;
#define cao cao_naf
 double cao = 0;
#define exptemp exptemp_naf
 double exptemp = 27;
#define erev erev_naf
 double erev = 55;
#define hexp hexp_naf
 double hexp = 1;
#define hq10 hq10_naf
 double hq10 = 5;
#define hbetaV0 hbetaV0_naf
 double hbetaV0 = -28;
#define hbetaB hbetaB_naf
 double hbetaB = -10;
#define hbetaA hbetaA_naf
 double hbetaA = 1;
#define hbflag hbflag_naf
 double hbflag = 2;
#define halphaV0 halphaV0_naf
 double halphaV0 = -58;
#define halphaB halphaB_naf
 double halphaB = -20;
#define halphaA halphaA_naf
 double halphaA = 0.07;
#define haflag haflag_naf
 double haflag = 1;
#define mexp mexp_naf
 double mexp = 3;
#define mq10 mq10_naf
 double mq10 = 5;
#define mbetaV0 mbetaV0_naf
 double mbetaV0 = -60;
#define mbetaB mbetaB_naf
 double mbetaB = -18;
#define mbetaA mbetaA_naf
 double mbetaA = 4;
#define mbflag mbflag_naf
 double mbflag = 1;
#define malphaV0 malphaV0_naf
 double malphaV0 = -35;
#define malphaB malphaB_naf
 double malphaB = -10;
#define malphaA malphaA_naf
 double malphaA = -0.1;
#define maflag maflag_naf
 double maflag = 3;
#define qq10 qq10_naf
 double qq10[2];
#define usetable usetable_naf
 double usetable = 1;
#define vmin vmin_naf
 double vmin = -100;
#define vmax vmax_naf
 double vmax = 100;
#define vrest vrest_naf
 double vrest = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_naf", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "erev_naf", "mV",
 "cao_naf", "mM",
 "cai_naf", "mM",
 "vmax_naf", "mV",
 "vmin_naf", "mV",
 "gmax_naf", "mho/cm2",
 "i_naf", "mA/cm^2",
 "g_naf", "mho/cm^2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "erev_naf", &erev_naf,
 "vrest_naf", &vrest_naf,
 "exptemp_naf", &exptemp_naf,
 "maflag_naf", &maflag_naf,
 "malphaA_naf", &malphaA_naf,
 "malphaB_naf", &malphaB_naf,
 "malphaV0_naf", &malphaV0_naf,
 "mbflag_naf", &mbflag_naf,
 "mbetaA_naf", &mbetaA_naf,
 "mbetaB_naf", &mbetaB_naf,
 "mbetaV0_naf", &mbetaV0_naf,
 "mq10_naf", &mq10_naf,
 "mexp_naf", &mexp_naf,
 "haflag_naf", &haflag_naf,
 "halphaA_naf", &halphaA_naf,
 "halphaB_naf", &halphaB_naf,
 "halphaV0_naf", &halphaV0_naf,
 "hbflag_naf", &hbflag_naf,
 "hbetaA_naf", &hbetaA_naf,
 "hbetaB_naf", &hbetaB_naf,
 "hbetaV0_naf", &hbetaV0_naf,
 "hq10_naf", &hq10_naf,
 "hexp_naf", &hexp_naf,
 "cao_naf", &cao_naf,
 "cai_naf", &cai_naf,
 "vmax_naf", &vmax_naf,
 "vmin_naf", &vmin_naf,
 "usetable_naf", &usetable_naf,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "Inf_naf", Inf_naf, 2,
 "Tau_naf", Tau_naf, 2,
 "qq10_naf", qq10_naf, 2,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"naf",
 "gmax_naf",
 0,
 "i_naf",
 "g_naf",
 0,
 "h_naf",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 7, _prop);
 	/*initialize range parameters*/
 	gmax = 0.035;
 	_prop->param = _p;
 	_prop->param_size = 7;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _naf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 7, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 naf /Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/naf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96489.0;
 static double R = 8.31441;
 static double *_t_Inf[2];
 static double *_t_Tau[2];
static int _reset;
static char *modelname = "Kevins Cvode modified Generalized Hodgkin-Huxley eqn Channel Model ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_mh(double);
static int mh(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_mh(double);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   mh ( _threadargscomma_ v ) ;
   Dh = ( - h + Inf [ 1 ] ) / Tau [ 1 ] ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 mh ( _threadargscomma_ v ) ;
 Dh = Dh  / (1. - dt*( ( ( - 1.0 ) ) / Tau[1] )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   mh ( _threadargscomma_ v ) ;
    h = h + (1. - exp(dt*(( ( - 1.0 ) ) / Tau[1])))*(- ( ( ( Inf[1] ) ) / Tau[1] ) / ( ( ( - 1.0 ) ) / Tau[1] ) - h) ;
   }
  return 0;
}
 static double _mfac_mh, _tmin_mh;
 static void _check_mh();
 static void _check_mh() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_maflag;
  static double _sav_malphaA;
  static double _sav_malphaB;
  static double _sav_malphaV0;
  static double _sav_mbflag;
  static double _sav_mbetaA;
  static double _sav_mbetaB;
  static double _sav_mbetaV0;
  static double _sav_exptemp;
  static double _sav_haflag;
  static double _sav_halphaA;
  static double _sav_halphaB;
  static double _sav_halphaV0;
  static double _sav_hbflag;
  static double _sav_hbetaA;
  static double _sav_hbetaB;
  static double _sav_hbetaV0;
  static double _sav_celsius;
  static double _sav_mq10;
  static double _sav_hq10;
  static double _sav_vrest;
  static double _sav_vmin;
  static double _sav_vmax;
  if (!usetable) {return;}
  if (_sav_maflag != maflag) { _maktable = 1;}
  if (_sav_malphaA != malphaA) { _maktable = 1;}
  if (_sav_malphaB != malphaB) { _maktable = 1;}
  if (_sav_malphaV0 != malphaV0) { _maktable = 1;}
  if (_sav_mbflag != mbflag) { _maktable = 1;}
  if (_sav_mbetaA != mbetaA) { _maktable = 1;}
  if (_sav_mbetaB != mbetaB) { _maktable = 1;}
  if (_sav_mbetaV0 != mbetaV0) { _maktable = 1;}
  if (_sav_exptemp != exptemp) { _maktable = 1;}
  if (_sav_haflag != haflag) { _maktable = 1;}
  if (_sav_halphaA != halphaA) { _maktable = 1;}
  if (_sav_halphaB != halphaB) { _maktable = 1;}
  if (_sav_halphaV0 != halphaV0) { _maktable = 1;}
  if (_sav_hbflag != hbflag) { _maktable = 1;}
  if (_sav_hbetaA != hbetaA) { _maktable = 1;}
  if (_sav_hbetaB != hbetaB) { _maktable = 1;}
  if (_sav_hbetaV0 != hbetaV0) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_mq10 != mq10) { _maktable = 1;}
  if (_sav_hq10 != hq10) { _maktable = 1;}
  if (_sav_vrest != vrest) { _maktable = 1;}
  if (_sav_vmin != vmin) { _maktable = 1;}
  if (_sav_vmax != vmax) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_mh =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_mh)/200.; _mfac_mh = 1./_dx;
   for (_i=0, _x=_tmin_mh; _i < 201; _x += _dx, _i++) {
    _f_mh(_x);
    for (_j = 0; _j < 2; _j++) { _t_Inf[_j][_i] = Inf[_j];
}    for (_j = 0; _j < 2; _j++) { _t_Tau[_j][_i] = Tau[_j];
}   }
   _sav_maflag = maflag;
   _sav_malphaA = malphaA;
   _sav_malphaB = malphaB;
   _sav_malphaV0 = malphaV0;
   _sav_mbflag = mbflag;
   _sav_mbetaA = mbetaA;
   _sav_mbetaB = mbetaB;
   _sav_mbetaV0 = mbetaV0;
   _sav_exptemp = exptemp;
   _sav_haflag = haflag;
   _sav_halphaA = halphaA;
   _sav_halphaB = halphaB;
   _sav_halphaV0 = halphaV0;
   _sav_hbflag = hbflag;
   _sav_hbetaA = hbetaA;
   _sav_hbetaB = hbetaB;
   _sav_hbetaV0 = hbetaV0;
   _sav_celsius = celsius;
   _sav_mq10 = mq10;
   _sav_hq10 = hq10;
   _sav_vrest = vrest;
   _sav_vmin = vmin;
   _sav_vmax = vmax;
  }
 }

 static int mh(double _lv){ _check_mh();
 _n_mh(_lv);
 return 0;
 }

 static void _n_mh(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_mh(_lv); return; 
}
 _xi = _mfac_mh * (_lv - _tmin_mh);
 if (isnan(_xi)) {
  for (_j = 0; _j < 2; _j++) { Inf[_j] = _xi;
}  for (_j = 0; _j < 2; _j++) { Tau[_j] = _xi;
}  return;
 }
 if (_xi <= 0.) {
 for (_j = 0; _j < 2; _j++) { Inf[_j] = _t_Inf[_j][0];
} for (_j = 0; _j < 2; _j++) { Tau[_j] = _t_Tau[_j][0];
} return; }
 if (_xi >= 200.) {
 for (_j = 0; _j < 2; _j++) { Inf[_j] = _t_Inf[_j][200];
} for (_j = 0; _j < 2; _j++) { Tau[_j] = _t_Tau[_j][200];
} return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 for (_j = 0; _j < 2; _j++) {double *_t = _t_Inf[_j]; Inf[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 for (_j = 0; _j < 2; _j++) {double *_t = _t_Tau[_j]; Tau[_j] = _t[_i] + _theta*(_t[_i+1] - _t[_i]);}
 }

 
static int  _f_mh (  double _lv ) {
   double _la , _lb , _lj ;
 qq10 [ 0 ] = pow( mq10 , ( ( celsius - exptemp ) / 10. ) ) ;
   qq10 [ 1 ] = pow( hq10 , ( ( celsius - exptemp ) / 10. ) ) ;
   {int  _lj ;for ( _lj = 0 ; _lj <= 1 ; _lj ++ ) {
     _la = alpha ( _threadargscomma_ _lv , ((double) _lj ) ) ;
     _lb = beta ( _threadargscomma_ _lv , ((double) _lj ) ) ;
     Inf [ _lj ] = _la / ( _la + _lb ) ;
     Tau [ _lj ] = 1. / ( _la + _lb ) / qq10 [ _lj ] ;
     if ( hexp  == 0.0 ) {
       Tau [ 1 ] = 1. ;
       Inf [ 1 ] = 1. ;
       }
     } }
    return 0; }
 
static void _hoc_mh(void) {
  double _r;
    _r = 1.;
 mh (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpha (  double _lv , double _lj ) {
   double _lalpha;
 double _lflag , _lA , _lB , _lV0 ;
 if ( _lj  == 1.0  && hexp  == 0.0 ) {
     _lalpha = 0.0 ;
     }
   else {
     if ( _lj  == 1.0 ) {
       _lA = halphaA ;
       _lB = halphaB ;
       _lV0 = halphaV0 + vrest ;
       _lflag = haflag ;
       }
     else {
       _lA = malphaA ;
       _lB = malphaB ;
       _lV0 = malphaV0 + vrest ;
       _lflag = maflag ;
       }
     if ( _lflag  == 1.0 ) {
       _lalpha = _lA * exp ( ( _lv - _lV0 ) / _lB ) ;
       }
     else if ( _lflag  == 2.0 ) {
       _lalpha = _lA / ( exp ( ( _lv - _lV0 ) / _lB ) + 1.0 ) ;
       }
     else if ( _lflag  == 3.0 ) {
       if ( _lv  == _lV0 ) {
         _lalpha = _lA * _lB ;
         }
       else {
         _lalpha = _lA * ( _lv - _lV0 ) / ( exp ( ( _lv - _lV0 ) / _lB ) - 1.0 ) ;
         }
       }
     }
   
return _lalpha;
 }
 
static void _hoc_alpha(void) {
  double _r;
   _r =  alpha (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double beta (  double _lv , double _lj ) {
   double _lbeta;
 double _lflag , _lA , _lB , _lV0 ;
 if ( _lj  == 1.0  && hexp  == 0.0 ) {
     _lbeta = 1.0 ;
     }
   else {
     if ( _lj  == 1.0 ) {
       _lA = hbetaA ;
       _lB = hbetaB ;
       _lV0 = hbetaV0 + vrest ;
       _lflag = hbflag ;
       }
     else {
       _lA = mbetaA ;
       _lB = mbetaB ;
       _lV0 = mbetaV0 + vrest ;
       _lflag = mbflag ;
       }
     if ( _lflag  == 1.0 ) {
       _lbeta = _lA * exp ( ( _lv - _lV0 ) / _lB ) ;
       }
     else if ( _lflag  == 2.0 ) {
       _lbeta = _lA / ( exp ( ( _lv - _lV0 ) / _lB ) + 1.0 ) ;
       }
     else if ( _lflag  == 3.0 ) {
       if ( _lv  == _lV0 ) {
         _lbeta = _lA * _lB ;
         }
       else {
         _lbeta = _lA * ( _lv - _lV0 ) / ( exp ( ( _lv - _lV0 ) / _lB ) - 1.0 ) ;
         }
       }
     }
   
return _lbeta;
 }
 
static void _hoc_beta(void) {
  double _r;
   _r =  beta (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
 {
   mh ( _threadargscomma_ v ) ;
   h = Inf [ 1 ] ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   mh ( _threadargscomma_ v ) ;
   g = gmax * Inf [ 0 ] * Inf [ 0 ] * Inf [ 0 ] * h ;
   i = g * ( v - erev ) ;
   ina = i ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  states();
 if(error){fprintf(stderr,"at line 100 in file naf.mod:\n  SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = h_columnindex;  _dlist1[0] = Dh_columnindex;
  for (_i=0; _i < 2; _i++) {  _t_Inf[_i] = makevector(201*sizeof(double)); }
  for (_i=0; _i < 2; _i++) {  _t_Tau[_i] = makevector(201*sizeof(double)); }
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/naf.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Wang, X.-J. and Buzsaki G. (1996) Gamma oscillations by synaptic\n"
  "inhibition in a hippocampal interneuronal network.  J. Neurosci. 16,\n"
  "6402-6413.\n"
  "ENDCOMMENT\n"
  ": $Id: naf.mod,v 1.7 2003/07/29 21:23:39 billl Exp $\n"
  "\n"
  "NEURON { SUFFIX naf }\n"
  "NEURON {  USEION na WRITE ina }\n"
  "ASSIGNED { ina }\n"
  "PARAMETER {\n"
  "	erev 		= 55.  (mV)\n"
  "	gmax 		= 0.035    (mho/cm2)\n"
  "        vrest           = 0.\n"
  "\n"
  "	exptemp		= 27\n"
  "	maflag 		= 3\n"
  "	malphaA 	= -0.1\n"
  "	malphaB		= -10.\n"
  "	malphaV0	= -35.\n"
  "	mbflag 		= 1\n"
  "	mbetaA 		= 4.\n"
  "	mbetaB		= -18.\n"
  "	mbetaV0		= -60.\n"
  "	mq10		= 5\n"
  "	mexp 		= 3\n"
  "\n"
  "	haflag 		= 1\n"
  "	halphaA 	= 0.07\n"
  "	halphaB		= -20\n"
  "	halphaV0	= -58.\n"
  "	hbflag 		= 2\n"
  "	hbetaA 		= 1.\n"
  "	hbetaB		= -10.\n"
  "	hbetaV0		= -28.\n"
  "	hq10		= 5\n"
  "	hexp 		= 1\n"
  "\n"
  "	cao                (mM)\n"
  "	cai                (mM)\n"
  "\n"
  "	celsius			   (degC)\n"
  "	dt 				   (ms)\n"
  "	v 			       (mV)\n"
  "\n"
  "	vmax 		= 100  (mV)\n"
  "	vmin 		= -100 (mV)\n"
  "} : end PARAMETER\n"
  "\n"
  ": $Id: naf.mod,v 1.7 2003/07/29 21:23:39 billl Exp $  \n"
  "TITLE Kevins Cvode modified Generalized Hodgkin-Huxley eqn Channel Model \n"
  "\n"
  "COMMENT\n"
  "\n"
  "Each channel has activation and inactivation particles as in the original\n"
  "Hodgkin Huxley formulation.  The activation particle mm and inactivation\n"
  "particle hh go from on to off states according to kinetic variables alpha\n"
  "and beta which are voltage dependent.\n"
  "Allows exponential, sigmoid and linoid forms (flags 0,1,2)\n"
  "See functions alpha() and beta() for details of parameterization\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	RANGE gmax, g, i\n"
  "	GLOBAL erev, Inf, Tau, vmin, vmax, vrest, qq10\n"
  "} : end NEURON\n"
  "\n"
  "CONSTANT {\n"
  "	  FARADAY = 96489.0	: Faraday's constant\n"
  "	  R= 8.31441		: Gas constant\n"
  "\n"
  "} : end CONSTANT\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "} : end UNITS\n"
  "\n"
  "ASSIGNED {\n"
  "	i (mA/cm^2)		\n"
  "	g (mho/cm^2)\n"
  "	Inf[2]		: 0 = m and 1 = h\n"
  "	Tau[2]		: 0 = m and 1 = h\n"
  "        qq10[2]\n"
  "} : end ASSIGNED \n"
  "\n"
  "STATE { h }\n"
  "\n"
  "INITIAL { \n"
  " 	mh(v)\n"
  "	h = Inf[1]\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "  SOLVE states METHOD cnexp\n"
  "  mh(v)\n"
  "  g = gmax * Inf[0]*Inf[0]*Inf[0] * h\n"
  "\n"
  "  i = g*(v-erev) \n"
  "  ina=i\n"
  "} : end BREAKPOINT\n"
  "\n"
  ": ASSIGNMENT PROCEDURES\n"
  ": Must be given by a user routines in parameters.multi\n"
  ": E.G.:\n"
  ":   PROCEDURE iassign () { i = g*(v-erev) ina=i }\n"
  ":   PROCEDURE iassign () { i = g*ghkca(v) ica=i }\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "\n"
  "DERIVATIVE states {\n"
  "	mh(v)\n"
  "	h' = (-h + Inf[1]) / Tau[1]\n"
  " }\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  ": NOTE : 0 = m and 1 = h\n"
  "PROCEDURE mh (v) {\n"
  "	LOCAL a, b, j\n"
  "	TABLE Inf, Tau DEPEND maflag, malphaA, malphaB, malphaV0, mbflag, mbetaA, mbetaB, mbetaV0, exptemp, haflag, halphaA, halphaB, halphaV0, hbflag, hbetaA, hbetaB, hbetaV0, celsius, mq10, hq10, vrest, vmin, vmax  FROM vmin TO vmax WITH 200\n"
  "\n"
  "	qq10[0] = mq10^((celsius-exptemp)/10.)	\n"
  "	qq10[1] = hq10^((celsius-exptemp)/10.)	\n"
  "\n"
  "	: Calculater Inf and Tau values for h and m\n"
  "	FROM j = 0 TO 1 {\n"
  "		a = alpha (v, j)\n"
  "		b = beta (v, j)\n"
  "\n"
  "		Inf[j] = a / (a + b)\n"
  "		Tau[j] = 1. / (a + b) / qq10[j]\n"
  "		if (hexp==0) { Tau[1] = 1. Inf[1] = 1.}\n"
  "	}\n"
  "} : end PROCEDURE mh (v)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "FUNCTION alpha(v,j) {\n"
  "  LOCAL flag, A, B, V0\n"
  "  if (j==1 && hexp==0) {\n"
  "	  alpha = 0\n"
  "  } else {\n"
  "\n"
  "     if (j == 1) {\n"
  "	  A = halphaA B = halphaB V0 = halphaV0+vrest flag = haflag\n"
  "     } else {\n"
  "	  A = malphaA B = malphaB V0 = malphaV0+vrest flag = maflag\n"
  "     }\n"
  "\n"
  "     if (flag == 1) { :  EXPONENTIAL\n"
  "	 alpha = A*exp((v-V0)/B)	\n"
  "     } else if (flag == 2) { :  SIGMOID\n"
  "	 alpha = A/(exp((v-V0)/B)+1)\n"
  "     } else if (flag == 3) { :  LINOID\n"
  "	 if(v == V0) {\n"
  "           alpha = A*B\n"
  "         } else {\n"
  "           alpha = A*(v-V0)/(exp((v-V0)/B)-1) }\n"
  "     }\n"
  "}\n"
  "} : end FUNCTION alpha (v,j)\n"
  "\n"
  ":-------------------------------------------------------------------\n"
  "FUNCTION beta (v,j) {\n"
  "  LOCAL flag, A, B, V0\n"
  "  if (j==1 && hexp==0) {\n"
  "	  beta = 1\n"
  "  } else {\n"
  "\n"
  "     if (j == 1) {\n"
  "	  A = hbetaA B = hbetaB V0 = hbetaV0+vrest flag = hbflag\n"
  "     } else {\n"
  "	  A = mbetaA B = mbetaB V0 = mbetaV0+vrest flag = mbflag\n"
  "     }\n"
  "\n"
  "    if (flag == 1) { :  EXPONENTIAL\n"
  "	 beta = A*exp((v-V0)/B)\n"
  "     } else if (flag == 2) { :  SIGMOID\n"
  "	 beta = A/(exp((v-V0)/B)+1)\n"
  "     } else if (flag == 3) { :  LINOID\n"
  "	 if(v == V0) {\n"
  "            beta = A*B \n"
  "         } else {\n"
  "            beta = A*(v-V0)/(exp((v-V0)/B)-1) }\n"
  "     }\n"
  "}\n"
  "} : end FUNCTION beta (v,j)\n"
  ;
#endif
