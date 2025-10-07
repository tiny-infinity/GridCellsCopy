/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__stellate_mech
#define _nrn_initial _nrn_initial__stellate_mech
#define nrn_cur _nrn_cur__stellate_mech
#define _nrn_current _nrn_current__stellate_mech
#define nrn_jacob _nrn_jacob__stellate_mech
#define nrn_state _nrn_state__stellate_mech
#define _net_receive _net_receive__stellate_mech 
#define rates rates__stellate_mech 
#define states states__stellate_mech 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnap_bar _p[0]
#define gnap_bar_columnindex 0
#define gnabar _p[1]
#define gnabar_columnindex 1
#define gkbar _p[2]
#define gkbar_columnindex 2
#define ghbar _p[3]
#define ghbar_columnindex 3
#define gna _p[4]
#define gna_columnindex 4
#define gh _p[5]
#define gh_columnindex 5
#define il _p[6]
#define il_columnindex 6
#define ih _p[7]
#define ih_columnindex 7
#define mna _p[8]
#define mna_columnindex 8
#define hna _p[9]
#define hna_columnindex 9
#define mnap _p[10]
#define mnap_columnindex 10
#define nk _p[11]
#define nk_columnindex 11
#define mhf _p[12]
#define mhf_columnindex 12
#define mhs _p[13]
#define mhs_columnindex 13
#define ena _p[14]
#define ena_columnindex 14
#define ek _p[15]
#define ek_columnindex 15
#define gk _p[16]
#define gk_columnindex 16
#define gnap _p[17]
#define gnap_columnindex 17
#define ina _p[18]
#define ina_columnindex 18
#define ik _p[19]
#define ik_columnindex 19
#define amna _p[20]
#define amna_columnindex 20
#define bmna _p[21]
#define bmna_columnindex 21
#define ahna _p[22]
#define ahna_columnindex 22
#define bhna _p[23]
#define bhna_columnindex 23
#define ank _p[24]
#define ank_columnindex 24
#define bnk _p[25]
#define bnk_columnindex 25
#define mnap_beta _p[26]
#define mnap_beta_columnindex 26
#define mnap_alpha _p[27]
#define mnap_alpha_columnindex 27
#define mhfinf _p[28]
#define mhfinf_columnindex 28
#define mhsinf _p[29]
#define mhsinf_columnindex 29
#define mhstau _p[30]
#define mhstau_columnindex 30
#define mhftau _p[31]
#define mhftau_columnindex 31
#define Dmna _p[32]
#define Dmna_columnindex 32
#define Dhna _p[33]
#define Dhna_columnindex 33
#define Dmnap _p[34]
#define Dmnap_columnindex 34
#define Dnk _p[35]
#define Dnk_columnindex 35
#define Dmhf _p[36]
#define Dmhf_columnindex 36
#define Dmhs _p[37]
#define Dmhs_columnindex 37
#define v _p[38]
#define v_columnindex 38
#define _g _p[39]
#define _g_columnindex 39
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_vtrap(void);
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_stellate_mech", _hoc_setdata,
 "rates_stellate_mech", _hoc_rates,
 "vtrap_stellate_mech", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_stellate_mech
 extern double vtrap( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define eh eh_stellate_mech
 double eh = -20;
#define el el_stellate_mech
 double el = -65;
#define gl gl_stellate_mech
 double gl = 0.0005;
#define hs_tau_input hs_tau_input_stellate_mech
 double hs_tau_input = 5.6;
#define hf_tau_input hf_tau_input_stellate_mech
 double hf_tau_input = 0.51;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "el_stellate_mech", "mV",
 "gl_stellate_mech", "S/cm2",
 "eh_stellate_mech", "mV",
 "hf_tau_input_stellate_mech", "ms",
 "hs_tau_input_stellate_mech", "ms",
 "gnap_bar_stellate_mech", "S/cm2",
 "gnabar_stellate_mech", "S/cm2",
 "gkbar_stellate_mech", "S/cm2",
 "ghbar_stellate_mech", "S/cm2",
 "gna_stellate_mech", "S/cm2",
 "gh_stellate_mech", "S/cm2",
 "il_stellate_mech", "mA/cm2",
 "ih_stellate_mech", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double hna0 = 0;
 static double mhs0 = 0;
 static double mhf0 = 0;
 static double mnap0 = 0;
 static double mna0 = 0;
 static double nk0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "el_stellate_mech", &el_stellate_mech,
 "gl_stellate_mech", &gl_stellate_mech,
 "eh_stellate_mech", &eh_stellate_mech,
 "hf_tau_input_stellate_mech", &hf_tau_input_stellate_mech,
 "hs_tau_input_stellate_mech", &hs_tau_input_stellate_mech,
 0,0
};
 static DoubVec hoc_vdoub[] = {
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
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"stellate_mech",
 "gnap_bar_stellate_mech",
 "gnabar_stellate_mech",
 "gkbar_stellate_mech",
 "ghbar_stellate_mech",
 0,
 "gna_stellate_mech",
 "gh_stellate_mech",
 "il_stellate_mech",
 "ih_stellate_mech",
 0,
 "mna_stellate_mech",
 "hna_stellate_mech",
 "mnap_stellate_mech",
 "nk_stellate_mech",
 "mhf_stellate_mech",
 "mhs_stellate_mech",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 40, _prop);
 	/*initialize range parameters*/
 	gnap_bar = 0.0005;
 	gnabar = 0.052;
 	gkbar = 0.011;
 	ghbar = 0.0015;
 	_prop->param = _p;
 	_prop->param_size = 40;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _stellate_mech_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 40, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 stellate_mech /Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/stellate_mech.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Stellate cells mechanism :acker";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[6], _dlist1[6];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dmna = amna * ( 1.0 - mna ) - bmna * mna ;
   Dhna = ahna * ( 1.0 - hna ) - bhna * hna ;
   Dmnap = mnap_alpha * ( 1.0 - mnap ) - mnap_beta * mnap ;
   Dnk = ank * ( 1.0 - nk ) - bnk * nk ;
   Dmhf = ( mhfinf - mhf ) / ( mhftau ) ;
   Dmhs = ( mhsinf - mhs ) / ( mhstau ) ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dmna = Dmna  / (1. - dt*( ( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ) )) ;
 Dhna = Dhna  / (1. - dt*( ( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ) )) ;
 Dmnap = Dmnap  / (1. - dt*( ( mnap_alpha )*( ( ( - 1.0 ) ) ) - ( mnap_beta )*( 1.0 ) )) ;
 Dnk = Dnk  / (1. - dt*( ( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ) )) ;
 Dmhf = Dmhf  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mhftau ) )) ;
 Dmhs = Dmhs  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( mhstau ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    mna = mna + (1. - exp(dt*(( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ))))*(- ( ( amna )*( ( 1.0 ) ) ) / ( ( amna )*( ( ( - 1.0 ) ) ) - ( bmna )*( 1.0 ) ) - mna) ;
    hna = hna + (1. - exp(dt*(( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ))))*(- ( ( ahna )*( ( 1.0 ) ) ) / ( ( ahna )*( ( ( - 1.0 ) ) ) - ( bhna )*( 1.0 ) ) - hna) ;
    mnap = mnap + (1. - exp(dt*(( mnap_alpha )*( ( ( - 1.0 ) ) ) - ( mnap_beta )*( 1.0 ))))*(- ( ( mnap_alpha )*( ( 1.0 ) ) ) / ( ( mnap_alpha )*( ( ( - 1.0 ) ) ) - ( mnap_beta )*( 1.0 ) ) - mnap) ;
    nk = nk + (1. - exp(dt*(( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ))))*(- ( ( ank )*( ( 1.0 ) ) ) / ( ( ank )*( ( ( - 1.0 ) ) ) - ( bnk )*( 1.0 ) ) - nk) ;
    mhf = mhf + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mhftau ))))*(- ( ( ( mhfinf ) ) / ( mhftau ) ) / ( ( ( ( - 1.0 ) ) ) / ( mhftau ) ) - mhf) ;
    mhs = mhs + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ( mhstau ))))*(- ( ( ( mhsinf ) ) / ( mhstau ) ) / ( ( ( ( - 1.0 ) ) ) / ( mhstau ) ) - mhs) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   amna = ( .1 ) * vtrap ( _threadargscomma_ - ( _lv + 23.0 ) , 10.0 ) ;
   bmna = 4.0 * exp ( - ( _lv + 48.0 ) / 18.0 ) ;
   ahna = 0.07 * exp ( - ( _lv + 37.0 ) / 20.0 ) ;
   bhna = 1.0 / ( exp ( - 0.1 * ( _lv + 7.0 ) ) + 1.0 ) ;
   mnap_alpha = 1.0 / ( 0.15 * ( exp ( - ( _lv + 38.0 ) / 6.5 ) + 1.0 ) ) ;
   mnap_beta = ( exp ( - ( _lv + 38.0 ) / 6.5 ) ) / ( 0.15 * ( exp ( - ( _lv + 38.0 ) / 6.5 ) + 1.0 ) ) ;
   ank = 0.01 * vtrap ( _threadargscomma_ - ( _lv + 27.0 ) , 10.0 ) ;
   bnk = 0.125 * exp ( - ( _lv + 37.0 ) / 80.0 ) ;
   mhfinf = 1.0 / ( 1.0 + exp ( ( _lv + 79.2 ) / 9.78 ) ) ;
   mhftau = ( ( hf_tau_input / ( ( exp ( ( _lv - 1.7 ) / 10.0 ) ) + exp ( - ( _lv + 340.0 ) / 52.0 ) ) ) + 1.0 ) ;
   mhsinf = 1.0 / ( 1.0 + exp ( ( _lv + 71.3 ) / 7.9 ) ) ;
   mhstau = ( ( hs_tau_input / ( ( exp ( ( _lv - 1.7 ) / 14.0 ) ) + exp ( - ( _lv + 260.0 ) / 43.0 ) ) ) + 1.0 ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 6;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 6; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  hna = hna0;
  mnap = mnap0;
  mna = mna0;
  mhs = mhs0;
  mhf = mhf0;
  nk = nk0;
 {
   rates ( _threadargscomma_ v ) ;
   mna = amna / ( amna + bmna ) ;
   hna = ahna / ( ahna + bhna ) ;
   mnap = mnap_alpha / ( mnap_alpha + mnap_beta ) ;
   nk = ank / ( ank + bnk ) ;
   mhf = mhfinf ;
   mhs = mhsinf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
  }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = ( gnabar * mna * mna * mna * hna ) ;
   gnap = ( gnap_bar * mnap ) ;
   ina = ( gna + gnap ) * ( v - ena ) ;
   gk = ( gkbar * nk * nk * nk * nk ) ;
   ik = gk * ( v - ek ) ;
   il = gl * ( v - el ) ;
   gh = ghbar * ( 0.65 * mhf + 0.35 * mhs ) ;
   ih = gh * ( v - eh ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;
 _current += ih;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  ena = _ion_ena;
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  }  }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = mna_columnindex;  _dlist1[0] = Dmna_columnindex;
 _slist1[1] = hna_columnindex;  _dlist1[1] = Dhna_columnindex;
 _slist1[2] = mnap_columnindex;  _dlist1[2] = Dmnap_columnindex;
 _slist1[3] = nk_columnindex;  _dlist1[3] = Dnk_columnindex;
 _slist1[4] = mhf_columnindex;  _dlist1[4] = Dmhf_columnindex;
 _slist1[5] = mhs_columnindex;  _dlist1[5] = Dmhs_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/stellate_mech.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Acker, Corey D., Nancy Kopell, and John A. White. \n"
  "\n"
  "\n"
  "Synchronization of Strongly \n"
  "Coupled Excitatory Neurons: Relating Network Behavior to Biophysics.\n"
  "\n"
  "\n"
  " \n"
  "Journal of Computational Neuroscience 15, no. 1 (July 1, 2003): 71\n"
  "\n"
  "\n"
  "90. \n"
  "https://doi.org/10.1023/A:1024474819512.\n"
  "ENDCOMMENT\n"
  "\n"
  "TITLE Stellate cells mechanism :acker\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (mV)=(millivolt)\n"
  "    (S) = (siemens)\n"
  "    (mA) = (milliamp)\n"
  "    \n"
  "}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX stellate_mech\n"
  "    USEION na READ ena WRITE ina \n"
  "    USEION k READ ek WRITE ik \n"
  "    NONSPECIFIC_CURRENT il \n"
  "    NONSPECIFIC_CURRENT ih \n"
  "    RANGE gnabar,gkbar,gnap_bar,ghbar,ena,ek,mhf,mhs,gh,gna\n"
  "    GLOBAL hf_tau_input,hs_tau_input\n"
  "    \n"
  "\n"
  "      \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    gnap_bar = 0.0005 (S/cm2)\n"
  "    gnabar = 0.052 (S/cm2) \n"
  "    gkbar = 0.011 (S/cm2)\n"
  "\n"
  "    ghbar = 0.0015 (S/cm2)\n"
  "    ena = 55 (mV) :reset by neuron, set after initialization\n"
  "    ek = -90 (mV) :reset by neuron, set after initialization\n"
  "    el = -65 (mV)\n"
  "    gl = 0.0005 (S/cm2)\n"
  "    eh = -20 (mV)\n"
  "    hf_tau_input=0.51 (ms)\n"
  "    hs_tau_input=5.6 (ms)\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "}   \n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "\n"
  "	    gna (S/cm2)\n"
  "	    gk (S/cm2)\n"
  "        gh (S/cm2)\n"
  "        gnap(S/cm2)\n"
  "        ina (mA/cm2)\n"
  "        ik (mA/cm2)\n"
  "        il (mA/cm2)\n"
  "        ih (mA/cm2)\n"
  "\n"
  "    amna (1/ms) bmna (1/ms) ahna (1/ms) bhna (1/ms) ank (1/ms) bnk (1/ms) mnap_beta(1/ms) mnap_alpha (1/ms)\n"
  "    mhfinf mhsinf  \n"
  "    mhstau (ms)  mhftau (ms)   \n"
  "} \n"
  "\n"
  "\n"
  "STATE {\n"
  "    mna hna mnap nk mhf mhs \n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT { \n"
  "    SOLVE states METHOD cnexp\n"
  "\n"
  "    gna = (gnabar*mna*mna*mna*hna)\n"
  "\n"
  "    gnap=(gnap_bar*mnap)\n"
  "	ina = (gna+gnap)*(v - ena)\n"
  "\n"
  "    gk = (gkbar*nk*nk*nk*nk)\n"
  "	ik = gk*(v - ek)      \n"
  "    il = gl*(v - el)\n"
  "    gh = ghbar*(0.65*mhf+0.35*mhs)\n"
  "    ih = gh*(v-eh)\n"
  "    \n"
  "\n"
  "\n"
  "\n"
  "}\n"
  "INITIAL {\n"
  "    rates(v)\n"
  "    mna = amna/(amna+bmna)\n"
  "    hna = ahna/(ahna+bhna)\n"
  "    mnap = mnap_alpha/(mnap_alpha+mnap_beta)\n"
  "\n"
  "    nk = ank/(ank+bnk)\n"
  "\n"
  "    mhf = mhfinf\n"
  "    mhs = mhsinf\n"
  "    \n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates(v)\n"
  "\n"
  "    mna' = amna*(1-mna) - bmna*mna\n"
  "    hna' = ahna*(1-hna) - bhna*hna\n"
  "    mnap' = mnap_alpha*(1-mnap) - mnap_beta*mnap\n"
  "    nk' = ank*(1-nk) - bnk*nk\n"
  "\n"
  "    mhf' = (mhfinf-mhf)/(mhftau)\n"
  "    mhs' = (mhsinf-mhs)/(mhstau)\n"
  "    \n"
  "\n"
  "\n"
  "}\n"
  "UNITSOFF\n"
  "PROCEDURE rates(v (mV)) {\n"
  "    amna = (.1)*vtrap(-(v+23),10)\n"
  "    bmna = 4*exp(-(v+48)/18)\n"
  "    ahna = 0.07*exp(-(v+37)/20)\n"
  "    bhna = 1/(exp(-0.1*(v+7))+1)\n"
  "    mnap_alpha=1/(0.15*(exp(-(v+38)/6.5)+1))\n"
  "    mnap_beta = (exp(-(v+38)/6.5))/(0.15*(exp(-(v+38)/6.5)+1))\n"
  "    ank = 0.01*vtrap(-(v+27),10)\n"
  "    bnk = 0.125 * exp(-(v+37)/80)\n"
  "    mhfinf = 1/(1+exp((v+79.2)/9.78))\n"
  "    mhftau = ((hf_tau_input / ((exp((v-1.7)/10)) + exp(-(v+340)/52))) + 1) :0.51\n"
  "    mhsinf = 1/(1+exp((v+71.3)/7.9))\n"
  "    mhstau = ((hs_tau_input/ ((exp((v-1.7)/14)) + exp(-(v+260)/43))) + 1 ) :5.6\n"
  "    \n"
  "    \n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{\n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
