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
 
#define nrn_init _nrn_init__i_theta
#define _nrn_initial _nrn_initial__i_theta
#define nrn_cur _nrn_cur__i_theta
#define _nrn_current _nrn_current__i_theta
#define nrn_jacob _nrn_jacob__i_theta
#define nrn_state _nrn_state__i_theta
#define _net_receive _net_receive__i_theta 
 
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
#define Amp _p[0]
#define Amp_columnindex 0
#define dc _p[1]
#define dc_columnindex 1
#define c _p[2]
#define c_columnindex 2
#define itheta _p[3]
#define itheta_columnindex 3
#define v _p[4]
#define v_columnindex 4
#define _g _p[5]
#define _g_columnindex 5
 
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
 static void _hoc_dc_to_freq(void);
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
 "setdata_i_theta", _hoc_setdata,
 "dc_to_freq_i_theta", _hoc_dc_to_freq,
 0, 0
};
#define dc_to_freq dc_to_freq_i_theta
 extern double dc_to_freq( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define omega omega_i_theta
 double omega = 0.01;
#define phi phi_i_theta
 double phi = 0;
#define theta_vel_tuning theta_vel_tuning_i_theta
 double theta_vel_tuning = 0;
#define vthresh vthresh_i_theta
 double vthresh = -80;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "vthresh_i_theta", "mV",
 "omega_i_theta", "1/ms",
 "Amp_i_theta", "S/cm2",
 "itheta_i_theta", "mA/cm2",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vthresh_i_theta", &vthresh_i_theta,
 "omega_i_theta", &omega_i_theta,
 "phi_i_theta", &phi_i_theta,
 "theta_vel_tuning_i_theta", &theta_vel_tuning_i_theta,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"i_theta",
 "Amp_i_theta",
 "dc_i_theta",
 "c_i_theta",
 0,
 "itheta_i_theta",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 6, _prop);
 	/*initialize range parameters*/
 	Amp = 0.0001;
 	dc = 0;
 	c = 0;
 	_prop->param = _p;
 	_prop->param_size = 6;
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _itheta_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 6, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 i_theta /Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/itheta.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Theta for interneurons";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
double dc_to_freq ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _ldc_to_freq;
 if ( _ly  == 0.0 ) {
     _ldc_to_freq = omega ;
     }
   else {
     _ldc_to_freq = ( 8.03 ) * _lx + 0.02556 ;
     }
   
return _ldc_to_freq;
 }
 
static void _hoc_dc_to_freq(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  dc_to_freq ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{

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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   itheta = Amp * sin ( 2.0 * 3.14 * dc_to_freq ( _threadargscomma_ dc , theta_vel_tuning ) * ( t + phi ) ) + c ;
   }
 _current += itheta;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
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

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/govindrnair/TheoNeuroLab/GridCellsCond/mod/itheta.mod";
static const char* nmodl_file_text = 
  "TITLE Theta for interneurons\n"
  "\n"
  "\n"
  "\n"
  "UNITS {\n"
  "    (mV)=(millivolt)\n"
  "    (S) = (siemens)\n"
  "    (mA) = (milliamp)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX i_theta\n"
  "    NONSPECIFIC_CURRENT itheta \n"
  "    RANGE Amp,dc,c\n"
  "    GLOBAL vthresh,omega,phi,theta_vel_tuning\n"
  "    \n"
  "\n"
  "      \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    Amp = 1e-4 (S/cm2)\n"
  "    vthresh = -80 (mV)\n"
  "    omega = 0.01 (1/ms)\n"
  "    phi=0\n"
  "    dc=0\n"
  "    c=0\n"
  "    theta_vel_tuning=0}   \n"
  "\n"
  "ASSIGNED {\n"
  "        v (mV)\n"
  "        itheta (mA/cm2)\n"
  "\n"
  "} \n"
  "\n"
  "\n"
  "\n"
  "\n"
  "BREAKPOINT { \n"
  "\n"
  "	itheta = Amp*sin(2*3.14*dc_to_freq(dc,theta_vel_tuning)*(t+phi)) + c\n"
  "    \n"
  "\n"
  "}\n"
  "\n"
  "FUNCTION dc_to_freq(x,y) { \n"
  "  if (y==0) {\n"
  "    dc_to_freq=omega\n"
  "  } else {\n"
  "        dc_to_freq=(8.03)*x+0.02556\n"
  "}\n"
  "}\n"
  ;
#endif
