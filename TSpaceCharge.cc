/////////////////////////////////////////////////////////////////////////////
//
// FILE NAME
//    TSpaceCharge.cc
//
// AUTHOR
//    John Galambos
//    ORNL, jdg@ornl.gov
//
// CREATED
//    5/20/98
//
//  DESCRIPTION
//    Module descriptor file for the Transverse Space Charge module.
//   This module contains source for the TSpaceCharge related info.
//
//  REVISION HISTORY
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
// include files
//
/////////////////////////////////////////////////////////////////////////////

#include "TSpaceCharge.h"
#include "Node.h"
#include "MapBase.h"
#include "TransMapHead.h"
#include "MacroPart.h"
#include "RealMat.h"
#include "ComplexMat.h"
#include "SCLTypes.h"
#include "StreamType.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fftw.h>
#include <rfftw.h>
#include "mpi.h"
#include "Boundary.h"
#include "FastMults.h"
#include "SC_Ellipt.h"
#include "SC_Ellipt_SBS.h"
#include "SC_CylCoords.h"
#include "SC_UprightEnv.h"
#include "SC_TiltEnv.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// STATIC DEFINITIONS
//
///////////////////////////////////////////////////////////////////////////
#define FORCETAG1 33
#define FORCETAG2 34
#define GRIDTAG1 31
#define GRIDTAG2 32

Array(ObjectPtr) TSpaceChargePointers;
extern Array(ObjectPtr) mParts, nodes, tMaps;
extern Array(ObjectPtr) syncP;
Real rClassical = 1.534698e-18; // Classical p radius

// Prototypes, globals for parallel routines:

Integer syncTSCBins(Real &xGridMin, Real &xGridMax,
                 Real &yGridMin, Real &yGridMax, Real &phiMin,
                 Real &phiMax, Integer &nMacros);
Void syncTSCRho1(Integer &nPoints);
Vector(Real) globalRho, localRho;

///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//    TransSC
//
// INHERITANCE RELATIONSHIPS
//    TransSC -> Node* -> Object
//
// USING/CONTAINING RELATIONSHIPS
//    None.
//
// DESCRIPTION
//    This is a base class for transverse  space charge implementations.
//
// PUBLIC MEMBERS
//   _nMacrosMin  Minimum number of macro particles before using this node.
//   _eps         Smoothing parameter for force calculation, fraction of
//                the bin size, i.e. _eps = 1 => smooth over entire bin
//                If _eps<0, use |_eps| as absoulte smoothing length [mm]
// PIC stuff:
//   _nXBins      Reference to the number of horizontal bins.
//   _nYBins      Reference to the number of vertical bins.
//   _lkick       Longitudinal length over which kick is applied (m)
//   _xGrid       The horizontal grid values relative to the closed orbit (mm)
//   _yGrid       The vertical grid values relative to the closed orbit (mm)
//   _phisc       The potential function over the grid
//   _fscx        The horizontal space charge force over the grid
//   _fscy        The vertical space charge force over the grid
//   _lambda      The average line charge / # of real particles (1/m)
//   _perveance   The perveance of the beam (using average line charge)
//   _resid       Residual of least squares fit for beam pipe potential.
//
// PROTECTED MEMBERS
//  None
// PRIVATE MEMBERS
//    None.
//
///////////////////////////////////////////////////////////////////////////

class TransSC : public Node {
  Declare_Standard_Members(TransSC, Node);
public:
    TransSC(const String &n, const Integer &order,
                 Integer &nXBins, Integer &nYBins, Real &length,
                 const Integer &nMacrosMin, Real &eps) :
            Node(n, 0., order), _nXBins(nXBins), _nYBins(nYBins),
            _lkick(length), _nMacrosMin(nMacrosMin), _eps(eps)

        {
                _xGrid.resize(_nXBins+1);
                _yGrid.resize(_nYBins+1);
                _phisc.resize(_nXBins+1, _nYBins+1);
                _fscx.resize(_nXBins+1, _nYBins+1);
                _fscy.resize(_nXBins+1, _nYBins+1);
        }

  ~TransSC()
  {
  }
    virtual Void nameOut(String &wname) { wname = _name; }
    virtual Void _updatePartAtNode(MacroPart &mp)=0;
    virtual Void _nodeCalculator(MacroPart &mp)=0;

    Integer _nXBins;
    Integer _nYBins;
    Real _lkick;
    Real _eps;
    const Integer _nMacrosMin;
    Real _perveance, _lambda;
    static Matrix(Real) _phisc, _fscx, _fscy;
    static Vector(Real) _xGrid, _yGrid;
    Real _resid;
 protected:
    Real _dx, _dy;
    Real _xGridMin, _xGridMax, _yGridMin, _yGridMax;
    Real _LFactor;
};

///////////////////////////////////////////////////////////////////////////
//
// NON-INLINE PROTECTED MEMBER FUNCTIONS FOR CLASS TransSC
//
///////////////////////////////////////////////////////////////////////////
Define_Standard_Members(TransSC, Node);
Matrix(Real) TransSC::_phisc, TransSC::_fscx, TransSC::_fscy;
Vector(Real) TransSC::_xGrid, TransSC::_yGrid;














///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//
//   UprightEnvSC
//
// INHERITANCE RELATIONSHIPS
//   UprightEnvSC -> TransSC -> Node* -> Object
//
// USING/CONTAINING RELATIONSHIPS
//   None.
//
// DESCRIPTION
//   This is a class for upright envelope evolution and space charge
//   calculation with beam pipe and uniform elliptical space charge
//   distribution. It uses a least squares solution to satisfy the beam
//   pipe conducting wall boundary condition.
//
//   The parameters BP1, BP2, BP3, BP4 have different meanings for different
//   shapes:
//   shape = "Ellipse" BP1 is 1/2 horizontal size, BP2 - 1/2 vertical size
//   shape = "Circle" BP1 - radius
//   shape = "Rectangle" BP1 - minX, BP2 - maxX, BP3 - minY, BP4 - maxY
//   Shape should be "Ellipse","Rectangle", "Circle", "None"
//   Sizes in [mm]
//
///////////////////////////////////////////////////////////////////////////

class UprightEnvSC : public TransSC
{
  Declare_Standard_Members(UprightEnvSC, TransSC);

  public:
  UprightEnvSC(const String &n, const Integer &order,
               Real &length, Real &nEnvelope,
               Real &Emitx, Real &Emity, 
               Real &Etax, Real &sigmaE,
               const String &BPShape,
               Real &BP1, Real &BP2, Real &BP3, Real &BP4,
               Integer &BPPoints, Integer &BPModes,
               Real &BPResid);
  ~UprightEnvSC();
  Void _updatePartAtNode(MacroPart &mp);
  Void _nodeCalculator(MacroPart &mp);
  Void _Invert(double** M, double** MINV, int n);

  static double _nEnvelope;
  static double _Emitx;
  static double _Emity;
  static double _Etax;
  static double _sigmaE;

  static String _BPShape;
  static int _BPPoints, _BPModes;
  static double   _BPResid;
  static double   _BPrnorm;
  static double*  _BPx;
  static double*  _BPy;
  static double** _BPPhiH1;
  static double** _BPPhiH2;

  static double _BPr, _BPa, _BPb, _BPxMin, _BPxMax, _BPyMin, _BPyMax;

  private:
  static UprightEnvCalculator* _UprightEnv;
  static int _initialized;
  static double _EnvTune_x, _EnvTune_y;
  static double _x4m, _y4m, _x3m, _y3m;
  static double _x2m, _y2m, _x1m, _y1m;
  static int _offx, _offy;
};

///////////////////////////////////////////////////////////////////////////
//
// Static MEMBERS FOR CLASS UprightEnvSC
//
///////////////////////////////////////////////////////////////////////////

  double UprightEnvSC::_nEnvelope;
  double UprightEnvSC::_Emitx;
  double UprightEnvSC::_Emity;
  double UprightEnvSC::_Etax;
  double UprightEnvSC::_sigmaE;

  String UprightEnvSC::_BPShape;
  int UprightEnvSC::_BPPoints;
  int UprightEnvSC::_BPModes;
  double   UprightEnvSC::_BPResid;
  double   UprightEnvSC::_BPrnorm;

  double*  UprightEnvSC::_BPx;
  double*  UprightEnvSC::_BPy;
  double** UprightEnvSC::_BPPhiH1;
  double** UprightEnvSC::_BPPhiH2;

  double UprightEnvSC::_BPr;
  double UprightEnvSC::_BPa;
  double UprightEnvSC::_BPb;
  double UprightEnvSC::_BPxMin;
  double UprightEnvSC::_BPxMax;
  double UprightEnvSC::_BPyMin;
  double UprightEnvSC::_BPyMax;

  UprightEnvCalculator* UprightEnvSC::_UprightEnv;

  int UprightEnvSC::_initialized  = 0;
  double UprightEnvSC::_EnvTune_x = 0.0;
  double UprightEnvSC::_EnvTune_y = 0.0;
  double UprightEnvSC::_x4m       = 0.0;
  double UprightEnvSC::_y4m       = 0.0;
  double UprightEnvSC::_x3m       = 0.0;
  double UprightEnvSC::_y3m       = 0.0;
  double UprightEnvSC::_x2m       = 0.0;
  double UprightEnvSC::_y2m       = 0.0;
  double UprightEnvSC::_x1m       = 0.0;
  double UprightEnvSC::_y1m       = 0.0;
  int UprightEnvSC::_offx = 10;
  int UprightEnvSC::_offy = 10;


///////////////////////////////////////////////////////////////////////////
//
// NON-INLINE PROTECTED MEMBER FUNCTIONS FOR CLASS UprightEnvSC
//
///////////////////////////////////////////////////////////////////////////

Define_Standard_Members(UprightEnvSC, TransSC);

UprightEnvSC::UprightEnvSC(const String &n, const Integer &order,
                           Real &length, Real &nEnvelope,
                           Real &Emitx, Real &Emity, 
                           Real &Etax, Real &sigmaE,
                           const String &BPShape,
                           Real &BP1, Real &BP2, Real &BP3, Real &BP4,
                           Integer &BPPoints, Integer &BPModes,
                           Real &BPResid) :
                           TransSC(n, order, BPModes, BPModes, length,
                           BPModes, BPResid)
{
  if(_initialized) return;

  _nEnvelope = double(nEnvelope);
  _Emitx     = double(Emitx);
  _Emity     = double(Emity);
  _Etax      = double(Etax);
  _sigmaE    = double(sigmaE);

  _BPShape  = BPShape;
  _BPPoints = 4 * (int(BPPoints) / 4);
  _BPModes  = int(BPModes);
  if(_BPShape == "None")
  {
    _BPPoints = 0;
    _BPModes  = 0;
  }
  _BPResid  = double(BPResid);
  _BPx      = new double[_BPPoints];
  _BPy      = new double[_BPPoints];
  _BPPhiH1  = new double*[_BPPoints];
  _BPPhiH2  = new double*[2 * _BPModes + 1];

  int i;

  for(i = 0; i < _BPPoints; i++)
  {
    _BPPhiH1[i] = new double[2 * _BPModes + 1];
  }

  for(i = 0; i < (2 * _BPModes + 1); i++)
  {
    _BPPhiH2[i] = new double[2 * _BPModes + 1];
  }

  double* rad;
  double* th;
  rad = new double[_BPPoints];
  th  = new double[_BPPoints];

  if(_BPShape == "Circle")
  {
    _BPr = double(BP1);
    _BPrnorm = _BPr;
    int i;
    double dtheta = 2. * Consts::pi / _BPPoints;
    for (i = 0; i < _BPPoints; i++)
    {
      rad[i] = _BPr;
      th[i]   = i * dtheta;
      _BPx[i] = _BPr * cos(th[i]);
      _BPy[i] = _BPr * sin(th[i]);
    }
  }

  if(_BPShape == "Ellipse")
  {
    _BPa = double(BP1);
    _BPb = double(BP2);
    double BPrsq = _BPa * _BPb;
    double BPasq = _BPa * _BPa;
    double BPbsq = _BPb * _BPb;
    _BPrnorm = pow(BPrsq, 0.5);
    int BPPO4 = _BPPoints / 4;

    double ds = 2. * double(Consts::pi) * _BPrnorm / _BPPoints;
    double rsq, term1, term2, dtheta, theta;
    double resid = 1.0;

    while(resid > 1.0e-08 || resid < -1.0e-08)
    {
      term1    = BPbsq;
      term2    = 0.0;
      rsq      = BPasq;
      dtheta   = ds / pow(rsq * (1.0 + term2 / (term1 * term1)), 0.5);
      theta    = dtheta / 2.;
      th[0] = 0.;
      for(i = 0; i < BPPO4; i++)
      {
        double sn2 = sin(theta) * sin(theta);
        term1 = BPbsq + (BPasq - BPbsq) * sn2;
        term2 = (BPbsq - BPasq) * (BPbsq - BPasq) * sn2 * (1.0 - sn2);
        rsq   = BPrsq * BPrsq / term1;
        dtheta = ds / pow(rsq * (1.0 + term2 / (term1 * term1)), 0.5);
        th[i + 1] = th[i] + dtheta;
        theta += dtheta;
      }
      resid = th[BPPO4] - double(Consts::pi) / 2.;
      ds *= double(Consts::pi) / (2. * th[BPPO4]);
    }

    int i1, i2, i3;

    for(i = 0; i < BPPO4; i++)
    {
      i1 = i;
      i2 = 2 * BPPO4 - i1;
      th[i2] = double(Consts::pi) - th[i1];
    }

    for (i = 1; i < 2 * BPPO4; i++)
    {
      i1 = i;
      i3 = 2 * BPPO4 + i1;
      th[i3] = double(Consts::pi) + th[i1];
    }

    for(i = 0; i < _BPPoints; i++)
    {
      term1   = BPbsq + (BPasq - BPbsq) * sin(th[i]) * sin(th[i]);
      rad[i]  = pow(BPrsq * BPrsq / term1, 0.5);
      _BPx[i] = rad[i] * cos(th[i]);
      _BPy[i] = rad[i] * sin(th[i]);
    }
  }

  if(_BPShape == "Rectangle")
  {
    _BPxMin = double(BP1);
    _BPxMax = double(BP2);
    _BPyMin = double(BP3);
    _BPyMax = double(BP4);
    double rsq = ((_BPxMax - _BPxMin) * (_BPxMax - _BPxMin) +
                 (_BPyMax - _BPyMin) * (_BPyMax - _BPyMin)) / 4.;
    _BPrnorm   = pow(rsq, 0.5);
    int BPPO4  = _BPPoints / 4;
    double dx  = (_BPxMax - _BPxMin) / BPPO4;
    double dy  = (_BPyMax - _BPyMin) / BPPO4;
    int i1, i2, i3, i4;
    for (i = 0; i < BPPO4; i++)
    {
      i1 = i;
      _BPx[i1]   = _BPxMax;
      _BPy[i1]   = _BPyMin + i * dy;
      rad[i1]    = pow(_BPx[i1] * _BPx[i1] + _BPy[i1] * _BPy[i1], 0.5);
      th[i1]     = atan2(_BPy[i1], _BPx[i1]);

      i2 = BPPO4 + i1;
      _BPx[i2]   = _BPxMax - i * dx;
      _BPy[i2]   = _BPyMax;
      rad[i2]    = pow(_BPx[i2] * _BPx[i2] + _BPy[i2] * _BPy[i2], 0.5);
      th[i2]     = atan2(_BPy[i2], _BPx[i2]);

      i3 = BPPO4 + i2;
      _BPx[i3]   = _BPxMin;
      _BPy[i3]   = _BPyMax - i * dy;
      rad[i3]    = pow(_BPx[i3] * _BPx[i3] + _BPy[i3] * _BPy[i3], 0.5);
      th[i3]     = atan2(_BPy[i3], _BPx[i3]);

      i4 = BPPO4 + i3;
      _BPx[i4]   = _BPxMin + i * dx;
      _BPy[i4]   = _BPyMin;
      rad[i4]    = pow(_BPx[i4] * _BPx[i4] + _BPy[i4] * _BPy[i4], 0.5);
      th[i4]     = atan2(_BPy[i4], _BPx[i4]);
    }
  }

  if(_BPShape != "None")
  {
    double rfac, rj;
    int j, jc, js;

    for(i = 0; i < _BPPoints; i++)
    {
      rfac = rad[i] / _BPrnorm;
      rj = 1.0;
      _BPPhiH1[i][0] = 1.0;
      for (j = 1; j <= _BPModes; j++)
      {
        jc = 2 * j - 1;
        js = 2 * j;
        rj *= rfac;
        _BPPhiH1[i][jc] = rj * cos(j * th[i]);
        _BPPhiH1[i][js] = rj * sin(j * th[i]);
      }
    }

    double** TEMP;
    TEMP = new double*[2 * _BPModes + 1];
    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      TEMP[i] = new double[2 * _BPModes + 1];
    }

    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      for(j = 0; j < (2 * _BPModes + 1); j++)
      {
        int k;
        TEMP[i][j] = 0.0;
        for(k = 0; k < _BPPoints; k++)
        {
          TEMP[i][j] += _BPPhiH1[k][i] * _BPPhiH1[k][j];
        }
      }
    }

    _Invert(TEMP, _BPPhiH2, 2 * _BPModes + 1);

    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      delete TEMP[i];
    }
    delete TEMP;
  }

  delete rad;
  delete th;

  /*
  for(i = 0; i < _BPPoints; i++)
  {
    cout << i << "  " << _BPx[i] << "  " << _BPy[i] << "\n";
  }
  */

  double ZLength = double(length);

  _UprightEnv = new UprightEnvCalculator(_Emitx, _Emity,
                                         _Etax, _sigmaE,
                                         _BPPoints, _BPModes,
                                         _BPrnorm, _BPResid,
                                         _BPx , _BPy,
                                         _BPPhiH1, _BPPhiH2);

  _initialized = 1;
}

UprightEnvSC::~UprightEnvSC()
{
  delete _BPx;
  delete _BPy;
  int i;
  for(i = 0; i < _BPPoints; i++)
  {
    delete _BPPhiH1[i];
  }
  for(i = 0; i < 2 * _BPModes + 1; i++)
  {
    delete _BPPhiH2[i];
  }
  delete _BPPhiH1;
  delete _BPPhiH2;
  delete _UprightEnv;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    UprightEnvSC::NodeCalculator
//
// DESCRIPTION
//    Kicks the envelope and sets up boundary force multipliers
//
// PARAMETERS
//    None.
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

Void UprightEnvSC::_nodeCalculator(MacroPart &mp)
{
  int i_flag;
  MPI_Initialized(&i_flag);

  if(i_flag)
  {
    MPI_Allreduce(&mp._nMacros, &mp._globalNMacros, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  else
  {
    mp._globalNMacros = mp._nMacros;
  }

  double _phiMax_tmp_local = mp._phiMax;
  double _phiMin_tmp_local = mp._phiMin;
  double _phiMax_tmp_global;
  double _phiMin_tmp_global;

  if(i_flag)
  {
    MPI_Allreduce(&_phiMax_tmp_local, &_phiMax_tmp_global, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&_phiMin_tmp_local, &_phiMin_tmp_global, 1,
                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else
  {
    mp._globalNMacros  = mp._nMacros;
    _phiMax_tmp_global = _phiMax_tmp_local;
    _phiMin_tmp_global = _phiMin_tmp_local;
  }

  double lambda = double(_nEnvelope) * Ring::harmonicNumber / Ring::lRing;

  double Intensity = 1000000. * 2.0 * lambda *
                     mp._syncPart._charge * mp._syncPart._charge *
                     rClassical /
                    (mp._syncPart._betaSync *
                     mp._syncPart._betaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync * mp._syncPart._mass);

  _UprightEnv->_KickEnvelope(mp, Intensity, double(_lkick));
  if(_BPShape != "None")
  {
    _UprightEnv->_BoundaryCoeffs();
  }

  OFstream fio("Envelope.out", ios::app);

  if(((_x2m - mp._x(1)) * (_x2m - _x4m) > 0.0) && _offx > 5)
  {
    _EnvTune_x += 1.0;
    _offx = 0;
  }
  if(((_y2m - mp._y(1)) * (_y2m - _y4m) > 0.0) && _offy > 5)
  {
    _EnvTune_y += 1.0;
    _offy = 0;
  }
  _x4m = _x3m;
  _y4m = _y3m;
  _x3m = _x2m;
  _y3m = _y2m;
  _x2m = _x1m;
  _y2m = _y1m;
  _x1m = mp._x(1);
  _y1m = mp._y(1);
  _offx += 1;
  _offy += 1;

  fio << Ring::nTurnsDone << "\t" << _position << "\t"
      << Ring::nTurnsDone * Ring::lRing + _position << "\t";
  fio << mp._x(1)   << "\t" << mp._xp(1)  << "\t"
      << mp._y(1)   << "\t" << mp._yp(1)  << "\t"
      << _EnvTune_x << "\t" << _EnvTune_y << "\n";

  fio.close();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    UprightEnvSC::_updatePartAtNode
//
// DESCRIPTION
//    Calls the specified local calculator for an operation on
//    a MacroParticle with a UprightEnvSC. The Transverse space
//    charge kick is added to each macro particle here.
//
// PARAMETERS
//    None.
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

Void UprightEnvSC::_updatePartAtNode(MacroPart &mp)
{

  double lambda = double(_nEnvelope) * Ring::harmonicNumber / Ring::lRing;

  double Intensity = 1000000. * 2.0 * lambda *
                     mp._syncPart._charge * mp._syncPart._charge *
                     rClassical /
                    (mp._syncPart._betaSync *
                     mp._syncPart._betaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync * mp._syncPart._mass);

  _UprightEnv->_ApplyForce(mp, Intensity, double(_lkick));
  _UprightEnv->_Clear();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   UprightEnvSC::_Invert
//
// DESCRIPTION
//   Inverts nxn matrix M returning result as MINV
//
///////////////////////////////////////////////////////////////////////////

Void UprightEnvSC::_Invert(double** M, double** MINV, int n)
{
  int i, j, k, ipvt;
  double Max, Coeff;

  double* Temp;
  Temp = new double[n];

  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      MINV[i][j] = 0.0;
    }
    MINV[i][i] = 1.0;
  }

  for(j = 0; j < n; j++)
  {
    Max = 0.0;
    ipvt = j;

    for(i = j; i < n; i++)
    {
      if(M[i][j] > Max)
      {
        Max = M[i][j];
        ipvt = i;
      }
      if(-M[i][j] > Max)
      {
        Max = -M[i][j];
        ipvt = i;
      }
    }

    if(ipvt != j)
    {
      for(k = 0; k < n; k++)
      {
        Temp[k]       = M[ipvt][k];
        M[ipvt][k]    = M[j][k];
        M[j][k]       = Temp[k];
        Temp[k]       = MINV[ipvt][k];
        MINV[ipvt][k] = MINV[j][k];
        MINV[j][k]    = Temp[k];
      }
    }

    for(i = 0; i < n; i++)
    {
      if(i != j)
      {
        Coeff = M[i][j] / M[j][j];
        for(k = j; k < n; k++)
        {
          M[i][k] -= Coeff * M[j][k];
        }
        M[i][j] = 0.0;
        for(k = 0; k < n; k++)
        {
          MINV[i][k] -= Coeff * MINV[j][k];
        }
      }
    }

    Coeff = M[j][j];
    for(k = j; k < n; k++)
    {
      M[j][k] /= Coeff;
    }
    for(k = 0; k < n; k++)
    {
      MINV[j][k] /= Coeff;
    }
  }

  delete Temp;
}














///////////////////////////////////////////////////////////////////////////
//
// CLASS NAME
//
//   TiltEnvSC
//
// INHERITANCE RELATIONSHIPS
//   TiltEnvSC -> TransSC -> Node* -> Object
//
// USING/CONTAINING RELATIONSHIPS
//   None.
//
// DESCRIPTION
//   This is a class for upright envelope evolution and space charge
//   calculation with beam pipe and uniform elliptical space charge
//   distribution. It uses a least squares solution to satisfy the beam
//   pipe conducting wall boundary condition.
//
//   The parameters BP1, BP2, BP3, BP4 have different meanings for different
//   shapes:
//   shape = "Ellipse" BP1 is 1/2 horizontal size, BP2 - 1/2 vertical size
//   shape = "Circle" BP1 - radius
//   shape = "Rectangle" BP1 - minX, BP2 - maxX, BP3 - minY, BP4 - maxY
//   Shape should be "Ellipse","Rectangle", "Circle", "None"
//   Sizes in [mm]
//
///////////////////////////////////////////////////////////////////////////

class TiltEnvSC : public TransSC
{
  Declare_Standard_Members(TiltEnvSC, TransSC);

  public:
  TiltEnvSC(const String &n, const Integer &order,
            Real &length, Real &nEnvelope,
            const String &BPShape,
            Real &BP1, Real &BP2, Real &BP3, Real &BP4,
            Integer &BPPoints, Integer &BPModes,
            Real &BPResid);
  ~TiltEnvSC();
  Void _updatePartAtNode(MacroPart &mp);
  Void _nodeCalculator(MacroPart &mp);
  Void _Invert(double** M, double** MINV, int n);

  static double _nEnvelope;

  static String _BPShape;
  static int _BPPoints, _BPModes;
  static double   _BPResid;
  static double   _BPrnorm;
  static double*  _BPx;
  static double*  _BPy;
  static double** _BPPhiH1;
  static double** _BPPhiH2;

  static double _BPr, _BPa, _BPb, _BPxMin, _BPxMax, _BPyMin, _BPyMax;

  private:
  static TiltEnvCalculator* _TiltEnv;
  static int _initialized;
  static double _EnvTune_a, _EnvTune_b;
  static double _axy4m, _bxy4m, _axy3m, _bxy3m;
  static double _axy2m, _bxy2m, _axy1m, _bxy1m;
  static int _offa, _offb;
};

///////////////////////////////////////////////////////////////////////////
//
// Static MEMBERS FOR CLASS TiltEnvSC
//
///////////////////////////////////////////////////////////////////////////

  double TiltEnvSC::_nEnvelope;

  String TiltEnvSC::_BPShape;
  int TiltEnvSC::_BPPoints;
  int TiltEnvSC::_BPModes;
  double   TiltEnvSC::_BPResid;
  double   TiltEnvSC::_BPrnorm;

  double*  TiltEnvSC::_BPx;
  double*  TiltEnvSC::_BPy;
  double** TiltEnvSC::_BPPhiH1;
  double** TiltEnvSC::_BPPhiH2;

  double TiltEnvSC::_BPr;
  double TiltEnvSC::_BPa;
  double TiltEnvSC::_BPb;
  double TiltEnvSC::_BPxMin;
  double TiltEnvSC::_BPxMax;
  double TiltEnvSC::_BPyMin;
  double TiltEnvSC::_BPyMax;

  TiltEnvCalculator* TiltEnvSC::_TiltEnv;

  int TiltEnvSC::_initialized = 0;
  double TiltEnvSC::_EnvTune_a = 0.0;
  double TiltEnvSC::_EnvTune_b = 0.0;
  double TiltEnvSC::_axy4m     = 0.0;
  double TiltEnvSC::_bxy4m     = 0.0;
  double TiltEnvSC::_axy3m     = 0.0;
  double TiltEnvSC::_bxy3m     = 0.0;
  double TiltEnvSC::_axy2m     = 0.0;
  double TiltEnvSC::_bxy2m     = 0.0;
  double TiltEnvSC::_axy1m     = 0.0;
  double TiltEnvSC::_bxy1m     = 0.0;
  int TiltEnvSC::_offa = 10;
  int TiltEnvSC::_offb = 10;


///////////////////////////////////////////////////////////////////////////
//
// NON-INLINE PROTECTED MEMBER FUNCTIONS FOR CLASS TiltEnvSC
//
///////////////////////////////////////////////////////////////////////////

Define_Standard_Members(TiltEnvSC, TransSC);

TiltEnvSC::TiltEnvSC(const String &n, const Integer &order,
                     Real &length, Real &nEnvelope,
                     const String &BPShape,
                     Real &BP1, Real &BP2, Real &BP3, Real &BP4,
                     Integer &BPPoints, Integer &BPModes,
                     Real &BPResid) :
                     TransSC(n, order, BPModes, BPModes, length,
                     BPModes, BPResid)
{
  if(_initialized) return;

  _nEnvelope = double(nEnvelope);

  _BPShape  = BPShape;
  _BPPoints = 4 * (int(BPPoints) / 4);
  _BPModes  = int(BPModes);
  if(_BPShape == "None")
  {
    _BPPoints = 0;
    _BPModes  = 0;
  }
  _BPResid  = double(BPResid);
  _BPx      = new double[_BPPoints];
  _BPy      = new double[_BPPoints];
  _BPPhiH1  = new double*[_BPPoints];
  _BPPhiH2  = new double*[2 * _BPModes + 1];

  int i;

  for(i = 0; i < _BPPoints; i++)
  {
    _BPPhiH1[i] = new double[2 * _BPModes + 1];
  }

  for(i = 0; i < (2 * _BPModes + 1); i++)
  {
    _BPPhiH2[i] = new double[2 * _BPModes + 1];
  }

  double* rad;
  double* th;
  rad = new double[_BPPoints];
  th  = new double[_BPPoints];

  if(_BPShape == "Circle")
  {
    _BPr = double(BP1);
    _BPrnorm = _BPr;
    int i;
    double dtheta = 2. * Consts::pi / _BPPoints;
    for (i = 0; i < _BPPoints; i++)
    {
      rad[i] = _BPr;
      th[i]   = i * dtheta;
      _BPx[i] = _BPr * cos(th[i]);
      _BPy[i] = _BPr * sin(th[i]);
    }
  }

  if(_BPShape == "Ellipse")
  {
    _BPa = double(BP1);
    _BPb = double(BP2);
    double BPrsq = _BPa * _BPb;
    double BPasq = _BPa * _BPa;
    double BPbsq = _BPb * _BPb;
    _BPrnorm = pow(BPrsq, 0.5);
    int BPPO4 = _BPPoints / 4;

    double ds = 2. * double(Consts::pi) * _BPrnorm / _BPPoints;
    double rsq, term1, term2, dtheta, theta;
    double resid = 1.0;

    while(resid > 1.0e-08 || resid < -1.0e-08)
    {
      term1    = BPbsq;
      term2    = 0.0;
      rsq      = BPasq;
      dtheta   = ds / pow(rsq * (1.0 + term2 / (term1 * term1)), 0.5);
      theta    = dtheta / 2.;
      th[0] = 0.;
      for(i = 0; i < BPPO4; i++)
      {
        double sn2 = sin(theta) * sin(theta);
        term1 = BPbsq + (BPasq - BPbsq) * sn2;
        term2 = (BPbsq - BPasq) * (BPbsq - BPasq) * sn2 * (1.0 - sn2);
        rsq   = BPrsq * BPrsq / term1;
        dtheta = ds / pow(rsq * (1.0 + term2 / (term1 * term1)), 0.5);
        th[i + 1] = th[i] + dtheta;
        theta += dtheta;
      }
      resid = th[BPPO4] - double(Consts::pi) / 2.;
      ds *= double(Consts::pi) / (2. * th[BPPO4]);
    }

    int i1, i2, i3;

    for(i = 0; i < BPPO4; i++)
    {
      i1 = i;
      i2 = 2 * BPPO4 - i1;
      th[i2] = double(Consts::pi) - th[i1];
    }

    for (i = 1; i < 2 * BPPO4; i++)
    {
      i1 = i;
      i3 = 2 * BPPO4 + i1;
      th[i3] = double(Consts::pi) + th[i1];
    }

    for(i = 0; i < _BPPoints; i++)
    {
      term1   = BPbsq + (BPasq - BPbsq) * sin(th[i]) * sin(th[i]);
      rad[i]  = pow(BPrsq * BPrsq / term1, 0.5);
      _BPx[i] = rad[i] * cos(th[i]);
      _BPy[i] = rad[i] * sin(th[i]);
    }
  }

  if(_BPShape == "Rectangle")
  {
    _BPxMin = double(BP1);
    _BPxMax = double(BP2);
    _BPyMin = double(BP3);
    _BPyMax = double(BP4);
    double rsq = ((_BPxMax - _BPxMin) * (_BPxMax - _BPxMin) +
                 (_BPyMax - _BPyMin) * (_BPyMax - _BPyMin)) / 4.;
    _BPrnorm   = pow(rsq, 0.5);
    int BPPO4  = _BPPoints / 4;
    double dx  = (_BPxMax - _BPxMin) / BPPO4;
    double dy  = (_BPyMax - _BPyMin) / BPPO4;
    int i1, i2, i3, i4;
    for (i = 0; i < BPPO4; i++)
    {
      i1 = i;
      _BPx[i1]   = _BPxMax;
      _BPy[i1]   = _BPyMin + i * dy;
      rad[i1]    = pow(_BPx[i1] * _BPx[i1] + _BPy[i1] * _BPy[i1], 0.5);
      th[i1]     = atan2(_BPy[i1], _BPx[i1]);

      i2 = BPPO4 + i1;
      _BPx[i2]   = _BPxMax - i * dx;
      _BPy[i2]   = _BPyMax;
      rad[i2]    = pow(_BPx[i2] * _BPx[i2] + _BPy[i2] * _BPy[i2], 0.5);
      th[i2]     = atan2(_BPy[i2], _BPx[i2]);

      i3 = BPPO4 + i2;
      _BPx[i3]   = _BPxMin;
      _BPy[i3]   = _BPyMax - i * dy;
      rad[i3]    = pow(_BPx[i3] * _BPx[i3] + _BPy[i3] * _BPy[i3], 0.5);
      th[i3]     = atan2(_BPy[i3], _BPx[i3]);

      i4 = BPPO4 + i3;
      _BPx[i4]   = _BPxMin + i * dx;
      _BPy[i4]   = _BPyMin;
      rad[i4]    = pow(_BPx[i4] * _BPx[i4] + _BPy[i4] * _BPy[i4], 0.5);
      th[i4]     = atan2(_BPy[i4], _BPx[i4]);
    }
  }

  if(_BPShape != "None")
  {
    double rfac, rj;
    int j, jc, js;

    for(i = 0; i < _BPPoints; i++)
    {
      rfac = rad[i] / _BPrnorm;
      rj = 1.0;
      _BPPhiH1[i][0] = 1.0;
      for (j = 1; j <= _BPModes; j++)
      {
        jc = 2 * j - 1;
        js = 2 * j;
        rj *= rfac;
        _BPPhiH1[i][jc] = rj * cos(j * th[i]);
        _BPPhiH1[i][js] = rj * sin(j * th[i]);
      }
    }

    double** TEMP;
    TEMP = new double*[2 * _BPModes + 1];
    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      TEMP[i] = new double[2 * _BPModes + 1];
    }

    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      for(j = 0; j < (2 * _BPModes + 1); j++)
      {
        int k;
        TEMP[i][j] = 0.0;
        for(k = 0; k < _BPPoints; k++)
        {
          TEMP[i][j] += _BPPhiH1[k][i] * _BPPhiH1[k][j];
        }
      }
    }

    _Invert(TEMP, _BPPhiH2, 2 * _BPModes + 1);

    for(i = 0; i < (2 * _BPModes + 1); i++)
    {
      delete TEMP[i];
    }
    delete TEMP;
  }

  delete rad;
  delete th;

  /*
  for(i = 0; i < _BPPoints; i++)
  {
    cout << i << "  " << _BPx[i] << "  " << _BPy[i] << "\n";
  }
  */

  double ZLength = double(length);

  _TiltEnv = new TiltEnvCalculator(_BPPoints, _BPModes,
                                   _BPrnorm, _BPResid,
                                   _BPx , _BPy,
                                   _BPPhiH1, _BPPhiH2);

  _initialized = 1;
}

TiltEnvSC::~TiltEnvSC()
{
  delete _BPx;
  delete _BPy;
  int i;
  for(i = 0; i < _BPPoints; i++)
  {
    delete _BPPhiH1[i];
  }
  for(i = 0; i < 2 * _BPModes + 1; i++)
  {
    delete _BPPhiH2[i];
  }
  delete _BPPhiH1;
  delete _BPPhiH2;
  delete _TiltEnv;
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    TiltEnvSC::NodeCalculator
//
// DESCRIPTION
//    Kicks the envelope and sets up boundary force multipliers
//
// PARAMETERS
//    None.
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

Void TiltEnvSC::_nodeCalculator(MacroPart &mp)
{
  int i_flag;
  MPI_Initialized(&i_flag);

  if(i_flag)
  {
    MPI_Allreduce(&mp._nMacros, &mp._globalNMacros, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  else
  {
    mp._globalNMacros = mp._nMacros;
  }

  double _phiMax_tmp_local = mp._phiMax;
  double _phiMin_tmp_local = mp._phiMin;
  double _phiMax_tmp_global;
  double _phiMin_tmp_global;

  if(i_flag)
  {
    MPI_Allreduce(&_phiMax_tmp_local, &_phiMax_tmp_global, 1,
                  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&_phiMin_tmp_local, &_phiMin_tmp_global, 1,
                  MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  else
  {
    mp._globalNMacros  = mp._nMacros;
    _phiMax_tmp_global = _phiMax_tmp_local;
    _phiMin_tmp_global = _phiMin_tmp_local;
  }

  double lambda = double(_nEnvelope) * Ring::harmonicNumber / Ring::lRing;

  double Intensity = 1000000. * 2.0 * lambda *
                     mp._syncPart._charge * mp._syncPart._charge *
                     rClassical /
                    (mp._syncPart._betaSync *
                     mp._syncPart._betaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync * mp._syncPart._mass);

  _TiltEnv->_KickEnvelope(mp, Intensity, double(_lkick));
  if(_BPShape != "None")
  {
    _TiltEnv->_BoundaryCoeffs();
  }

  double c1sq, c2sq;
  double tiltxy,   axy,   bxy;
  double tiltxxp,  axxp,  bxxp;
  double tiltyyp,  ayyp,  byyp;
  double tiltxpyp, axpyp, bxpyp;
  double tiltxyp,  axyp,  bxyp;
  double tiltyxp,  ayxp,  byxp;

  _TiltEnv->_getEllipsDat(mp._x(1), mp._y(1), mp._x(2), mp._y(2),
                          tiltxy, c1sq, c2sq);
  axy = pow(c1sq, 0.5);
  bxy = pow(c2sq, 0.5);

  if(((_axy2m - axy) * (_axy2m -_axy4m) > 0.0) && _offa > 5)
  {
    _EnvTune_a += 1.0;
    _offa = 0;
  }
  if(((_bxy2m - bxy) * (_bxy2m -_bxy4m) > 0.0) && _offb > 5)
  {
    _EnvTune_b += 1.0;
    _offb = 0;
  }
  _axy4m = _axy3m;
  _bxy4m = _bxy3m;
  _axy3m = _axy2m;
  _bxy3m = _bxy2m;
  _axy2m = _axy1m;
  _bxy2m = _bxy1m;
  _axy1m = axy;
  _bxy1m = bxy;
  _offa += 1;
  _offb += 1;

  _TiltEnv->_getEllipsDat(mp._x(1), mp._xp(1), mp._x(2), mp._xp(2),
                          tiltxxp, c1sq, c2sq);
  axxp = pow(c1sq, 0.5);
  bxxp = pow(c2sq, 0.5);

  _TiltEnv->_getEllipsDat(mp._y(1), mp._yp(1), mp._y(2), mp._yp(2),
                          tiltyyp, c1sq, c2sq);
  ayyp = pow(c1sq, 0.5);
  byyp = pow(c2sq, 0.5);

  _TiltEnv->_getEllipsDat(mp._xp(1), mp._yp(1), mp._xp(2), mp._yp(2),
                          tiltxpyp, c1sq, c2sq);
  axpyp = pow(c1sq, 0.5);
  bxpyp = pow(c2sq, 0.5);

  _TiltEnv->_getEllipsDat(mp._x(1), mp._yp(1), mp._x(2), mp._yp(2),
                          tiltxyp, c1sq, c2sq);
  axyp = pow(c1sq, 0.5);
  bxyp = pow(c2sq, 0.5);

  _TiltEnv->_getEllipsDat(mp._y(1), mp._xp(1), mp._y(2), mp._xp(2),
                          tiltyxp, c1sq, c2sq);
  ayxp = pow(c1sq, 0.5);
  byxp = pow(c2sq, 0.5);

  OFstream fio("Envelope.out", ios::app);

  fio << Ring::nTurnsDone << "\t" << _position << "\t"
      << Ring::nTurnsDone * Ring::lRing + _position << "\t";
  fio << mp._x(1)   << "\t" << mp._x(2)  << "\t"
      << mp._y(1)   << "\t" << mp._y(2)  << "\t"
      << mp._xp(1)  << "\t" << mp._xp(2) << "\t"
      << mp._yp(1)  << "\t" << mp._yp(2) << "\t"
      << tiltxy     << "\t" << axy       << "\t" << bxy   << "\t"
      << tiltxxp    << "\t" << axxp      << "\t" << bxxp  << "\t"
      << tiltyyp    << "\t" << ayyp      << "\t" << byyp  << "\t"
      << tiltxpyp   << "\t" << axpyp     << "\t" << bxpyp << "\t"
      << tiltxyp    << "\t" << axyp      << "\t" << bxyp  << "\t"
      << tiltyxp    << "\t" << ayxp      << "\t" << byxp  << "\t"
      << _EnvTune_a << "\t" << _EnvTune_b << "\n";

  fio.close();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//    TiltEnvSC::_updatePartAtNode
//
// DESCRIPTION
//    Calls the specified local calculator for an operation on
//    a MacroParticle with a TiltEnvSC. The Transverse space
//    charge kick is added to each macro particle here.
//
// PARAMETERS
//    None.
//
// RETURNS
//    Nothing.
//
///////////////////////////////////////////////////////////////////////////

Void TiltEnvSC::_updatePartAtNode(MacroPart &mp)
{

  double lambda = double(_nEnvelope) * Ring::harmonicNumber / Ring::lRing;

  double Intensity = 1000000. * 2.0 * lambda *
                     mp._syncPart._charge * mp._syncPart._charge *
                     rClassical /
                    (mp._syncPart._betaSync *
                     mp._syncPart._betaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync *
                     mp._syncPart._gammaSync * mp._syncPart._mass);

  _TiltEnv->_ApplyForce(mp, Intensity, double(_lkick));
  _TiltEnv->_Clear();
}

///////////////////////////////////////////////////////////////////////////
//
// NAME
//   TiltEnvSC::_Invert
//
// DESCRIPTION
//   Inverts nxn matrix M returning result as MINV
//
///////////////////////////////////////////////////////////////////////////

Void TiltEnvSC::_Invert(double** M, double** MINV, int n)
{
  int i, j, k, ipvt;
  double Max, Coeff;

  double* Temp;
  Temp = new double[n];

  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      MINV[i][j] = 0.0;
    }
    MINV[i][i] = 1.0;
  }

  for(j = 0; j < n; j++)
  {
    Max = 0.0;
    ipvt = j;

    for(i = j; i < n; i++)
    {
      if(M[i][j] > Max)
      {
        Max = M[i][j];
        ipvt = i;
      }
      if(-M[i][j] > Max)
      {
        Max = -M[i][j];
        ipvt = i;
      }
    }

    if(ipvt != j)
    {
      for(k = 0; k < n; k++)
      {
        Temp[k]       = M[ipvt][k];
        M[ipvt][k]    = M[j][k];
        M[j][k]       = Temp[k];
        Temp[k]       = MINV[ipvt][k];
        MINV[ipvt][k] = MINV[j][k];
        MINV[j][k]    = Temp[k];
      }
    }

    for(i = 0; i < n; i++)
    {
      if(i != j)
      {
        Coeff = M[i][j] / M[j][j];
        for(k = j; k < n; k++)
        {
          M[i][k] -= Coeff * M[j][k];
        }
        M[i][j] = 0.0;
        for(k = 0; k < n; k++)
        {
          MINV[i][k] -= Coeff * MINV[j][k];
        }
      }
    }

    Coeff = M[j][j];
    for(k = j; k < n; k++)
    {
      M[j][k] /= Coeff;
    }
    for(k = 0; k < n; k++)
    {
      MINV[j][k] /= Coeff;
    }
  }

  delete Temp;
}
