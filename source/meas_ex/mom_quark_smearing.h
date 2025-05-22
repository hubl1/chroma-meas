// -*- C++ -*-
// $Id: gaus_quark_smearing.h,v 3.2 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color std::vector and propagator
 */

#ifndef __mom_quark_smearing_h__
#define __mom_quark_smearing_h__

#include "chroma_config.h"

#include "meas/smear/quark_smearing.h"
#include "actions/boson/operator/klein_gord.h"

namespace Chroma {
//! Name and registration
/*! @ingroup smear */
namespace MomQuarkSmearingEnv {
bool registerAll();

//! Return the name
std::string getName();

//! Params for Gauss quark smearing
/*! @ingroup smear */
struct Params {
  Params() {
  }
  Params(XMLReader& in, const std::string& path);
  void writeXML(XMLWriter& in, const std::string& path) const;

  Real wvf_param; /*!< Smearing width */
  int wvfIntPar; /*!< Number of smearing hits */
  int no_smear_dir; /*!< No smearing in this direction */
  int smear_dir; /*!< Only apply smearing in this direction */
  int smear_tr_size; /*!< the size of operator in the transverse direction */
  double xi; /*!< the anisotropic lattice asymmetric ratio */
  bool qudaSmearingP;
  

  // !!! THIS CODE NOW ONLY WORKS FOR MOMENTUM-INJECTED SMEARING !!! //
  multi1d<Real> mom; /*<! quasi-momentum */
};

//! Gaussian quark smearing
/*! @ingroup smear
 *
 * Gaussian quark smearing object
 */
 
 
template<typename T>
class QuarkSmear: public QuarkSmearing<T> {
public:

  //! Full constructor
  QuarkSmear(const Params& p) :
        params(p) {
  }

  void gausSmear_plane(const multi1d<LatticeColorMatrix>& u, T& chi, Real width,
      int ItrGaus, int j_decay, int j_smear, int tr_size, double xi) const{
    T psi;
  
    Real ftmp = -(width * width) / Real(4 * ItrGaus);
    /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
    /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */
    Real ftmpi = Real(1) / ftmp;
    Real ftmp2 = Real(2 * Nd) + ftmpi;
  
    if (j_decay < Nd) ftmp2 = Real(2 * Nd - 2 ) + ftmpi;
    if (j_decay < Nd-1 && j_decay>=0) ftmp2 = ftmp2 +  2*xi*xi-2;
    if(j_smear>=0) ftmp2=2+ftmpi;

    multi1d<Real> mom(Nd);
    for(int mu=0;mu<Nd;mu++) mom[mu]=0.0;
  
    for (int n = 0; n < ItrGaus; ++n) {
      psi = chi * ftmp;
      chi = psi * ftmp2;
      for (int mu = 0; mu < Nd; ++mu)
        if ((mu != j_decay&&j_smear<0)||mu==j_smear) {
          double fac=(j_decay < Nd-1 && j_decay>=0 && mu==Nd-1)?xi*xi:1;
          chi -= fac*(u[mu] * shift(psi, FORWARD, mu)
                 + shift(adj(u[mu]) * psi, BACKWARD, mu));
        }
    }
    
    if(tr_size>1)
    {
       T vp1,vm1;
       for (int mu =0; mu<Nd; mu++)
       if(mu!=j_smear && mu!=j_decay)
       {
           vp1=chi;vm1=chi;
           for (int n = 1; n < tr_size; ++n)
           {
               psi=shift(vp1, FORWARD, mu);vp1=psi;
               psi=shift(vm1, BACKWARD, mu);vm1=psi;
               chi+=vp1;
               chi+=vm1;
           }
       }
    }
      
  }    

  //! Smear the quark
//  void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const ;
   void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const{
      StopWatch Timer;
      Timer.reset();
      Timer.start();
      multi1d<LatticeColorMatrix> links_single(Nd);
      const Real twopi = 6.283185307179586476925286;
      if(params.wvfIntPar>0)
      for(int mu=0; mu < Nd; mu++) 
      {
          Real k = params.mom[mu]*twopi/Layout::lattSize()[mu];
          QDPIO::cout << mu << ":" << cos(k) << "+" << -sin(k) << "I;"  << std::endl;
          links_single[mu] = u[mu]*cmplx(cos(k),-sin(k));
      }
        
      LatticeComplex corr=trace(adj(quark)*quark);
      double value=QDP::toDouble(sum(corr).elem().elem().elem().real());
      QDPIO::cout << "Norm2 is (before smearing)" << value << std::endl;
      Timer.stop();
      QDPIO::cout << "time to prepare = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

      if(params.wvfIntPar!=0)
          gausSmear_plane(links_single, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir,params.smear_dir, params.smear_tr_size, params.xi);
      corr=trace(adj(quark)*quark);
      value=QDP::toDouble(sum(corr).elem().elem().elem().real());
      QDPIO::cout << "Norm2 is (after smearing)" << value << std::endl;
      fflush(stdout);
   }
  
private:
  //! Hide partial constructor
  QuarkSmear() {
  }

private:
  Params params; /*!< smearing params */
};

}  // end namespace

//! Reader
/*! @ingroup smear */
void read(XMLReader& xml, const std::string& path,
    MomQuarkSmearingEnv::Params& param);

//! Writer
/*! @ingroup smear */
void write(XMLWriter& xml, const std::string& path,
    const MomQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
