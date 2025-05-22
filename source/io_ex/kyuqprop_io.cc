// $Id: kyuqprop_io.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $
/*!
 * @file
 * @brief  Read/write a Kentucky quark propagator
 */

#include "chromabase.h"
#include "io/kyuqprop_io.h"
#include "util/ferm/transf.h"
#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"

namespace Chroma {

//! Read a Kentucky quark propagator
/*!
 * \ingroup io
 *
 * \param q          propagator ( Modify )
 * \param file       path ( Read )
 */    

SpinMatrixD KYToDRMat()
{

  SpinMatrixD U = zero;
  RealD     foo = RealD(1) / sqrt(RealD(2));
  ComplexD  one = cmplx( foo,RealD(0));
  ComplexD mone = cmplx(-foo,RealD(0));

  pokeSpin(U,  one, 0, 1);
  pokeSpin(U, mone, 0, 3);
  pokeSpin(U, mone, 1, 0);
  pokeSpin(U,  one, 1, 2);
  pokeSpin(U,  one, 2, 1);
  pokeSpin(U,  one, 2, 3);
  pokeSpin(U, mone, 3, 0);
  pokeSpin(U, mone, 3, 2);


  return U;
}

void writeKYUQprop2(LatticePropagator& q, const std::string& file)
{
  START_CODE();
  
  BinaryFileWriter bin(file);
/*  
  for(int site=0; site < Layout::vol(); ++site)
  {
     multi1d<int> coord = crtesn(site, Layout::lattSize());
     write(bin, q, coord);
  }
 */ 
   SpinMatrixD U = KYToDRMat();
   LatticePropagatorD q_w;
   LatticeFermionD f;
   LatticeColorVectorD  cv;
   multi2d<LatticeComplexD> vec(Nc,Ns);
   LatticeRealD tmp;
   
   q_w = adj(U) * q * U;
  
  for(int src_spin=0; src_spin < Ns; ++src_spin)
  for(int src_color=0; src_color < Nc; ++src_color)
  {
      PropToFerm(q_w,f,src_color, src_spin);
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
      {
         cv=peekSpin(f,snk_spin);
         for(int snk_color=0; snk_color < Nc; ++snk_color)
            vec(snk_color,snk_spin)=peekColor(cv,snk_color);
      }
      
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
      for(int snk_color=0; snk_color < Nc; ++snk_color)
      {
          tmp=real(vec(snk_color,snk_spin));
          write(bin,tmp);
      }

      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
      for(int snk_color=0; snk_color < Nc; ++snk_color)
      {
          tmp=imag(vec(snk_color,snk_spin));
          write(bin,tmp);
      }
       
  }
  bin.close();
  
  END_CODE();

}


void readKYUQprop2(LatticePropagator& q, const std::string& file)
{
  START_CODE();

  if (Nc != 3)
  {
    QDPIO::cerr << "readKYUQprop - only supports Nc=3" << std::endl;
    QDP_abort(1);
  }

  if (Ns != 4)
  {
    QDPIO::cerr << "readKYUQprop - only supports Ns=4" << std::endl;
    QDP_abort(1);
  }
  
  QDPIO::cout << "file_name = " << file << std::endl;

  BinaryFileReader bin(file);

  /* KY Indices: 
     x,y,z,t,snk_col,snk_spin,ri,src_col,src_spin 
     x is fastest (Fortran Order)
  */
  // KYU always uses 64 bits
  LatticePropagatorD q_old;
  multi2d<LatticeFermionD> source(Nc,Ns);
  multi2d<LatticeRealD> re(Nc,Ns);
  multi2d<LatticeRealD> im(Nc,Ns);
  LatticeFermionD f;
  LatticeColorVectorD  cv;

  for(int src_spin=0; src_spin < Ns; ++src_spin)
    for(int src_color=0; src_color < Nc; ++src_color)
    {
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  read(bin, re(snk_color,snk_spin));
	}

      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  read(bin, im(snk_color,snk_spin));
	}

      // Stuff into a fermion
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
      {
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  pokeColor(cv, 
		    cmplx(re(snk_color,snk_spin), im(snk_color,snk_spin)),
		    snk_color);
	}

	pokeSpin(f, cv, snk_spin);
      }

      // Hold temporarily in a multi2d - will need to rearrange src_spin later
      source(src_color,src_spin) = f;
    }

  bin.close();

  for(int src_dirac=0; src_dirac < Nd; ++src_dirac)
  for(int src_color=0; src_color < Nc; ++src_color)
    FermToProp(LatticeFermionD(source(src_color,src_dirac)),
	       q_old, src_color, src_dirac);
  // Now that we have read the prop, need to change the spin basis
  SpinMatrixD U = KYToDRMat();

  // And finally...
   q = U * q_old * adj(U);
  

#if 1
  {
    // HACK - compute rho correlator like thingy
    ComplexD  i   = cmplx(RealD(0), RealD(1));
    ComplexD mi   = -i;
    ComplexD mone = cmplx(-RealD(1), RealD(0));

    SpinMatrixD gamma_1 = zero;
    pokeSpin(gamma_1, mi, 0, 3);
    pokeSpin(gamma_1, mi, 1, 2);
    pokeSpin(gamma_1,  i, 2, 1);
    pokeSpin(gamma_1,  i, 3, 0);

    SpinMatrixD gamma_5 = zero;
    pokeSpin(gamma_5, mone, 0, 2);
    pokeSpin(gamma_5, mone, 1, 3);
    pokeSpin(gamma_5, mone, 2, 0);
    pokeSpin(gamma_5, mone, 3, 1);

 
    LatticeComplex ps_rho = trace(adj(gamma_5 * q_old * gamma_5) * gamma_1 * q_old * gamma_1);

    LatticeComplex dr_rho = trace(adj(Gamma(15) * q * Gamma(15)) * Gamma(1) * q * Gamma(1));

    QDPIO::cout << "readKYUQprop: norm2(ps_rho) = " << RealD(norm2(ps_rho)) << std::endl;
    QDPIO::cout << "readKYUQprop: norm2(dr_rho) = " << RealD(norm2(dr_rho)) << std::endl;
    QDPIO::cout << "readKYUQprop: norm2(ps_rho - dr_rho) = " << RealD(norm2(ps_rho-dr_rho)) << std::endl;
    QDPIO::cout << "readKYUQprop: norm2(ps_rho + dr_rho) = " << RealD(norm2(ps_rho+dr_rho)) << std::endl;

#if 1
    {
      // Keep a copy of the phases with NO momenta
      int j_decay = Nd-1;
      SftMom phases(0, true, j_decay);

      multi1d<DComplex> ps_hsum = sumMulti(ps_rho, phases.getSet());
      multi1d<DComplex> dr_hsum = sumMulti(dr_rho, phases.getSet());

      // Initial group
      XMLFileWriter xml_out("kyuqprop.xml");

      push(xml_out, "kyuqprop");
      write(xml_out, "norm_diff", RealD(norm2(ps_rho-dr_rho)));
      write(xml_out, "ps_hsum", ps_hsum);
      write(xml_out, "dr_hsum", ps_hsum);
      pop(xml_out);
      
      xml_out.close();
    }
#endif

  }
#endif


  END_CODE();
}

}  // end namespace Chroma
