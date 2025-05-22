// $Id: gaus_quark_smearing.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color std::vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "mom_quark_smearing.h"
#include "meas/smear/gaus_smear.h"
#include "util/ferm/transf.h"

namespace Chroma {

// Read parameters
void read(XMLReader& xml, const std::string& path,
    MomQuarkSmearingEnv::Params& param) {
  MomQuarkSmearingEnv::Params tmp(xml, path);
  param = tmp;
}

//! Parameters for running code
void write(XMLWriter& xml, const std::string& path,
    const MomQuarkSmearingEnv::Params& param) {
  param.writeXML(xml, path);
}

//! Hooks to register the class
namespace MomQuarkSmearingEnv {
namespace {
//! Callback function
QuarkSmearing<LatticePropagator>* createProp(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticePropagator>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeStaggeredPropagator>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeFermion>* createFerm(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeFermion>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeColorVector>* createColorVec(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeColorVector>(Params(xml_in, path));
}

//! Local registration flag
bool registered = false;

//! Name to be used
const std::string name = "MOM_GAUSSIAN";
}

//! Return the name
std::string getName() {
  return name;
}

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= Chroma::ThePropSmearingFactory::Instance().registerObject(name,
      createProp);
    success &= Chroma::TheStagPropSmearingFactory::Instance().registerObject(
      name, createStagProp);
    success &= Chroma::TheFermSmearingFactory::Instance().registerObject(name,
      createFerm);
    success &= Chroma::TheColorVecSmearingFactory::Instance().registerObject(
      name, createColorVec);
    registered = true;
  }
  return success;
}

//! Parameters for running code
Params::Params(XMLReader& xml, const std::string& path) {
  XMLReader paramtop(xml, path);

  read(paramtop, "wvf_param", wvf_param);
  read(paramtop, "wvfIntPar", wvfIntPar);
  read(paramtop, "no_smear_dir", no_smear_dir);
  if( paramtop.count("xi") > 0 ){
      read(paramtop, "xi", xi);
  }
  else
      xi=1.0;
  if( paramtop.count("smear_dir") > 0 ){
      read(paramtop, "smear_dir", smear_dir);
      if( paramtop.count("smear_tr_size") > 0 )
          read(paramtop, "smear_tr_size", smear_tr_size);
      else
          smear_tr_size=1;
  }
  else
  {
      smear_tr_size=0;
      smear_dir=-1;
  }

  read(paramtop, "mom", mom);
  if (mom.size() != Nd) {
    QDPIO::cerr << name << ": wrong size of mom array: expected length=" << Nd
    << std::endl;
    QDP_abort(1);
  }
}

//! Parameters for running code
void Params::writeXML(XMLWriter& xml, const std::string& path) const {
  push(xml, path);

  write(xml, "wvf_kind", MomQuarkSmearingEnv::getName());
  write(xml, "wvf_param", wvf_param);
  write(xml, "wvfIntPar", wvfIntPar);
  write(xml, "no_smear_dir", no_smear_dir);
  write(xml, "xi", xi);
  write(xml, "mom", mom);

  pop(xml);
}

}  // end namespace
}  // end namespace Chroma

