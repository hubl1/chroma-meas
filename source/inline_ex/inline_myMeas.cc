#include "inline_myMeas.h"

#include "../io_ex/io_general_class.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "meas/glue/mesplq.h"
#include "meas/hadron/BuildingBlocks_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/make_xml_file.h"
#include "util/ferm/transf.h"
#include "util/info/proginfo.h"

namespace Chroma
{

// Environment in which the measurement lives (holds params and such)
namespace InlineMyMeasIOGEnv
{
namespace
{
// Function to register with a factory
AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path)
{
    return new InlineMyMeasIOG(InlineMyMeasIOGParams(xml_in, path));
}

//! Local registration flag
bool registered = false;
}  // namespace

// The name of the measurement for the XML file
const std::string name = "My_Measurements";

//! Register all the factories
bool registerAll()
{
    bool success = true;
    if (!registered)
    {
        success &= CreateFermStateEnv::registerAll();
        success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
        registered = true;
    }
    return success;
}
}  // end namespace InlineMyMeasIOGEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

//! Reader for parameters
void read(XMLReader& xml, const std::string& path, InlineMyMeasIOGParams::Param_t& param)
{
    XMLReader paramtop(xml, path);
    read(paramtop, "cfg_serial", param.cfg_serial);
    read(paramtop, "t_src", param.t_src);
    read(paramtop, "hadrons", param.hadrons);
    read(paramtop, "l_prop", param.l_prop);
    read(paramtop, "s_prop", param.s_prop);
    read(paramtop, "c_prop", param.c_prop);
    read(paramtop, "file_name", param.file_name);
}

//! Writer for parameters
void write(XMLWriter& xml, const std::string& path, const InlineMyMeasIOGParams::Param_t& param)
{
    push(xml, path);

    write(xml, "cfg_serial", param.cfg_serial);
    write(xml, "t_src", param.t_src);
    write(xml, "hadrons", param.hadrons);
    write(xml, "l_prop", param.l_prop);
    write(xml, "s_prop", param.s_prop);
    write(xml, "c_prop", param.c_prop);
    write(xml, "file_name", param.file_name);

    pop(xml);
}

// Construct params from XML
InlineMyMeasIOGParams::InlineMyMeasIOGParams(XMLReader& xml_in, const std::string& path)
{
    try
    {
        XMLReader paramtop(xml_in, path);
        frequency = 1;
        read(paramtop, "Param", param);  // Read in the parameters
    }
    catch (const std::string& e)
    {
        QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
        QDP_abort(1);
    }
}

// Write out the parameters we constructed
void InlineMyMeasIOGParams::write(XMLWriter& xml_out, const std::string& path)
{
    push(xml_out, path);

    Chroma::write(xml_out, "Param", param);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
}

// Set up the XML and invoke func, which does the acual work
void InlineMyMeasIOG::operator()(unsigned long update_no, XMLWriter& xml_out)
{
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
        std::string xml_file = makeXMLFileName(params.xml_file, update_no);

        push(xml_out, "BuildingBlocks");
        write(xml_out, "update_no", update_no);
        write(xml_out, "xml_file", xml_file);
        pop(xml_out);

        XMLFileWriter xml(xml_file);
        func(update_no, xml);
    }
    else
    {
        func(update_no, xml_out);
    }
}

//------------------------------------------------------------------------------
// Real work done here
void InlineMyMeasIOG::func(unsigned long update_no, XMLWriter& xml_out)
{
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    //--------------------------------------------------------------------------
    // Start building the output XML
    push(xml_out, "BuildingBlocks");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " BuildingBlocks" << std::endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i = 1; i < Nd; ++i)
    {
        QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out);  // Print out basic program info
    push(xml_out, "Output_version");
    write(xml_out, "out_version", 2);
    pop(xml_out);

    //--------------------------------------------------------------------------
    // Grab propagator
    LatticePropagator L_prop, S_prop, C_prop;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << std::endl;

    try
    {
        // Grab the forward propagator
        L_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.param.l_prop);
        S_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.param.s_prop);
        C_prop = TheNamedObjMap::Instance().getData<LatticePropagator>(params.param.c_prop);
    }
    catch (const std::string& e)
    {
        QDPIO::cerr << InlineMyMeasIOGEnv::name << ": propagators: error message: " << e << std::endl;
        QDP_abort(1);
    }
    QDPIO::cout << "All propagators successfully parsed" << std::endl;

    //--------------------------------------------------------------------------
    // Grab the hadrons
    multi1d<std::string> hadron_list = params.param.hadrons;
    for (int i = 0; i < hadron_list.size(); i++)
    {
        QDPIO::cout << "Calculate hadrons include " << hadron_list[i] << std::endl;
    }

    QDPIO::cout << "Total hadron number: " << hadron_list.size() << std::endl;

    int operator_no = hadron_list.size();
    int j_decay     = 3;

    LatticeReal Phases = 1.;

    // Keep a copy of the phases with no momenta
    SftMom     phases_nomom(0, true, j_decay);
    const Set& timeslice = phases_nomom.getSet();
    int        tlen      = Layout::lattSize()[j_decay];

    general_data_base res(params.param.file_name.c_str());
    res.add_dimension(dim_conf, 1, &params.param.cfg_serial);
    res.add_dimension(dim_operator, operator_no);
    res.add_dimension(dim_t, tlen);
    res.add_dimension(dim_complex, 2);
    if (Layout::primaryNode())
        res.initialize();

    SpinMatrix g_one = 1.0, g_mone = -1.0;
    SpinMatrix Cg5 = Gamma(5) * g_one;    // C g_5 = C gamma_5 = gamma^1 gamma^3 = Gamma(5)
    SpinMatrix Cgx = Gamma(11) * g_one;   // Cg^x=g^2g^4g^1=Gamma(2+8+1)=Gamma(11)
    SpinMatrix Cgy = Gamma(8) * g_mone;   // Cg^y=g^2g^4g^2=-g^4=-Gamma(8)
    SpinMatrix Cgz = Gamma(14) * g_mone;  // Cg^z=g^2g^4g^3=-g^2g^3g^4=-Gamma(14)

    int                nt = Layout::lattSize()[3];
    int                t0 = params.param.t_src;
    SpinMatrix         prj_p(0.5 * (g_one + (g_one * Gamma(8))));
    SpinMatrix         prj_m(0.5 * (g_one - (g_one * Gamma(8))));
    LatticeSpinMatrix  T_unpol  = where(((Layout::latticeCoordinate(3)) - t0 + nt) % nt < nt / 2, prj_p, prj_m);
    LatticeSpinMatrix  N_unpol  = where(((Layout::latticeCoordinate(3)) - t0 + nt) % nt < nt / 2, prj_m, prj_p);
    QDP::LatticeDouble aaa      = Layout::latticeCoordinate(3);
    LatticeSpinMatrix  m_factor = where(QDP::fabs(aaa - t0) < nt / 2, g_one, g_mone);

    int offset = 0;
    for (int i = 0; i < operator_no; i++)
    {
        LatticeComplex corr = zero;

        if (hadron_list[i] == "PION")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(15) * L_prop * Gamma(15));
        }
        else if (hadron_list[i] == "A0")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * L_prop);
        }
        else if (hadron_list[i] == "ETA_S")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(15) * S_prop * Gamma(15));
        }
        else if (hadron_list[i] == "KAON")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(15) * L_prop * Gamma(15));
        }
        else if (hadron_list[i] == "PION_A4")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(15) * Gamma(8) * L_prop * Gamma(15));
        }
        else if (hadron_list[i] == "ETA_S_A4")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(15) * Gamma(8) * S_prop * Gamma(15));
        }
        else if (hadron_list[i] == "KAON_A4")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(15) * Gamma(8) * L_prop * Gamma(15));
        }
        else if (hadron_list[i] == "PROTON")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "test")
        {
            corr = LatticeComplex(
                                  trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
        }
        else if (hadron_list[i] == "NUCLEON_STAR")
        {
            corr = LatticeComplex(trace(m_factor * N_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))) +
                                  trace(m_factor * N_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * L_prop))));
            corr = -corr;
        }
        else if (hadron_list[i] == "SIGMA+-")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * L_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "SIGMA_UDS")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))) );
        }
        else if (hadron_list[i] == "SIGMA_USD")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))) );
        }
        else if (hadron_list[i] == "SIGMA_USD_UDS")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
        }
        else if (hadron_list[i] == "LAMBDA")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cg5, Cg5 * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
            corr -= LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "LAMBDA_STAR")
        {
            corr += LatticeComplex(4 * trace(m_factor * N_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * N_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * N_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * N_unpol * traceColor(S_prop * quarkContract13(L_prop * Cg5, Cg5 * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * N_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
            corr -= LatticeComplex(2 * trace(m_factor * N_unpol * traceColor(L_prop * quarkContract13(S_prop * Cg5, Cg5 * L_prop))));
            corr = -corr;
        }
        else if (hadron_list[i] == "XI")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
        }
        else if (hadron_list[i] == "LAMBDA_C")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cg5, Cg5 * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cg5, Cg5 * C_prop))));
            corr -= LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "SIGMA_C")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cg5, Cg5 * L_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "XI_C")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cg5, Cg5 * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cg5, Cg5 * L_prop))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cg5, Cg5 * C_prop))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cg5, Cg5 * C_prop))));
            corr -= LatticeComplex(1 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cg5, Cg5 * S_prop))));
            corr -= LatticeComplex(1 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "XI_C_s")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))) );
        }
        else if (hadron_list[i] == "XI_C_old")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))) );
            corr += LatticeComplex( trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cg5, Cg5 * S_prop))));
        }
        else if (hadron_list[i] == "XI_C_s2")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * C_prop)))) );
        }
        else if (hadron_list[i] == "XI_C_PRIME")
        {
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cg5, Cg5 * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cg5, Cg5 * L_prop))));
        }
        else if (hadron_list[i] == "XI_C_PRIME_old")
        {
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cg5, Cg5 * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cg5, Cg5 * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cg5, Cg5 * C_prop))));
        }
        else if (hadron_list[i] == "XI_CC")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cg5, Cg5 * C_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cg5, Cg5 * C_prop))));
        }
        else if (hadron_list[i] == "OMEGA_C")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cg5, Cg5 * S_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cg5, Cg5 * S_prop))));
        }
        else if (hadron_list[i] == "OMEGA_CC")
        {
            corr = LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cg5, Cg5 * C_prop)))) +
                                  trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cg5, Cg5 * C_prop))));
        }
        else if (hadron_list[i] == "SIGMA_STAR+_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgx, Cgx * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "SIGMA_STAR+_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgy, Cgy * L_prop))));
        }
        else if (hadron_list[i] == "SIGMA_STAR+_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgz, Cgz * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "SIGMA_STAR_C_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgx, Cgx * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "SIGMA_STAR_C_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgy, Cgy * L_prop))));
        }
        else if (hadron_list[i] == "SIGMA_STAR_C_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * L_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * L_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgz, Cgz * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR+_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgx, Cgx * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR+_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgy, Cgy * S_prop))));
        }
        else if (hadron_list[i] == "XI_STAR+_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgz, Cgz * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR_CC_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgx, Cgx * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR_CC_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgy, Cgy * C_prop))));
        }
        else if (hadron_list[i] == "XI_STAR_CC_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgz, Cgz * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_STAR_C_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgx, Cgx * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_STAR_C_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgy, Cgy * S_prop))));
        }
        else if (hadron_list[i] == "OMEGA_STAR_C_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * S_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * S_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgz, Cgz * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_STAR_CC_x")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgx, Cgx * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_STAR_CC_y")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgy, Cgy * C_prop))));
        }
        else if (hadron_list[i] == "OMEGA_STAR_CC_z")
        {
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * C_prop)))));
            corr += LatticeComplex(2 * trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * C_prop)))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(4 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgz, Cgz * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR_C_x")
        {
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgx, Cgx * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgx, Cgx * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgx, Cgx * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * L_prop)))));
            corr *= -1;
        }
        else if (hadron_list[i] == "XI_STAR_C_y")
        {
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgy, Cgy * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgy, Cgy * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgy, Cgy * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * L_prop)))));
        }
        else if (hadron_list[i] == "XI_STAR_C_z")
        {
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * S_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(S_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(L_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(S_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * C_prop)))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(C_prop * Cgz, Cgz * S_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(L_prop * Cgz, Cgz * C_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(C_prop * Cgz, Cgz * L_prop))));
            corr += LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * L_prop)))));
            corr *= -1;
        }
        else if (hadron_list[i] == "DELTA+_x")
        {
            corr = 6 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgx, Cgx * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "DELTA+_y")
        {
            corr = 6 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgy, Cgy * L_prop))));
        }
        else if (hadron_list[i] == "DELTA+_z")
        {
            corr = 6 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgz, Cgz * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "DELTA++_x")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgx, Cgx * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgx, Cgx * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "DELTA++_y")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgy, Cgy * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgy, Cgy * L_prop))));
        }
        else if (hadron_list[i] == "DELTA++_z")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(L_prop * traceSpin(quarkContract13(L_prop * Cgz, Cgz * L_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(L_prop * quarkContract13(L_prop * Cgz, Cgz * L_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_x")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(S_prop * Cgx, Cgx * S_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgx, Cgx * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_y")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(S_prop * Cgy, Cgy * S_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgy, Cgy * S_prop))));
        }
        else if (hadron_list[i] == "OMEGA_z")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(S_prop * traceSpin(quarkContract13(S_prop * Cgz, Cgz * S_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(S_prop * quarkContract13(S_prop * Cgz, Cgz * S_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_CCC_x")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(C_prop * Cgx, Cgx * C_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgx, Cgx * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_CCC_y")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(C_prop * Cgy, Cgy * C_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgy, Cgy * C_prop))));
        }
        else if (hadron_list[i] == "OMEGA_CCC_z")
        {
            corr = 2 * LatticeComplex(trace(m_factor * T_unpol * traceColor(C_prop * traceSpin(quarkContract13(C_prop * Cgz, Cgz * C_prop)))) +
                                      2 * trace(m_factor * T_unpol * traceColor(C_prop * quarkContract13(C_prop * Cgz, Cgz * C_prop))));
            corr *= -1;
        }
        else if (hadron_list[i] == "RHO_x")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(1) * L_prop * Gamma(1));
        }
        else if (hadron_list[i] == "RHO_y")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(2) * L_prop * Gamma(2));
        }
        else if (hadron_list[i] == "RHO_z")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(4) * L_prop * Gamma(4));
        }
        else if (hadron_list[i] == "K_STAR_x")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(1) * S_prop * Gamma(1));
        }
        else if (hadron_list[i] == "K_STAR_y")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(2) * S_prop * Gamma(2));
        }
        else if (hadron_list[i] == "K_STAR_z")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(4) * S_prop * Gamma(4));
        }
        else if (hadron_list[i] == "PHI_x")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(1) * S_prop * Gamma(1));
        }
        else if (hadron_list[i] == "PHI_y")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(2) * S_prop * Gamma(2));
        }
        else if (hadron_list[i] == "PHI_z")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(4) * S_prop * Gamma(4));
        }
        else if (hadron_list[i] == "JPSI_x")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(1) * C_prop * Gamma(1));
        }
        else if (hadron_list[i] == "JPSI_y")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(2) * C_prop * Gamma(2));
        }
        else if (hadron_list[i] == "JPSI_z")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(4) * C_prop * Gamma(4));
        }
        else if (hadron_list[i] == "CHI_C1_x")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(1) * Gamma(15) * C_prop * Gamma(1) * Gamma(15));
        }
        else if (hadron_list[i] == "CHI_C1_y")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(2) * Gamma(15) * C_prop * Gamma(2) * Gamma(15));
        }
        else if (hadron_list[i] == "CHI_C1_z")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(4) * Gamma(15) * C_prop * Gamma(4) * Gamma(15));
        }
        else if (hadron_list[i] == "H_C_xy")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(1) * Gamma(2) * C_prop * Gamma(1) * Gamma(2));
        }
        else if (hadron_list[i] == "H_C_yz")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(2) * Gamma(4) * C_prop * Gamma(2) * Gamma(4));
        }
        else if (hadron_list[i] == "H_C_zx")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(4) * Gamma(1) * C_prop * Gamma(4) * Gamma(1));
        }
        else if (hadron_list[i] == "D")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(15) * C_prop * Gamma(15));
        }
        else if (hadron_list[i] == "ETA_C")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * Gamma(15) * C_prop * Gamma(15));
        }
        else if (hadron_list[i] == "D_S")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(15) * C_prop * Gamma(15));
        }
        else if (hadron_list[i] == "D0_STAR")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * C_prop);
        }
        else if (hadron_list[i] == "DS0_STAR")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * C_prop);
        }
        else if (hadron_list[i] == "CHI_C0")
        {
            corr = trace(adj(Gamma(15) * C_prop * Gamma(15)) * C_prop);
        }
        else if (hadron_list[i] == "D_STAR_x")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(1) * C_prop * Gamma(1));
        }
        else if (hadron_list[i] == "D_STAR_y")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(2) * C_prop * Gamma(2));
        }
        else if (hadron_list[i] == "D_STAR_z")
        {
            corr = trace(adj(Gamma(15) * L_prop * Gamma(15)) * Gamma(4) * C_prop * Gamma(4));
        }
        else if (hadron_list[i] == "DS_STAR_x")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(1) * C_prop * Gamma(1));
        }
        else if (hadron_list[i] == "DS_STAR_y")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(2) * C_prop * Gamma(2));
        }
        else if (hadron_list[i] == "DS_STAR_z")
        {
            corr = trace(adj(Gamma(15) * S_prop * Gamma(15)) * Gamma(4) * C_prop * Gamma(4));
        }
        else if (hadron_list[i] == "OMEGA_QA_x")  // sss, interpolator C\gamma^i, if i=x,
                                                  // C\gamma^x=g^2g^4g^1=Gamma(11); if i=y,
                                                  // -Gamma(8); if i=z, -Gamma(14)
        {
            LatticePropagator sink_prop = S_prop;
            SpinMatrix        Cgi       = Cgx;

            LatticePropagator qc_type1 = quarkContract13(sink_prop * Cgi, Cgi * sink_prop);
            LatticePropagator qc_type2 = quarkContract13(sink_prop, Cgi * sink_prop * Cgi);
            LatticePropagator qc_type3 = quarkContract13(sink_prop, Cgi * sink_prop);

            corr = LatticeComplex(trace(m_factor * traceSpin(qc_type1) * sink_prop) + trace(m_factor * sink_prop * qc_type1) + trace(m_factor * traceSpin(qc_type2) * sink_prop) +
                                  trace(m_factor * sink_prop * Cgi * qc_type3) + trace(m_factor * sink_prop * transposeSpin(qc_type2)) + trace(m_factor * sink_prop * Cgi * transposeSpin(qc_type3)));
            corr *= -1;
        }
        else if (hadron_list[i] == "OMEGA_QA_y")  // sss, interpolator C\gamma^i, if i=x,
                                                  // C\gamma^x=g^2g^4g^1=Gamma(11); if i=y,
                                                  // -Gamma(8); if i=z, -Gamma(14)
        {
            LatticePropagator sink_prop = S_prop;
            SpinMatrix        Cgi       = Cgy;

            LatticePropagator qc_type1 = quarkContract13(sink_prop * Cgi, Cgi * sink_prop);
            LatticePropagator qc_type2 = quarkContract13(sink_prop, Cgi * sink_prop * Cgi);
            LatticePropagator qc_type3 = quarkContract13(sink_prop, Cgi * sink_prop);

            corr = LatticeComplex(trace(m_factor * traceSpin(qc_type1) * sink_prop) + trace(m_factor * sink_prop * qc_type1) + trace(m_factor * traceSpin(qc_type2) * sink_prop) +
                                  trace(m_factor * sink_prop * Cgi * qc_type3) + trace(m_factor * sink_prop * transposeSpin(qc_type2)) + trace(m_factor * sink_prop * Cgi * transposeSpin(qc_type3)));
        }
        else if (hadron_list[i] == "OMEGA_QA_z")  // sss, interpolator C\gamma^i, if i=x,
                                                  // C\gamma^x=g^2g^4g^1=Gamma(11); if i=y,
                                                  // -Gamma(8); if i=z, -Gamma(14)
        {
            LatticePropagator sink_prop = S_prop;
            SpinMatrix        Cgi       = Cgz;

            LatticePropagator qc_type1 = quarkContract13(sink_prop * Cgi, Cgi * sink_prop);
            LatticePropagator qc_type2 = quarkContract13(sink_prop, Cgi * sink_prop * Cgi);
            LatticePropagator qc_type3 = quarkContract13(sink_prop, Cgi * sink_prop);

            corr = LatticeComplex(trace(m_factor * traceSpin(qc_type1) * sink_prop) + trace(m_factor * sink_prop * qc_type1) + trace(m_factor * traceSpin(qc_type2) * sink_prop) +
                                  trace(m_factor * sink_prop * Cgi * qc_type3) + trace(m_factor * sink_prop * transposeSpin(qc_type2)) + trace(m_factor * sink_prop * Cgi * transposeSpin(qc_type3)));
            corr *= -1;
        }
        else
        {
            QDPIO::cerr << "Unknown hadron name: " << hadron_list[i] << std::endl;
            QDP_abort(1);
        }

        multi1d<DComplex> hsum = sumMulti(Phases * corr, timeslice);
        if (Layout::primaryNode()){
            QDPIO::cout << "===== hadron_list[" << i << "] = " << hadron_list[i] << " =====" << std::endl;
            for (int t = 0; t < tlen; ++t)
            {
                res.data[offset * tlen * 2 + 2 * t]     = QDP::toDouble(hsum[(t + params.param.t_src) % tlen].elem().elem().elem().real());
                res.data[offset * tlen * 2 + 2 * t + 1] = QDP::toDouble(hsum[(t + params.param.t_src) % tlen].elem().elem().elem().imag());
                     QDPIO::cout << "t = " << t
                    << "  data_re = "
                    << res.data[offset * tlen * 2 + 2 * t]
                    << "  data_im = "
                    << res.data[offset * tlen * 2 + 2 * t + 1]
                    << std::endl;
            }
        offset++;
    }
  }

    if (Layout::primaryNode())
        res.save();

    snoop.stop();
    QDPIO::cout << InlineMyMeasIOGEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;

    QDPIO::cout << InlineMyMeasIOGEnv::name << ": ran successfully" << std::endl;

    END_CODE();
}  // end of InlineMyMeasIOG::func

};  // end of namespace Chroma
