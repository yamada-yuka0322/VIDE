#include "config.hpp"
#include "fortran.hpp"

static void customCosmotoolHandler()
{
  try {
    if (PyErr_Occurred())
      ;
    else
      throw;
  } catch (const CosmoTool::InvalidArgumentException& e) {
    PyErr_SetString(PyExc_ValueError, "Invalid argument");
  } catch (const CosmoTool::NoSuchFileException& e) {
    PyErr_SetString(PyExc_IOError, "No such file");
  } catch (const CosmoTool::InvalidUnformattedAccess& e) {
    PyErr_SetString(PyExc_RuntimeError, "Invalid fortran unformatted access");
  } catch (const CosmoTool::Exception& e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
  } catch (const std::bad_alloc& exn) {
    PyErr_SetString(PyExc_MemoryError, exn.what());
  } catch (const std::bad_cast& exn) {
    PyErr_SetString(PyExc_TypeError, exn.what());
  } catch (const std::domain_error& exn) {
    PyErr_SetString(PyExc_ValueError, exn.what());
  } catch (const std::invalid_argument& exn) {
    PyErr_SetString(PyExc_ValueError, exn.what());
  } catch (const std::ios_base::failure& exn) {
    // Unfortunately, in standard C++ we have no way of distinguishing EOF
    // from other errors here; be careful with the exception mask
    PyErr_SetString(PyExc_IOError, exn.what());
  } catch (const std::out_of_range& exn) {
    // Change out_of_range to IndexError
    PyErr_SetString(PyExc_IndexError, exn.what());
  } catch (const std::overflow_error& exn) {
    PyErr_SetString(PyExc_OverflowError, exn.what());
  } catch (const std::range_error& exn) {
    PyErr_SetString(PyExc_ArithmeticError, exn.what());
  } catch (const std::underflow_error& exn) {
    PyErr_SetString(PyExc_ArithmeticError, exn.what());
  } catch (const std::exception& exn) {
    PyErr_SetString(PyExc_RuntimeError, exn.what());
  } catch(...) {
    PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
  }
}
