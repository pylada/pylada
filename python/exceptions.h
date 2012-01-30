#ifndef LADA_PYTHON_EXCEPTIONS_H
#define LADA_PYTHON_EXCEPTIONS_H
#include "LaDaConfig.h"

#include <boost/python/object.hpp>
#include <boost/python/handle.hpp>
#include <boost/python/str.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include <root_exceptions.h>

namespace LaDa
{
  namespace error
  {
    //! Attribute error thrown explicitely by lada.
    struct AttributeError: virtual input, virtual pyerror {};
    //! Key error thrown explicitely by lada.
    struct KeyError: virtual input, virtual out_of_range, virtual pyerror {};
    //! Value error thrown explicitely by lada.
    struct ValueError: virtual input, virtual pyerror {};
    //! Index error thrown explicitely by lada.
    struct IndexError: virtual input, virtual pyerror {};
    //! Argument error thrown explicitely by lada.
    struct TypeError: virtual input, virtual pyerror {};
    //! Not implemented error thrown explicitely by lada.
    struct NotImplementedError: virtual internal, virtual pyerror {};
    //! Subclasses python's ImportError.
    struct ImportError: virtual internal, virtual pyerror {};
    //! Subclasses python's ImportError.
    struct InternalError: virtual internal, virtual pyerror {};
  }

  namespace python
  {
    // Class to declare and register c++ to python exceptions.
    template<class T> 
      class PyException
      {
        public:
          PyException() {}
          boost::python::object const &initialize( std::string const &_name, 
                                                   std::string const &_doc, 
                 boost::python::tuple const &_bases
                   = boost::python::make_tuple(boost::python::object(
                       boost::python::handle<>(
                         PyExc_StandardError))) ) 
          {
            static bool is_first = true;
            if(is_first)
            {
              namespace bp = boost::python;
              bp::dict d;
              d["__doc__"] = bp::str(_doc);
              d["name"] = bp::str(_name);
              name_ = _name;
              doc_ = _doc;
              exception_ = bp::object(bp::handle<>(bp::borrowed(
                PyErr_NewExceptionWithDoc(&name_[0], &doc_[0], _bases.ptr(), d.ptr()) ) ) );
              is_first = false;
            }
            return exception_;
          }

          void operator()(T const &_e) const
          {
            std::string message = boost::diagnostic_information(_e);
            if(    name_[0] == 'a' or name_[0] == 'e' or name_[0] == 'i'
                or name_[0] == 'o' or name_[0] == 'u' or name_[0] == 'y' )
              message += "Encountered an " + name_ + " error.";
            else 
              message += "Encountered a " + name_ + " error.";
            PyErr_SetString(exception_.ptr(), message.c_str());
          }

          static void throw_error(std::string const &_message)
          {
            PyErr_SetString(exception_.ptr(), _message.c_str());
            boost::python::throw_error_already_set();
          };

          static boost::python::object const & exception() { return exception_; }

        private:
          std::string name_;
          std::string doc_;
          //! static exception object. Act as global.
          static boost::python::object exception_;
      };
    template<class T> boost::python::object PyException<T>::exception_ = boost::python::object();
#   ifdef LADA_PYERROR
#     error LADA_PYERROR already  defined. 
#   endif
#   ifdef LADA_PYERROR_FORMAT
#     error LADA_PYERROR already  defined. 
#   endif
#   ifdef LADA_PYTHROW
#     error LADA_PYERROR already  defined. 
#   endif
    //! \def LADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with the interpreter, but no c++ exception.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define LADA_PYERROR(EXCEPTION, MESSAGE) \
      PyErr_SetString( ::LaDa::python::PyException< ::LaDa::error::EXCEPTION >::exception().ptr(), \
                       MESSAGE )
    //! \def LADA_PYERROR(EXCEPTION, MESSAGE)
    //!      Raises a python exception with a formatted message, but no c++ exception.
    //!      For formatting, see PyErr_Format from the python C API.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
#   define LADA_PYERROR_FORMAT(EXCEPTION, MESSAGE, OTHER) \
      PyErr_Format( ::LaDa::python::PyException< ::LaDa::error::EXCEPTION >::exception().ptr(), \
                    MESSAGE, OTHER )
    //! \def LADA_PYTHROW(EXCEPTION, MESSAGE)
    //!      Raises a boost exception where EXCEPTION is stored as pyexcetp and MESSAGE as string.
    //!      EXCEPTION should be an unqualified declared in python/exceptions.h.
    //!      This macro makes it easy to catch all thrown python exceptions in a single statement.
#   define LADA_PYTHROW(EXCEPTION, MESSAGE)                                                   \
    {                                                                                         \
      PyObject * const exception                                                              \
         = ::LaDa::python::PyException< ::LaDa::error::EXCEPTION >::exception().ptr();        \
      PyErr_SetString(exception, MESSAGE);                                                    \
      BOOST_THROW_EXCEPTION( ::LaDa::error::EXCEPTION()                                       \
                             << ::LaDa::error::pyexcept(exception)                            \
                             << ::LaDa::error::string(MESSAGE) );                             \
    }
  }
}
# endif 