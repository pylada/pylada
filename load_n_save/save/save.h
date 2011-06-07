#ifndef LADA_LOADNSAVE_SAVE_SAVE_H
#define LADA_LOADNSAVE_SAVE_SAVE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../tree/tree.h"
#include "../xpr/section.h"
#include "../string_type.h"
#include "../access.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace save
    {
      //! Links text with operators for loading purposes.
      class Save 
      {
        public:
          //! Verbose or quiet output.
          bool verbose;
     
          //! Constructor.
          Save() : verbose(true) {}
          //! Copy Constructor.
          Save( Save const &_c ) : verbose(_c.verbose) {}

          //! Returns an in put/output tree.
          boost::shared_ptr<tree::Base> operator()( xpr::Section const& _sec,
                                                    version_type _version = 0u ) const;
     
          //! Loads an archive.
          virtual bool is_loading() const { return false; } 

        protected:
          //! Class for parsing a section.
          class Section;
      };

      class Save :: Section : public parser_base::Section
      {
        public:
          //! Constructor.
          Section   (boost::shared_ptr<tree::Section> const &_tree, version_type _version)
                  : tree_(_tree), version_(_version) {}
          //! Copy constructor.
          Section(Section & _c) : tree_(_c.tree_), version_(_c.version_) {}

          bool operator()(xpr::Section const& _sec) const
            { return operator&(_sec); }
          virtual bool operator&( xpr::Section const& _sec ) const
            { return _sec.parse(*this, version_); }

          //! Saves an archive.
          virtual bool is_loading() const { return false; } 

          //! Parses content.
          virtual bool content( xpr::Section const & _sec, t_String const &_name = "" ) const;
          //! Double dispatch.
          virtual bool regular( xpr::Section const &_sec, xpr::regular_data const&_data ) const;
          //! Returns the current tree;
          boost::shared_ptr<tree::Section> tree() const { return tree_; }
          //! Starts recurrence;
          virtual void start_recurrence() const {};
          //! Increments recurrence.
          virtual void step_recurrence() const {};
          //! Starts recurrence;
          virtual void stop_recurrence() const {};
          //! Whether actually parsing or simply looking at grammar.
          virtual bool grammar_only() const {return false;}

        protected:
          //! Parses content of an expression range;
          template<class T> bool content_(T const &_range) const;
          //! Reference to the tree being built.
          boost::shared_ptr<tree::Section> tree_;
          //! Version of this archive.
          version_type version_;
      };

    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa


#endif
