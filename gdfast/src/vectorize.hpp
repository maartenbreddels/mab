#pragma once 

#include <boost/python/make_function.hpp>
#include <boost/python/def_visitor.hpp>
#include <boost/python/signature.hpp>
#include <boost/mpl/at.hpp>
#include <pyublas/numpy.hpp>

namespace gd { 


template <class F, class R>
struct vectorizer
{
    vectorizer(F fn) : fn(fn) {}

    template <class A0, class A1>
    pyublas::numpy_vector<A1> operator()(A0& a0, pyublas::numpy_vector<A1>& a1)
    {
		pyublas::numpy_vector<A1> outputvec(a1.size());
		//printf("size = %d\n", a1.size());
		A1* input = a1.data().begin();
		A1* input_end = a1.data().end();
		A1* output = outputvec.data().begin();
		while(input != input_end) {
			//printf("%f %f\n", *input, (a0.*fn)(*input));
			*output++ = (a0.*fn)(*input++);
		}
        return outputvec;
		//return 0;
    }
    F fn;
};

template <class F>
struct visitor : boost::python::def_visitor<visitor<F> >
{
    visitor(F fn)
      : fn(fn)
    {}

    template <class Class, class Options, class Signature>
    void visit_aux(
        Class& cl, char const* name
      , Options const& options, Signature const&) const
    {
	    typedef typename boost::mpl::at_c<Signature,0>::type	return_type;
	    typedef typename boost::mpl::at_c<Signature,2>::type	call_type;
	    typedef typename boost::mpl::at_c<Signature, 1>::type	call_class;
	    typedef pyublas::numpy_vector<return_type>				return_type_vec;
	    typedef pyublas::numpy_vector<call_type>				call_type_vec;
//float t = signature;
		//const boost::mpl::vector3<return_type_vec, gd::JeansAnisotropicConstant&, return_type_vec> sig_vec;
	    //boost::mpl::vector3<double, cls, double> sig_veqqc();
		//float t = signature;
	    //boost::mpl::vector3<return_type_vec, gd::JeansAnisotropicConstant&, return_type_vec> sig_vec;
	    //boost::mpl::vector3<double, gd::JeansAnisotropicConstant&, double> sig_test;
	    boost::mpl::vector3<return_type_vec, call_class, call_type_vec> sig_vec;
	    boost::mpl::vector3<return_type, call_class, call_type> sig_test;
	    
        cl.def(
            name
          , boost::python::make_function(
                fn
              , options.policies()
              , options.keywords()
              , sig_test
            )
        );
		cl.def(
            name
          , boost::python::make_function(
                vectorizer<F, return_type>(fn)
              , options.policies()
              , options.keywords()
              , sig_vec
            )
        );
    }

    template <class Class, class Options>
    void visit(Class& cl, char const* name, Options const& options) const
    {
        this->visit_aux(
            cl, name, options
          , boost::python::detail::get_signature(fn, (typename Class::wrapped_type*)0)
        );
    }

    F fn;
};

// Member function adaptor that releases and aqcuires the GIL
// around the function call.
template <class F>
visitor<F> vectorize(F fn)
{
    return visitor<F>(fn);
}

}






