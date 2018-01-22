/*
 * STRUMPACK -- STRUctured Matrices PACKage, Copyright (c) 2014, The Regents of
 * the University of California, through Lawrence Berkeley National Laboratory
 * (subject to receipt of any required approvals from the U.S. Dept. of Energy).
 * All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Technology Transfer Department at TTD@lbl.gov.
 *
 * NOTICE. This software is owned by the U.S. Department of Energy. As such, the
 * U.S. Government has been granted for itself and others acting on its behalf a
 * paid-up, nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and display publicly.
 * Beginning five (5) years after the date permission to assert copyright is
 * obtained from the U.S. Department of Energy, and subject to any subsequent five
 * (5) year renewals, the U.S. Government is granted for itself and others acting
 * on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
 * Software to reproduce, prepare derivative works, distribute copies to the
 * public, perform publicly and display publicly, and to permit others to do so.
 *
 * Developers: Pieter Ghysels, Francois-Henry Rouet, Xiaoye S. Li.
 *             (Lawrence Berkeley National Lab, Computational Research Division).
 *
 */
#include <iostream>
#include "StrumpackSparseSolver.hpp"
#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;
using namespace strumpack;
typedef double numType;
typedef int intType;

namespace sparse_solver {

  //template<typename numType, typename intType> 
  struct my_solver{
    my_solver(intType n)//int argc, char* argv[])
    {
      StrumpackSparseSolver<numType,intType> spss;
      spss.options().set_mc64job(0);
      spss.options().set_reordering_method(ReorderingStrategy::GEOMETRIC);
      //spss.options().set_from_command_line(argc, argv);

      intType N = n * n;
      intType nnz = 5 * N - 4 * n;
      CSRMatrix<numType,intType> A(N, nnz);
      intType* col_ptr = A.get_ptr();
      intType* row_ind = A.get_ind();
      numType* val = A.get_val();

      nnz = 0;
      col_ptr[0] = 0;
      for (intType row=0; row<n; row++) {
        for (intType col=0; col<n; col++) {
          intType ind = col+n*row;
          val[nnz] = 4.0;
          row_ind[nnz++] = ind;
          if (col > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-1; } // left
          if (col < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+1; } // right
          if (row > 0)  { val[nnz] = -1.0; row_ind[nnz++] = ind-n; } // up
          if (row < n-1){ val[nnz] = -1.0; row_ind[nnz++] = ind+n; } // down
          col_ptr[ind+1] = nnz;
        }
      }
      //A.set_symmetric_sparsity();

      std::vector<numType> b(N, numType(1.)), x(N, numType(0.));

      spss.set_csr_matrix(N, col_ptr, row_ind, val, true);
      spss.reorder(n, n);
      // spss.factor();   // not really necessary, called if needed by solve
      spss.solve(b.data(), x.data());

      // just a check, system is already solved, so solving again
      // with the solution as initial guess should stop immediately
      //spss.solve(b.data(), x.data(), true);
      b_res.swap(b);
      x_res.swap(x);

      std::cout << "# COMPONENTWISE SCALED RESIDUAL = " << A.max_scaled_residual(x_res.data(), b_res.data()) << std::endl;
    }
    std::vector<numType> b_res;
    std::vector<numType> x_res;
  };

  void export_strumpack_solver()
  {
    typedef return_value_policy<return_by_value> rbv;
    //templated constructor for int and size_t flex arrays of id
    class_<sparse_solver::my_solver>("sparse_solver", init<int>())
    .add_property("x",make_getter(&sparse_solver::my_solver::x_res, rbv()))
    .add_property("b",make_getter(&sparse_solver::my_solver::b_res, rbv()))
    ;
  }

  BOOST_PYTHON_MODULE(sparse_solver_ext)
  {
    //def("sparse_solver", sparse_solver, args("n"), "My docstring");
    export_strumpack_solver();
  }
}//sparse_solver
