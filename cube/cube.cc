// Includes
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <deal.II/grid/grid_out.h>

using namespace dealii;
//////////////////////////////////////////////////////////////////////////////
class Poisson
{
public:
	Poisson();
	void run();
private:
	void make_grid();
	void setup_system();
	void assemble_system();
	void solve();
	void optput_results() const;
	
	Triangulation<2>      triangulation;
	FE_Q<2>				  fe;
	DoFHandler<2>		  dof_handler;
	
	SparsityPattern       sparsity_pattern;
	SparseMatrix<double>  system_matrix;
	
	Vector<double>        solution;
	Vector<double>        system_rhs;
};

Poisson::Poisson():fe(1), dof_handler(triangulation) {}
//////////////////////////////////////////////////////////////////////////////
void Poisson::make_grid()
{
	Point<2> center(0,0);
	GridGenerator::hyper_cube(triangulation, -1, 1);
	triangulation.refine_global(5);
	
	std::ofstream out("grid.svg");
    GridOut       grid_out;
    grid_out.write_svg(triangulation, out);
}
//////////////////////////////////////////////////////////////////////////////
void Poisson::setup_system()
{
	dof_handler.distribute_dofs(fe);
	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);
	
	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}
//////////////////////////////////////////////////////////////////////////////
void Poisson::assemble_system()
{
	QGauss<2> quadrature_formula(fe.degree+1);
	FEValues<2> fe_values(fe, quadrature_formula,
						 update_values | update_gradients | update_JxW_values);
						 
	const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
	
	const double F=1.0;
	
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);
	
	std::vector<types::global_dof_index>   local_dof_indices(dofs_per_cell);
	
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs    = 0;
		
		for (const unsigned int q_index : fe_values.quadrature_point_indices())
		{
			for (const unsigned int i : fe_values.dof_indices())
				for (const unsigned int j: fe_values.dof_indices())
					cell_matrix(i,j) +=
						(fe_values.shape_grad(i, q_index)*
						 fe_values.shape_grad(j, q_index)*
						 fe_values.JxW(q_index));
		
			for (const unsigned int i : fe_values.dof_indices())
				cell_rhs(i) +=(fe_values.shape_value(i, q_index) *
				F *
				fe_values.JxW(q_index));
		}
	cell -> get_dof_indices(local_dof_indices);
	
	for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

	for (const unsigned int i : fe_values.dof_indices())
		system_rhs(local_dof_indices[i]) += cell_rhs(i);
	} 
	
	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler,
											 0,
											 Functions::ZeroFunction<2>(),
											 boundary_values);
											 
	MatrixTools::apply_boundary_values(boundary_values,
									   system_matrix,
									   solution,
									   system_rhs);
}
//////////////////////////////////////////////////////////////////////////////
void Poisson::solve()
{
	SolverControl		       solver_control(1000, 1e-6 * system_rhs.l2_norm());
	SolverCG<Vector<double>>   solver(solver_control);
	
	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
	
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches();
    std::ofstream output1("solution.vtk");
    std::ofstream output2("solution.gpl");
    data_out.write_vtk(output1);
    data_out.write_gnuplot(output2);
}
//////////////////////////////////////////////////////////////////////////////
void Poisson::run()
{
	make_grid();
	setup_system();
	assemble_system();
	solve();
}
////////////////////////////////////////////////////////////////////////////////
int main()
{
	Poisson poisson_problem;
	poisson_problem.run();
}

