#pragma once

#include "ArbitraryQuadrature.h"
#include "MooseTypes.h"

#include "libmesh/fe_base.h"

/**
 * Helper class that recomputes a field variable value at a physical point on
 * the fly.
 *
 * The intention is to be able to quickly reinit only a single variable at a
 * single point without modifying Assembly.
 */
class OnTheFlyVariableValue
{
public:
  /**
   * Constructor.
   * @param dim Mesh dimension
   * @param fe_type The FEType of the variable
   * @param dof_values The variable dof values
   * @param current_elem The current element (should probably be from Assembly)
   */
  OnTheFlyVariableValue(const unsigned int dim,
                        const libMesh::FEType fe_type,
                        const VariableValue & dof_values,
                        const Elem * const & current_elem);

  /**
   * Compute the field value in the current element at a point
   * @param point The point
   * @return The value
   */
  Real operator()(const Point & point);

protected:
  /// Mesh dimension
  const unsigned int _dim;
  /// FE for reinit
  std::unique_ptr<FEBase> _fe;
  /// Arbitrary quadrature for FE
  ArbitraryQuadrature _qrule;
  /// Current dof values
  const VariableValue & _dof_values;
  /// Current element
  const Elem * const & _current_elem;
};
