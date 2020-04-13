#include "OnTheFlyVariableValue.h"

// MOOSE includes
#include "MooseArray.h"

// libMesh includes
#include "libmesh/elem.h"

OnTheFlyVariableValue::OnTheFlyVariableValue(const unsigned int dim,
                                             const libMesh::FEType fe_type,
                                             const VariableValue & dof_values,
                                             const Elem * const & current_elem)
  : _dim(dim),
    _fe(FEBase::build(_dim, fe_type)),
    _qrule(_dim),
    _dof_values(dof_values),
    _current_elem(current_elem)
{
  _fe->attach_quadrature_rule(&_qrule);
  _fe->get_phi();
}

Real
OnTheFlyVariableValue::operator()(const Point & point)
{
  mooseAssert(_current_elem->contains_point(point), "Doesn't contain point");

  const auto ref_point = FEMap::inverse_map(_dim, _current_elem, point);
  _qrule.setPoints({ref_point});

  _fe->reinit(_current_elem);
  const auto & phi = _fe->get_phi();
  mooseAssert(_dof_values.size() == phi.size(), "Size mismatch");

  Real value = 0;
  for (unsigned int i = 0; i < _fe->n_shape_functions(); ++i)
    value += phi[i][0] * _dof_values[i];

  return value;
}
