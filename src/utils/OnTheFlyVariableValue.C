#include "OnTheFlyVariableValue.h"

// MOOSE includes
#include "MooseArray.h"

// libMesh includes
#include "libmesh/elem.h"

OnTheFlyVariableValue::OnTheFlyVariableValue(const unsigned int dim,
                                             const libMesh::FEType fe_type,
                                             const VariableValue & dof_values,
                                             const Elem * const & current_elem,
                                             const std::set<VariableValueType> & value_types)
  : _dim(dim),
    _fe(FEBase::build(_dim, fe_type)),
    _qrule(_dim),
    _dof_values(dof_values),
    _current_elem(current_elem),
    _has_value(value_types.count(VALUE)),
    _has_grad_value(value_types.count(GRAD_VALUE))
{
  _fe->attach_quadrature_rule(&_qrule);
  if (value_types.count(VALUE))
    _fe->request_phi();

  if (_has_value)
    _fe->request_phi();
  if (_has_grad_value)
    _fe->request_dphi();
}

Real
OnTheFlyVariableValue::value(const Point & point)
{
  mooseAssert(_current_elem->contains_point(point), "Doesn't contain point");

  if (!_has_value)
    mooseError("OnTheFlyVariableValue does not have value requested");

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

RealGradient
OnTheFlyVariableValue::gradValue(const Point & point)
{
  mooseAssert(_current_elem->contains_point(point), "Doesn't contain point");
  if (!_has_grad_value)
    mooseError("OnTheFlyVariableValue does not have gradient requested");

  const auto ref_point = FEMap::inverse_map(_dim, _current_elem, point);
  _qrule.setPoints({ref_point});

  _fe->reinit(_current_elem);
  const auto & dphi = _fe->get_dphi();
  mooseAssert(_dof_values.size() == dphi.size(), "Size mismatch");

  RealGradient value = 0;
  for (unsigned int i = 0; i < _fe->n_shape_functions(); ++i)
    value.add_scaled(dphi[i][0], _dof_values[i]);

  return value;
}
