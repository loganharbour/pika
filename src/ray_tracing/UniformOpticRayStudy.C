#include "UniformOpticRayStudy.h"

#include "MooseUtils.h"

registerMooseObject("PikaApp", UniformOpticRayStudy);

InputParameters
UniformOpticRayStudy::validParams()
{
  auto params = OpticRayStudyBase::validParams();

  params.addRequiredParam<std::vector<BoundaryName>>("boundary", "The list of boundary IDs");
  params.addRequiredParam<Point>("direction", "g");

  params.addRequiredParam<std::vector<Real>>("intensity", "intensity");

  params.addRequiredParam<std::vector<std::vector<Real>>>("splits", "g");

  return params;
}

UniformOpticRayStudy::UniformOpticRayStudy(const InputParameters & parameters)
  : OpticRayStudyBase(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _direction(getParam<Point>("direction")),
    _intensity(getParam<std::vector<Real>>("intensity")),
    _splits(getParam<std::vector<std::vector<Real>>>("splits"))
{
  if (!isRectangularDomain())
    mooseError("Optic ray study ", name(), " does not support non-rectangular domains");
  if (_intensity.size() != _num_energy_groups)
    paramError("intesntiy",
               "Should be of the same size as the number of energy groups, ",
               _num_energy_groups);
  if (_splits.size() != _bnd_ids.size())
    paramError("splits", "Should be the same length as boundary");
  for (unsigned int i = 0; i < _splits.size(); ++i)
    if (_mesh.dimension() == 2 && _splits[i].size() != 1)
      paramError("splits",
                 "The split at index ",
                 i,
                 " has ",
                 _splits[i].size(),
                 " entries, but it should have only 1 in 2D");
    else if (_mesh.dimension() == 3 && _splits[i].size() != 2)
      paramError("splits",
                 "The split at index ",
                 i,
                 " has ",
                 _splits[i].size(),
                 " entries, but it should have 2 in 3D");
}

void
UniformOpticRayStudy::buildBoundingBoxes()
{
  _local_bboxes.resize(_bnd_ids.size(), BoundingBox());
  _bboxes.resize(_bnd_ids.size(), BoundingBox());
  _widths.resize(_bnd_ids.size());
  _dims.resize(_bnd_ids.size());

  // First build the bounding box for each of the boundaries we need
  for (const auto & bnd_elem : *_mesh.getBoundaryElementRange())
  {
    const Elem * elem = bnd_elem->_elem;
    if (elem->processor_id() != processor_id())
      continue;

    const auto side = bnd_elem->_side;
    const auto bnd_id = bnd_elem->_bnd_id;

    for (unsigned int i_bnd = 0; i_bnd < _bnd_ids.size(); ++i_bnd)
      if (bnd_id == _bnd_ids[i_bnd])
      {
        _local_bboxes[i_bnd].union_with(elem->side_ptr(side)->loose_bounding_box());
        break;
      }
  }

  for (unsigned int i_bnd = 0; i_bnd < _bnd_ids.size(); ++i_bnd)
  {
    auto & local_bbox = _local_bboxes[i_bnd];
    auto & bbox = _bboxes[i_bnd];
    std::vector<Real> & widths = _widths[i_bnd];
    std::vector<unsigned int> & dims = _dims[i_bnd];
    std::vector<Real> & splits = _splits[i_bnd];

    // Decide on the global bounding box for this boundary
    bbox = local_bbox;
    _communicator.min(bbox.min());
    _communicator.max(bbox.max());

    // Find the dimensions in the bounding box that aren't the same
    for (unsigned int dim = 0; dim < 3; ++dim)
      if (!MooseUtils::absoluteFuzzyEqual(bbox.first(dim), bbox.second(dim)))
      {
        dims.push_back(dim);
        widths.push_back(std::abs(bbox.first(dim) - bbox.second(dim)));
      }

    // Make sure we're planar!
    if ((_mesh.dimension() == 2 && dims.size() != 1) ||
        (_mesh.dimension() == 3 && dims.size() != 2))
      mooseError("In optic study '",
                 name(),
                 "', the requested boundary with ID ",
                 _bnd_ids[i_bnd],
                 " is not planar.\n\nRequested boundaries must be planar to use this object.");

    // If 2D, we're going to add some dummy dimensions so that we can use the same two for loops to
    // fill the Rays regardless of 2D or 3D
    if (_mesh.dimension() == 2)
    {
      dims.push_back(dims[0] == 0 ? 1 : 0);
      widths.push_back(0);
      splits.push_back(1);
    }
  }
}

void
UniformOpticRayStudy::defineRays()
{
  buildBoundingBoxes();

  unsigned int created = 0;

  // To be filled with the starting ID for the Rays for a particular boundary
  RayID start_id = 0;

  for (unsigned int i_bnd = 0; i_bnd < _bnd_ids.size(); ++i_bnd)
  {
    // The ID that we will actually increment as we generated rays
    RayID id = start_id;

    const auto & local_bbox = _local_bboxes[i_bnd];
    const auto & bbox = _bboxes[i_bnd];
    const std::vector<Real> & widths = _widths[i_bnd];
    const std::vector<unsigned int> & dims = _dims[i_bnd];
    const std::vector<Real> & splits = _splits[i_bnd];

    // Start with the lower bounds for the starting point, and change
    // individual dimensions as needed
    Point start_point(bbox.first);

    // The width associate with each region that a Ray is a part of
    const Real di = widths[0] / splits[0];
    const Real dj = widths[1] / splits[1];
    // The factor to multiply by the intensity (is either the length
    // or area associated with a ray)
    const Real energy_factor = di * (dj > 0 ? dj : 1);

    for (unsigned int i = 0; i < splits[0]; ++i)
    {
      start_point(dims[0]) = bbox.first(dims[0]) + (i + 0.5) * di;
      for (unsigned int j = 0; j < splits[1]; ++j)
      {
        start_point(dims[1]) = bbox.first(dims[1]) + (j + 0.5) * dj;

        // If our local box doesn't contain it, we absolutely don't have this ray
        if (!local_bbox.contains_point(start_point))
        {
          id++; // still increment the ID!
          continue;
        }

        ++created;

        // Define the ray and set the data!
        Ray & ray = defineRay(start_point, _direction, id++);
        for (unsigned int g = 0; g < _num_energy_groups; ++g)
          ray.data(_energy_data_index + g) = _intensity[g] * energy_factor;
        std::cerr << "set data " << ray.data(0) << std::endl;
      }
    }
  }
}
